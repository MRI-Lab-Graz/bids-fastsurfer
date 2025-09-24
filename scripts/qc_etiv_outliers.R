#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

has_pkg <- function(p) requireNamespace(p, quietly=TRUE)

option_list <- list(
  make_option(c("-q", "--qdec"), type="character", help="Path to qdec.table.dat (TSV)"),
  make_option(c("-a", "--aseg"), type="character", help="Path to aseg.long.table (whitespace)"),
  make_option(c("-o", "--outdir"), type="character", default="qc", help="Output directory [default %default]"),
  make_option(c("--time-col"), type="character", default=NULL, help="Time column name (default: 'tp' if present)"),
  make_option(c("--id-col"), type="character", default=NULL, help="aseg ID column name (default auto: Measure.volume or Measure:volume)"),
  make_option(c("--etiv-col"), type="character", default=NULL, help="Column for TIV (default auto: EstimatedTotalIntraCranialVol/eTIV/to.eTIV)"),
  make_option(c("--roi-pattern"), type="character", default="Thalamus|Caudate|Putamen|Pallidum|Hippocampus|Amygdala|Accumbens\\.area|VentralDC", 
              help="Regex of ROI names to include (applied after global exclusions) [default %default]"),
  make_option(c("--exclude-pattern"), type="character", default="Ventricle|Choroid\\.plexus|vessel|WM\\.hypointensities|Cerebellum|Brain\\.Stem|Cortex", 
              help="Regex of ROI names to exclude [default %default]"),
  make_option(c("--no-labels"), action="store_true", default=FALSE, help="Do not label flagged points in plots"),
  make_option(c("--plot-top"), type="integer", default=12, help="Number of ROIs with most flags to plot [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$qdec) || is.null(opt$aseg)) stop("--qdec and --aseg are required")
if (!file.exists(opt$qdec)) stop(sprintf("qdec not found: %s", opt$qdec))
if (!file.exists(opt$aseg)) stop(sprintf("aseg not found: %s", opt$aseg))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
plot_dir <- file.path(opt$outdir, "plots"); dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

# Read tables
qdec <- tryCatch(read.delim(opt$qdec, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
if (is.null(qdec)) qdec <- read.table(opt$qdec, header=TRUE, stringsAsFactors=FALSE)
aseg <- read.table(opt$aseg, header=TRUE, stringsAsFactors=FALSE)

# Standardize qdec id column names
if ("fsid-base" %in% names(qdec)) names(qdec)[names(qdec)=="fsid-base"] <- "fsid_base"
if ("fsid.base" %in% names(qdec)) names(qdec)[names(qdec)=="fsid.base"] <- "fsid_base"

# Determine id column in aseg and derive fsid/fsid_base
id_col <- if (!is.null(opt$id_col) && opt$id_col %in% names(aseg)) opt$id_col else {
  if ("Measure.volume" %in% names(aseg)) "Measure.volume" else if ("Measure:volume" %in% names(aseg)) "Measure:volume" else names(aseg)[1]
}
ids <- as.character(aseg[[id_col]])
has_long <- grepl("\\.long\\.", ids, perl=TRUE)
aseg$fsid <- ifelse(has_long, sub("^(.*)\\.long\\..*$", "\\1", ids, perl=TRUE), ids)
aseg$fsid_base <- ifelse(has_long, sub("^.*\\.long\\.(.*)$", "\\1", ids, perl=TRUE), NA_character_)

# Time col
time_col <- if (!is.null(opt$`time-col`)) opt$`time-col` else if ("tp" %in% names(qdec)) "tp" else stop("time_col not specified and 'tp' not in qdec")

# Merge
if (!all(c("fsid","fsid_base") %in% names(qdec))) stop("qdec must have fsid and fsid_base")
if (!all(c("fsid","fsid_base") %in% names(aseg))) stop("aseg must include derived fsid/fsid_base")

# Coerce numeric-like qdec columns
coerce_numeric_like <- function(x) {
  if (!is.character(x)) return(x)
  stripped <- trimws(tolower(x))
  stripped[stripped %in% c("", "na", "n/a", "nan", "null")] <- NA_character_
  nums <- suppressWarnings(as.numeric(stripped))
  if (any(!is.na(stripped) & is.na(nums))) {
    x[is.na(stripped)] <- NA_character_
    return(x)
  }
  x[is.na(stripped)] <- NA
  nums
}
if (nrow(qdec) > 0) {
  for (nm in names(qdec)) qdec[[nm]] <- coerce_numeric_like(qdec[[nm]])
}

# Merge
dat <- merge(qdec, aseg, by=c("fsid","fsid_base"))
if (nrow(dat) == 0) stop("Merged data is empty")

# Determine eTIV column
etiv_candidates <- c(opt$`etiv-col`, "EstimatedTotalIntraCranialVol", "eTIV", "to.eTIV")
etiv_candidates <- unique(etiv_candidates[!is.na(etiv_candidates)])
etiv_col <- NULL
for (c in etiv_candidates) { if (!is.null(c) && c %in% names(dat)) { etiv_col <- c; break } }
if (is.null(etiv_col)) stop("Could not find eTIV column (try --etiv-col)")

# Baseline subset: min time per base
if (!time_col %in% names(dat)) stop(sprintf("time_col '%s' not found after merge", time_col))
baseline_time <- tapply(dat[[time_col]], dat$fsid_base, function(x) suppressWarnings(min(x, na.rm=TRUE)))
keep <- mapply(function(sb, tp) sb == dat$fsid_base & !is.na(dat$fsid_base) & dat[[time_col]] == tp, dat$fsid_base, baseline_time[dat$fsid_base])
# 'keep' is vectorized wrong, fix via indexing
keep <- rep(FALSE, nrow(dat))
for (sb in unique(dat$fsid_base)) {
  if (is.na(sb)) next
  tpmin <- baseline_time[[sb]]
  if (is.infinite(tpmin) || is.na(tpmin)) next
  idx <- which(dat$fsid_base == sb & dat[[time_col]] == tpmin)
  if (length(idx)) keep[idx[1]] <- TRUE
}
baseline <- dat[keep, , drop=FALSE]

# Identify ROI columns: everything not in keys/qdec/id/etiv and drop global measures
exclude <- unique(c("fsid","fsid_base", time_col, names(qdec), id_col, etiv_col))
brain_cols <- setdiff(names(baseline), exclude)
summary_patterns <- c(
  "^BrainSegVol$", "^BrainSegVolNotVent$", "^CortexVol$",
  "^lhCortexVol$", "^rhCortexVol$", "^CerebralWhiteMatterVol$",
  "^lhCerebralWhiteMatterVol$", "^rhCerebralWhiteMatterVol$",
  "^TotalGrayVol$", "^SubCortGrayVol$", "^SupraTentorialVol$",
  "^SupraTentorialVolNotVent$", "^MaskVol$", "eTIV", "to\\.eTIV$",
  "^EstimatedTotalIntraCranialVol$", "^WM\\.hypointensities$",
  "hypointensities$", "^CC_", "^Optic\\.Chiasm$", "^CSF$",
  "^X3rd\\.Ventricle$", "^X4th\\.Ventricle$", "^X5th\\.Ventricle$",
  "baseline_value"
)
keep_mask <- !Reduce(`|`, lapply(summary_patterns, function(p) grepl(p, brain_cols, perl=TRUE)))
brain_cols <- brain_cols[keep_mask]
# Apply user include/exclude patterns
if (!is.null(opt$`roi-pattern`) && nchar(opt$`roi-pattern`) > 0) {
  brain_cols <- brain_cols[grepl(opt$`roi-pattern`, brain_cols, perl=TRUE)]
}
if (!is.null(opt$`exclude-pattern`) && nchar(opt$`exclude-pattern`) > 0) {
  brain_cols <- brain_cols[!grepl(opt$`exclude-pattern`, brain_cols, perl=TRUE)]
}
if (!length(brain_cols)) stop("No ROI columns found for QC")

mad_z <- function(x) {
  x <- as.numeric(x)
  m <- median(x, na.rm=TRUE)
  s <- mad(x, center=m, constant=1.4826, na.rm=TRUE)
  if (isTRUE(is.na(s)) || s == 0) {
    mu <- mean(x, na.rm=TRUE); sdv <- sd(x, na.rm=TRUE)
    if (isTRUE(is.na(sdv)) || sdv == 0) return(rep(0, length(x)))
    return((x - mu) / sdv)
  }
  (x - m) / s
}

out_rows <- list()
counts <- list()

for (roi in brain_cols) {
  y <- suppressWarnings(as.numeric(baseline[[roi]]))
  etiv <- suppressWarnings(as.numeric(baseline[[etiv_col]]))
  ok <- is.finite(y) & is.finite(etiv)
  if (!any(ok)) next
  y_ok <- y[ok]; etiv_ok <- etiv[ok]
  if (length(unique(y_ok)) < 3) next
  fit <- try(lm(y_ok ~ etiv_ok), silent=TRUE)
  if (inherits(fit, "try-error")) next
  pred <- as.numeric(fitted(fit))
  res <- y_ok - pred
  z <- mad_z(res)
  flag <- abs(z) > 3
  if (any(flag)) {
    idxs <- which(ok)[flag]
    out_rows[[roi]] <- data.frame(
      roi=roi,
      fsid_base=baseline$fsid_base[idxs],
      tp=baseline[[time_col]][idxs],
      value=y[flag], etiv=etiv[flag], resid=res[flag], z_mad=z[flag],
      stringsAsFactors=FALSE
    )
    counts[[roi]] <- data.frame(roi=roi, n_flagged=sum(flag), stringsAsFactors=FALSE)
    # Plot
    dfp <- data.frame(etiv=etiv_ok, y=y_ok, z=z, fsid_base=baseline$fsid_base[ok])
    p <- ggplot(dfp, aes(x=etiv, y=y)) +
      geom_point(alpha=0.6, size=1.2) +
      geom_smooth(method="lm", se=FALSE, color="steelblue") +
      geom_point(data=subset(dfp, abs(z) > 3), aes(x=etiv, y=y), color="red", size=1.4)
    if (!isTRUE(opt$`no-labels`)) {
      flagged <- subset(dfp, abs(z) > 3)
      if (nrow(flagged)) {
        if (has_pkg("ggrepel")) {
          p <- p + ggrepel::geom_text_repel(data=flagged, aes(label=fsid_base), size=3, max.overlaps=50)
        } else {
          p <- p + geom_text(data=flagged, aes(label=fsid_base), size=3, nudge_y=0.02)
        }
      }
    }
    p <- p +
      labs(title=paste0(roi, " vs ", etiv_col, " (baseline)"), x=etiv_col, y=roi) +
      theme_minimal(base_size=11)
    ggsave(filename=file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", roi), "_scatter.png")), plot=p, width=5.5, height=4, dpi=120)
  }
}

outliers <- if (length(out_rows)) do.call(rbind, out_rows) else data.frame()
roi_counts <- if (length(counts)) do.call(rbind, counts) else data.frame(roi=character(0), n_flagged=integer(0))

# Subject-level counts
subj_counts <- NULL
if (nrow(outliers)) {
  subj_counts <- aggregate(roi ~ fsid_base, data=transform(outliers, roi=1L), FUN=sum)
  names(subj_counts) <- c("fsid_base", "n_flagged")
}

# Write outputs
write.csv(outliers, file=file.path(opt$outdir, "roi_etiv_outliers.csv"), row.names=FALSE)
write.csv(roi_counts, file=file.path(opt$outdir, "roi_flag_counts.csv"), row.names=FALSE)
if (!is.null(subj_counts)) write.csv(subj_counts, file=file.path(opt$outdir, "subject_flag_counts.csv"), row.names=FALSE)

# Select top ROIs to keep plots for (already saved all with flags)
if (nrow(roi_counts)) {
  ord <- roi_counts[order(-roi_counts$n_flagged), , drop=FALSE]
  top <- head(ord$roi, min(nrow(ord), opt$`plot-top`))
  kept <- file.path(plot_dir, paste0(gsub("[^A-Za-z0-9_-]", "_", top), "_scatter.png"))
  # Write a small summary file listing top plots
  con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat(sprintf("QC eTIV outliers (baseline)\n\nTotal ROIs flagged: %d\nTotal flagged points: %d\nTotal subjects flagged: %d\n\nTop ROIs by flags:\n",
              nrow(ord), sum(ord$n_flagged), if (!is.null(subj_counts)) nrow(subj_counts) else 0), file=con)
  for (r in top) {
    n <- ord$n_flagged[ord$roi==r]
    cat(sprintf("- %s: %d flags\n", r, n), file=con)
  }
}

msg("QC done. Outputs in %s\n", opt$outdir)
