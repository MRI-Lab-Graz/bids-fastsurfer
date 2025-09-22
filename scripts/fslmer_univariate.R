#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(fslmer)
  library(jsonlite)
})

# Flexible univariate LME analysis for aseg/aparc tables using fslmer
# Inputs:
#  - Qdec TSV (from scripts/generate_qdec.py)
#  - aseg/aparc table (from asegstats2table/aparcstats2table with --qdec-long)
#  - Model specification via R formula string (e.g., "~ tp*group")
#  - Random effects columns (Zcols); default: c(1,2) for intercept and first time-like term
#  - ROI column name to analyze (e.g., "Left.Hippocampus")
#  - Outputs: RDS of fit, text summary, optional CSV of merged data

option_list <- list(
  make_option(c("--config"), type="character", default=NULL, help="Path to JSON config file with arguments"),
  make_option(c("-q", "--qdec"), type="character", help="Path to qdec.table.dat (TSV)"),
  make_option(c("-a", "--aseg"), type="character", help="Path to aseg/aparc long table (TSV/whitespace)"),
  make_option(c("-r", "--roi"), type="character", help="ROI column name in aseg/aparc table (e.g., Left.Hippocampus)"),
  make_option(c("-f", "--formula"), type="character", default="~ tp", help="Model formula for fixed effects (e.g., '~ tp*group') [default %default]"),
  make_option(c("-z", "--zcols"), type="character", default="1,2", help="Comma-separated columns for random effects (e.g., '1,2' for intercept+time) [default %default]"),
  make_option(c("-c", "--contrast"), type="character", default=NULL, help="Comma-separated contrast vector (length must match ncol(X))"),
  make_option(c("-o", "--outdir"), type="character", default="fslmer_out", help="Output directory [default %default]"),
  make_option(c("--save-merged"), action="store_true", default=FALSE, help="Save merged dat as CSV"),
  make_option(c("--time-col"), type="character", default=NULL, help="Name of time variable in qdec; if NULL and 'tp' exists, uses 'tp'"),
  make_option(c("--id-col"), type="character", default=NULL, help="Override ID column in aseg (default autodetect 'Measure:volume')"),
  make_option(c("--print-cols"), action="store_true", default=FALSE, help="Print column names of qdec and aseg then exit"),
  make_option(c("--all-regions"), action="store_true", default=FALSE, help="Run analysis on all brain regions (subcortical volumes)"),
  make_option(c("--region-pattern"), type="character", default=NULL, help="Regex pattern to match ROI names (e.g., 'Hippocampus|Amygdala')"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)

opt <- parse_args(OptionParser(option_list=option_list))

msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (!is.null(opt$config)) {
  if (!file.exists(opt$config)) stop(sprintf("Config file not found: %s", opt$config))
  cfg <- jsonlite::fromJSON(opt$config, simplifyVector = TRUE)
  # Merge: config provides defaults, CLI overrides if provided
  merge_field <- function(name, default_val) {
    if (!is.null(opt[[name]])) return(opt[[name]])
    if (!is.null(cfg[[name]])) return(cfg[[name]])
    default_val
  }
  opt$qdec       <- merge_field("qdec", opt$qdec)
  opt$aseg       <- merge_field("aseg", opt$aseg)
  opt$roi        <- merge_field("roi", opt$roi)
  opt$formula    <- merge_field("formula", opt$formula)
  opt$zcols      <- merge_field("zcols", opt$zcols)
  opt$contrast   <- merge_field("contrast", opt$contrast)
  opt$outdir     <- merge_field("outdir", opt$outdir)
  opt$time_col   <- merge_field("time_col", opt$time_col)
  opt$id_col     <- merge_field("id_col", opt$id_col)
}

if (is.null(opt$qdec) || is.null(opt$aseg)) stop("qdec and aseg paths are required (via flags or config)")

# Check if we're doing multi-region analysis
multi_region <- isTRUE(opt$`all-regions`) || !is.null(opt$`region-pattern`)
if (!multi_region && is.null(opt$roi)) stop("ROI is required (e.g., Left.Hippocampus) — set via --roi, --all-regions, or --region-pattern")

if (!file.exists(opt$qdec)) stop(sprintf("qdec file not found: %s", opt$qdec))
if (!file.exists(opt$aseg)) stop(sprintf("aseg/aparc table not found: %s", opt$aseg))

# Read inputs
qdec <- tryCatch(read.delim(opt$qdec, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
if (is.null(qdec)) qdec <- read.table(opt$qdec, header=TRUE, stringsAsFactors=FALSE)
aseg <- read.table(opt$aseg, header=TRUE, stringsAsFactors=FALSE)

if (isTRUE(opt$print_cols)) {
  cat("qdec columns:\n"); print(names(qdec))
  cat("aseg/aparc columns:\n"); print(names(aseg))
  quit(status=0)
}

# Identify ID column in aseg
id_col <- opt$id_col
if (is.null(id_col)) {
  cand <- grep("^Measure[.:]volume$", names(aseg), ignore.case=TRUE, value=TRUE)
  if (length(cand) == 0) stop("Could not find 'Measure:volume' column in aseg/aparc table; use --id-col to override")
  id_col <- cand[[1]]
}

# Extract fsid and fsid.base
aseg$fsid      <- sub("\\.long\\..*", "",  aseg[[id_col]])
aseg$fsid.base <- sub(".*\\.long\\.",  "", aseg[[id_col]])

# Ensure time column
time_col <- opt$time_col
if (is.null(time_col)) {
  if ("Time.From.Baseline" %in% names(qdec)) time_col <- "Time.From.Baseline"
  else if ("tp" %in% names(qdec)) time_col <- "tp"
  else stop("Could not infer time column. Set --time-col (e.g., tp or Time.From.Baseline)")
}
if (!(time_col %in% names(qdec))) stop(sprintf("Time column '%s' not found in qdec", time_col))

# Handle different naming conventions for base ID column
base_id_col <- "fsid.base"
if ("fsid-base" %in% names(qdec) && !("fsid.base" %in% names(qdec))) {
  base_id_col <- "fsid-base"
  names(aseg)[names(aseg) == "fsid.base"] <- "fsid-base"
}

# Merge
dat <- merge(qdec, aseg, by=c(base_id_col, "fsid"))
if (nrow(dat) == 0) stop("Merged data is empty; check IDs and inputs")

# Order by subject then time (following fslmer documentation)
dat <- dat[order(dat[, "fsid-base"], dat[, time_col]), ]

# Create vector of observations per subject (following fslmer documentation) 
ni <- matrix(unname(table(dat[, "fsid-base"])), ncol=1)

# Function to analyze a single ROI
analyze_roi <- function(roi_name, dat, ni, opt, time_col) {
  msg <- function(...) if (!isTRUE(opt$quiet)) cat(sprintf(...))
  
  # Check if ROI exists in data
  if (!(roi_name %in% names(dat))) {
    # Try to adapt: replace '-' with '.' in ROI name (R converts hyphens to dots)
    roi2 <- gsub("-", ".", roi_name, fixed=TRUE)
    if (!(roi2 %in% names(dat))) {
      return(list(error = sprintf("ROI column '%s' not found", roi_name)))
    }
    roi_name <- roi2
  }
  
  msg("Analyzing ROI: %s\n", roi_name)
  
  # Extract ROI values as column vector (following fslmer documentation)
  Y <- matrix(dat[[roi_name]], ncol=1)
  
  # Add 'y' column to data frame for formula parsing
  dat$y <- dat[[roi_name]]
  
  # Design matrix from formula (fixed effects)
  form <- as.formula(opt$formula)
  X <- model.matrix(form, dat)
  
  # Random effects columns (support JSON arrays or comma-separated string)
  if (is.numeric(opt$zcols)) {
    zcols <- as.integer(opt$zcols)
  } else {
    zcols <- as.integer(strsplit(opt$zcols, ",")[[1]])
  }
  if (any(is.na(zcols))) return(list(error = "--zcols must be a comma-separated list of integers"))
  
  msg("Design matrix: %d x %d\n", nrow(X), ncol(X))
  msg("Random effects columns (Zcols): %s\n", paste(zcols, collapse=","))
  msg("ROI: %s; observations: %d; subjects: %d\n", roi_name, nrow(Y), length(ni))
  
  # Fit model
  tryCatch({
    stats <- lme_fit_FS(X, zcols, Y, ni)
    
    # Optional contrast (support JSON array or comma-separated string)
    F_C <- NULL
    if (!is.null(opt$contrast)) {
      if (is.numeric(opt$contrast)) {
        Cvec <- as.numeric(opt$contrast)
      } else {
        Cvec <- as.numeric(strsplit(opt$contrast, ",")[[1]])
      }
      if (length(Cvec) != ncol(X)) return(list(error = sprintf("Contrast length (%d) must match ncol(X) (%d)", length(Cvec), ncol(X))))
      C <- matrix(Cvec, nrow=1)
      F_C <- lme_F(stats, C)
    }
    
    return(list(
      roi = roi_name,
      stats = stats,
      F_C = F_C,
      X = X,
      Y = Y,
      zcols = zcols,
      Cvec = if(!is.null(opt$contrast)) Cvec else NULL
    ))
  }, error = function(e) {
    return(list(error = sprintf("Model fitting failed for %s: %s", roi_name, e$message)))
  })
}

# Determine which ROIs to analyze
if (multi_region) {
  # Get all potential brain region columns
  brain_cols <- names(dat)[!names(dat) %in% c("fsid-base", "fsid", time_col, names(qdec))]
  
  if (isTRUE(opt$`all-regions`)) {
    # Use all subcortical/brain regions (exclude summary measures)
    rois_to_analyze <- brain_cols[grepl("^(Left|Right)\\.", brain_cols) | 
                                  grepl("(Hippocampus|Amygdala|Thalamus|Caudate|Putamen|Pallidum|Accumbens|Ventricle|Stem)", brain_cols, ignore.case=TRUE)]
  } else {
    # Use pattern matching
    pattern <- opt$`region-pattern`
    rois_to_analyze <- brain_cols[grepl(pattern, brain_cols, ignore.case=TRUE, perl=TRUE)]
  }
  
  if (length(rois_to_analyze) == 0) {
    stop("No ROIs found matching criteria. Available brain regions: ", paste(head(brain_cols, 10), collapse=", "))
  }
  
  cat(sprintf("Found %d ROIs to analyze: %s\n", length(rois_to_analyze), paste(head(rois_to_analyze, 5), collapse=", ")))
  if (length(rois_to_analyze) > 5) cat("... and", length(rois_to_analyze) - 5, "more\n")
  
} else {
  # Single ROI analysis
  rois_to_analyze <- opt$roi
}

# Run analysis on all selected ROIs
results <- list()
failed_rois <- c()

for (roi in rois_to_analyze) {
  result <- analyze_roi(roi, dat, ni, opt, time_col)
  
  if (!is.null(result$error)) {
    cat("ERROR:", result$error, "\n")
    failed_rois <- c(failed_rois, roi)
    next
  }
  
  results[[roi]] <- result
}

if (length(results) == 0) {
  stop("No ROIs were successfully analyzed")
}

cat(sprintf("Successfully analyzed %d ROIs", length(results)))
if (length(failed_rois) > 0) {
  cat(sprintf(" (%d failed: %s)", length(failed_rois), paste(head(failed_rois, 3), collapse=", ")))
}
cat("\n")

# Outputs
outdir <- opt$outdir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

if (multi_region) {
  # Multi-region outputs: create summary tables and individual results
  
  # Create summary table of coefficients
  coef_summary <- do.call(rbind, lapply(names(results), function(roi) {
    r <- results[[roi]]
    data.frame(
      roi = roi,
      coef_name = rownames(r$stats$Bhat),
      beta = as.numeric(r$stats$Bhat),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(coef_summary, file=file.path(outdir, "lme_coefficients.csv"), row.names=FALSE)
  
  # Create summary table of F-tests (if contrasts were used)
  if (!is.null(opt$contrast)) {
    f_summary <- do.call(rbind, lapply(names(results), function(roi) {
      r <- results[[roi]]
      if (!is.null(r$F_C)) {
        data.frame(
          roi = roi,
          F = r$F_C$F,
          pval = r$F_C$pval,
          sgn = r$F_C$sgn,
          df = r$F_C$df,
          stringsAsFactors = FALSE
        )
      }
    }))
    write.csv(f_summary, file=file.path(outdir, "lme_F_test.csv"), row.names=FALSE)
  }
  
  # Save individual results
  roi_dir <- file.path(outdir, "individual_rois")
  dir.create(roi_dir, showWarnings=FALSE, recursive=TRUE)
  
  for (roi in names(results)) {
    r <- results[[roi]]
    safe_roi_name <- gsub("[^A-Za-z0-9_-]", "_", roi)
    
    saveRDS(r$stats, file=file.path(roi_dir, paste0(safe_roi_name, "_fit.rds")))
    write.table(data.frame(Bhat=I(r$stats$Bhat)), 
                file=file.path(roi_dir, paste0(safe_roi_name, "_Bhat.txt")), 
                quote=FALSE, sep="\t")
    
    if (!is.null(r$F_C)) {
      write.table(data.frame(F=r$F_C$F, pval=r$F_C$pval, sgn=r$F_C$sgn, df=r$F_C$df),
                  file=file.path(roi_dir, paste0(safe_roi_name, "_F_test.txt")), 
                  quote=FALSE, sep="\t", row.names=FALSE)
    }
  }
  
} else {
  # Single ROI outputs (backward compatibility)
  r <- results[[1]]
  saveRDS(r$stats, file=file.path(outdir, "fit.rds"))
  write.table(data.frame(Bhat=I(r$stats$Bhat)), file=file.path(outdir, "Bhat.txt"), quote=FALSE, sep="\t")
  
  if (!is.null(r$F_C)) {
    write.table(data.frame(F=r$F_C$F, pval=r$F_C$pval, sgn=r$F_C$sgn, df=r$F_C$df),
                file=file.path(outdir, "F_test.txt"), quote=FALSE, sep="\t", row.names=FALSE)
  }
}

if (isTRUE(opt$save_merged)) {
  write.csv(dat, file=file.path(outdir, "merged_data.csv"), row.names=FALSE)
}

# Save the resolved config for reproducibility
resolved <- list(
  qdec = opt$qdec,
  aseg = opt$aseg,
  roi = if(multi_region) paste(rois_to_analyze, collapse=",") else opt$roi,
  formula = opt$formula,
  zcols = if(multi_region) results[[1]]$zcols else results[[1]]$zcols,
  contrast = if (!is.null(opt$contrast)) results[[1]]$Cvec else NULL,
  outdir = outdir,
  time_col = time_col,
  id_col = opt$id_col,
  multi_region = multi_region,
  n_rois_analyzed = length(results),
  n_rois_failed = length(failed_rois)
)
jsonlite::write_json(resolved, file.path(outdir, "used_config.json"), pretty = TRUE, auto_unbox = TRUE)

msg("Done. Outputs written to %s\n", outdir)
