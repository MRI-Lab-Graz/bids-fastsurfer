#!/usr/bin/env Rscript
#
# Direct extension of the colleague's reported design: HEMISPHERE x TIME x
# GROUP x CES-D (baseline depression, continuous). Their original test was
# an rm-ANOVA restricted to the hippocampal tail; this generalizes it to
# every hippocampal subfield/composite (or amygdala nucleus/composite) via
# a mixed model instead, for two reasons: (1) LMM handles the continuous
# CES-D x categorical GROUP x categorical TIME x categorical HEMISPHERE
# four-way interaction without the sum-of-squares ambiguity classical
# rm-ANOVA has for covariate-by-within-factor terms, and (2) it uses all
# available data rather than requiring complete cases across all 3 waves.
#
#   log_volume ~ hemisphere * time_f * group3 * cesd_z + age_z + sex + etiv_z
#                + (1 | subject_id)
#
# GROUP is the 3-level solo/small-group/control factor (matching the
# colleague's design), not our usual pooled dance-vs-control binary.
# The 4-way interaction term is the primary test; FDR-corrected across
# every ROI tested, with targeted LOSO on any hit below the screen
# threshold -- this is a new, exploratory, multi-comparison sweep, so
# every one of this project's standing safeguards applies.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to a tidy volumes TSV"),
  make_option(c("-p", "--participants"), type="character", help="Path to participants TSV (subject_id, group, age, sex)"),
  make_option(c("-m", "--moderators"), type="character", help="Path to a moderators TSV with a CES-D/ADS column"),
  make_option(c("--cesd-col"), type="character", default="ads_score", help="Column name for the depression score [default %default]"),
  make_option(c("-o", "--outdir"), type="character", default="results/23_hem_time_group_cesd"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,Whole_hippocampal_head,Whole_hippocampal_body,Hippocampal_tail,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP"),
  make_option(c("--roi-column"), type="character", default="subfield"),
  make_option(c("--solo-pattern"), type="character", default="Single"),
  make_option(c("--group-pattern"), type="character", default="Group"),
  make_option(c("--control-value"), type="character", default="Control"),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (is.null(opt$moderators)) stop("--moderators is required")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$outdir, "models"), showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_col <- opt$`roi-column`

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file", roi_col))
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set")

agg <- aggregate(as.formula(paste("volume ~ subject_id + session + hemisphere +", roi_col)), data=roi_dat, FUN=sum)
names(agg)[names(agg) == roi_col] <- "subfield"
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

participants$group3 <- ifelse(grepl(opt$`solo-pattern`, participants$group), "Solo",
                        ifelse(grepl(opt$`group-pattern`, participants$group), "Group",
                        ifelse(participants$group == opt$`control-value`, "Control", NA_character_)))
if (any(is.na(participants$group3))) stop("Some subjects didn't classify into Solo/Group/Control")

moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!opt$`cesd-col` %in% names(moderators)) stop(sprintf("--cesd-col '%s' not found in --moderators", opt$`cesd-col`))
moderators <- moderators[, c("subject_id", opt$`cesd-col`)]
names(moderators)[2] <- "cesd_raw"
moderators$cesd_raw <- suppressWarnings(as.numeric(as.character(moderators$cesd_raw)))

meta <- merge(participants[, c(required_pcols, "group3")], moderators, by="subject_id", all.x=TRUE)
agg <- merge(agg, meta, by="subject_id")
agg <- agg[!is.na(agg$cesd_raw), , drop=FALSE]
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$group3 <- factor(agg$group3, levels=c("Control","Solo","Group"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$cesd_z <- as.numeric(scale(agg$cesd_raw))
agg$log_volume <- log(agg$volume)

msg("n subjects with CES-D/ADS available: %d\n", length(unique(agg$subject_id)))
msg("group3 sizes: %s\n", paste(names(table(unique(agg[,c("subject_id","group3")])$group3)),
                                  table(unique(agg[,c("subject_id","group3")])$group3), sep="=", collapse=", "))

base_formula <- log_volume ~ hemisphere * time_f * group3 * cesd_z + age_z + sex + etiv_z

fit_term_p <- function(model, term) {
  at <- tryCatch(anova(model), error=function(e) NULL)
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}
fourway_term <- "hemisphere:time_f:group3:cesd_z"

summary_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting HEM x TIME x GROUP x CES-D model for %s...\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_by_roi[[roi_name]] <- d_roi

  model <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d_roi, REML=TRUE), error=function(e) NULL)
  if (is.null(model)) { msg("  -> failed to fit\n"); next }

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_model.txt")))
  cat("ROI:", roi_name, "\n\n")
  print(summary(model))
  sink()

  p4 <- fit_term_p(model, fourway_term)
  msg("  4-way p = %s\n", if (is.na(p4)) "NA (term not estimable)" else sprintf("%.4f", p4))

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name, n_obs = nrow(d_roi), n_subjects = length(unique(d_roi$subject_id)),
    fourway_p = p4, stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (is.null(summary_df) || !nrow(summary_df)) stop("No ROI models were successfully fit")
summary_df$fourway_p_fdr <- p.adjust(summary_df$fourway_p, method="fdr")
write.csv(summary_df, file.path(opt$outdir, "summary.csv"), row.names=FALSE)
msg("\nWrote summary.csv (FDR-corrected across %d ROIs)\n", nrow(summary_df))

# ------------------------------------------------------------------------
# Targeted LOSO for hits below screen threshold
# ------------------------------------------------------------------------
screen <- summary_df[!is.na(summary_df$fourway_p) & summary_df$fourway_p < opt$`loso-screen-threshold`, ]
if (nrow(screen)) {
  msg("\nRunning targeted LOSO for %d ROI(s) with p < %.2f...\n", nrow(screen), opt$`loso-screen-threshold`)
  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat(sprintf("Targeted LOSO for HEM x TIME x GROUP x CES-D, p < %.2f\n\n", opt$`loso-screen-threshold`), file=con)
  loso_all <- list()
  for (roi_name in screen$roi) {
    d_roi <- change_by_roi[[roi_name]]
    full_p <- screen$fourway_p[screen$roi == roi_name]
    subjects <- unique(d_roi$subject_id)
    loso_p <- sapply(subjects, function(s) {
      d_sub <- d_roi[d_roi$subject_id != s, , drop=FALSE]
      if (length(unique(d_sub$group3)) < 3) return(NA_real_)
      fit <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d_sub, REML=TRUE), error=function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      fit_term_p(fit, fourway_term)
    })
    loso_p <- loso_p[!is.na(loso_p)]
    n_flips <- sum((loso_p < opt$alpha) != (full_p < opt$alpha))
    cat(sprintf("%s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance\n",
                roi_name, full_p, min(loso_p), max(loso_p), n_flips, length(loso_p)), file=con)
    loso_all[[roi_name]] <- data.frame(roi=roi_name, loso_p=loso_p)
  }
  write.csv(do.call(rbind, loso_all), file.path(opt$outdir, "loso_hits.csv"), row.names=FALSE)
  msg("Done. See %s/loso_summary.txt\n", opt$outdir)
} else {
  msg("\nNo ROI reached p < %.2f; nothing to LOSO-check.\n", opt$`loso-screen-threshold`)
}
