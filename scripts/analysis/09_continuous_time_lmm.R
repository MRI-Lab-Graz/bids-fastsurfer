#!/usr/bin/env Rscript
#
# Continuous-time longitudinal model: replaces the categorical 3-level
# session factor (ses-1/ses-2/ses-3) with REAL elapsed days since each
# subject's own baseline scan (from scripts/extract_scan_dates.py).
#
# This matters here specifically because nominal sessions do NOT correspond
# to fixed real time: in this study, "ses-2" ranges from 14 to 32 days post
# baseline, and "ses-3" from 40 to 56 days -- treating them as one fixed
# categorical timepoint throws away real information and adds noise (two
# subjects both labelled "ses-2" may differ by over two weeks of actual
# training exposure).
#
# Model per ROI: log(volume) ~ dance * days_since_baseline + hemisphere +
#   age_z + sex + etiv_z + (1 + days_since_baseline | subject_id)
# with fallback to (1 | subject_id) if the random-slope model is singular.
#
# The dance:days_since_baseline interaction tests whether the RATE of
# volume change per day differs between dance and control -- a genuine
# growth-rate test using all available real timing information, rather
# than the coarser "does the group gap differ between nominal sessions"
# test in 01_primary_lmm.R / 05_baseline_corrected_change.R.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-d", "--scan-dates"), type="character",
              help="Path to scan_dates.tsv (from extract_scan_dates.py)"),
  make_option(c("-o", "--outdir"), type="character", default="results/09_continuous_time_lmm",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (is.null(opt$`scan-dates`)) stop("--scan-dates is required (run scripts/extract_scan_dates.py first)")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))
if (!file.exists(opt$`scan-dates`)) stop(sprintf("scan-dates file not found: %s", opt$`scan-dates`))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$outdir, "models"), showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=roi_dat, FUN=sum)
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

msg("Loading scan dates from %s...\n", opt$`scan-dates`)
scan_dates <- read.delim(opt$`scan-dates`, header=TRUE, sep="\t", stringsAsFactors=FALSE)
scan_dates$days_since_baseline <- suppressWarnings(as.numeric(scan_dates$days_since_baseline))
scan_dates <- scan_dates[!is.na(scan_dates$days_since_baseline), c("subject_id","session","days_since_baseline")]

n_before <- nrow(agg)
agg <- merge(agg, scan_dates, by=c("subject_id","session"))
msg("Rows with a valid real timestamp: %d of %d (subjects missing ses-1 are dropped, matching how baseline-referenced analyses elsewhere in this pipeline already exclude them)\n", nrow(agg), n_before)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(agg, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(dance_groups, control_group), , drop=FALSE]
agg$dance <- factor(ifelse(agg$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)
# Days as a "per month" (30-day) unit so the interaction coefficient is
# interpretable as proportional change per month rather than per day.
agg$time_months <- agg$days_since_baseline / 30

n_subjects <- length(unique(agg$subject_id))
msg("Modelling data: %d rows, %d subjects, %d ROIs\n", nrow(agg), n_subjects, length(unique(agg$subfield)))

# ------------------------------------------------------------------------
# Per-ROI continuous-time model fitting
# ------------------------------------------------------------------------
base_formula <- log_volume ~ dance * time_months + hemisphere + age_z + sex + etiv_z

fit_one_roi <- function(d, roi_name) {
  fit_re <- function(re_formula) {
    full <- update(base_formula, paste(". ~ . +", re_formula))
    lmerTest::lmer(full, data=d, REML=TRUE)
  }
  model <- tryCatch(fit_re("(1 + time_months | subject_id)"), error=function(e) NULL)
  re_structure <- "random_slope"
  if (is.null(model) || isSingular(model, tol=1e-4)) {
    model <- tryCatch(fit_re("(1 | subject_id)"), error=function(e) NULL)
    re_structure <- "random_intercept_only"
  }
  if (is.null(model)) {
    warning(sprintf("Could not fit any model for ROI '%s'", roi_name))
    return(NULL)
  }
  list(model=model, re_structure=re_structure)
}

summary_rows <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting continuous-time model for ROI: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]

  fit_out <- fit_one_roi(d_roi, roi_name)
  if (is.null(fit_out)) next
  model <- fit_out$model

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_continuous_time_model.txt")))
  cat("ROI:", roi_name, "-- Random-effects structure:", fit_out$re_structure, "\n\n")
  print(summary(model))
  sink()

  anova_tab <- tryCatch(anova(model), error=function(e) NULL)
  interaction_row <- if (!is.null(anova_tab) && "dance:time_months" %in% rownames(anova_tab)) anova_tab["dance:time_months", ] else NULL

  # emmeans slope contrast: does the per-month rate of change differ
  # between dance and control? (emtrends gives the estimated slope of
  # log_volume on time_months, per group)
  trends <- tryCatch(emtrends(model, ~ dance, var="time_months"), error=function(e) NULL)
  slope_contrast <- if (!is.null(trends)) as.data.frame(contrast(trends, method="revpairwise")) else NULL
  if (!is.null(slope_contrast)) {
    write.csv(slope_contrast, file.path(opt$outdir, "models", paste0(roi_name, "_slope_contrast.csv")), row.names=FALSE)
  }
  slope_by_group <- if (!is.null(trends)) as.data.frame(trends) else NULL
  if (!is.null(slope_by_group)) {
    write.csv(slope_by_group, file.path(opt$outdir, "models", paste0(roi_name, "_slopes_by_group.csv")), row.names=FALSE)
  }

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    re_structure = fit_out$re_structure,
    n_obs = nrow(d_roi),
    interaction_F = if (!is.null(interaction_row)) interaction_row[["F value"]] else NA_real_,
    interaction_p = if (!is.null(interaction_row)) interaction_row[["Pr(>F)"]] else NA_real_,
    slope_diff = if (!is.null(slope_contrast) && nrow(slope_contrast)) slope_contrast$estimate[1] else NA_real_,
    slope_diff_p = if (!is.null(slope_contrast) && nrow(slope_contrast)) slope_contrast$p.value[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$interaction_p_fdr <- p.adjust(summary_df$interaction_p, method="fdr")
  summary_df$significant <- summary_df$interaction_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$interaction_p_fdr), ]
  write.csv(summary_df, file.path(opt$outdir, "continuous_time_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (FDR-corrected across %d ROIs) to continuous_time_summary.csv\n", nrow(summary_df))
  msg("'slope_diff' = (dance rate of change) - (control rate of change), in log-volume units per month.\n")
} else {
  warning("No ROI models were successfully fit; no summary written")
}

saveRDS(summary_rows, file.path(opt$outdir, "all_roi_models.rds"))
msg("Done.\n")
