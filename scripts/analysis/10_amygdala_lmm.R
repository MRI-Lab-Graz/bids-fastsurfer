#!/usr/bin/env Rscript
#
# Amygdala nucleus analysis: same confirmatory LMM design as
# 01_primary_lmm.R, applied to amygdala nuclei (from extract_amygdala_subfields.py
# -- the "hippo-amygdala" segmentation run produces these for free alongside
# the hippocampal subfields).
#
# No pre-registered a priori nucleus subset exists for this ROI set, so all
# nuclei are tested (FDR-corrected across them), with Whole_amygdala reported
# separately as the summary measure analogous to Whole_hippocampus.
#
#   log(volume) ~ dance * time_f + hemisphere + age_z + sex + etiv_z
#                 + (1 + time_numeric | subject_id)
# with fallback to (1 | subject_id) if singular.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to amygdala_tidy.tsv (from extract_amygdala_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/10_amygdala_lmm",
              help="Output directory [default %default]"),
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
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$outdir, "models"), showWarnings=FALSE, recursive=TRUE)

dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

msg("Loading tidy amygdala data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
tidy$is_composite <- as.logical(tidy$is_composite)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(tidy, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(dance_groups, control_group), , drop=FALSE]
agg$time_numeric <- as.numeric(factor(agg$session, levels=sort(unique(agg$session)))) - 1
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$dance <- factor(ifelse(agg$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

n_subjects <- length(unique(agg$subject_id))
msg("Modelling data: %d rows, %d subjects, %d nuclei\n", nrow(agg), n_subjects, length(unique(agg$nucleus)))

base_formula <- log_volume ~ dance * time_f + hemisphere + age_z + sex + etiv_z

fit_one_roi <- function(d, roi_name) {
  fit_re <- function(re_formula) {
    full <- update(base_formula, paste(". ~ . +", re_formula))
    lmerTest::lmer(full, data=d, REML=TRUE)
  }
  model <- tryCatch(fit_re("(1 + time_numeric | subject_id)"), error=function(e) NULL)
  re_structure <- "random_slope"
  if (is.null(model) || isSingular(model, tol=1e-4)) {
    model <- tryCatch(fit_re("(1 | subject_id)"), error=function(e) NULL)
    re_structure <- "random_intercept_only"
  }
  if (is.null(model)) {
    warning(sprintf("Could not fit any model for nucleus '%s'", roi_name))
    return(NULL)
  }
  list(model=model, re_structure=re_structure)
}

summary_rows <- list()
for (roi_name in unique(agg$nucleus)) {
  msg("Fitting model for nucleus: %s\n", roi_name)
  d_roi <- agg[agg$nucleus == roi_name, , drop=FALSE]
  fit_out <- fit_one_roi(d_roi, roi_name)
  if (is.null(fit_out)) next
  model <- fit_out$model

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_model.txt")))
  cat("Nucleus:", roi_name, "-- Random-effects structure:", fit_out$re_structure, "\n\n")
  print(summary(model))
  sink()

  anova_tab <- tryCatch(anova(model), error=function(e) NULL)
  interaction_row <- if (!is.null(anova_tab) && "dance:time_f" %in% rownames(anova_tab)) anova_tab["dance:time_f", ] else NULL

  em <- tryCatch(emmeans(model, ~ dance | time_f), error=function(e) NULL)
  dance_contrasts <- if (!is.null(em)) as.data.frame(contrast(em, method="revpairwise")) else NULL
  if (!is.null(dance_contrasts)) {
    write.csv(dance_contrasts, file.path(opt$outdir, "models", paste0(roi_name, "_dance_vs_control_by_time.csv")), row.names=FALSE)
  }
  last_time <- levels(agg$time_f)[length(levels(agg$time_f))]
  primary_row <- if (!is.null(dance_contrasts)) dance_contrasts[dance_contrasts$time_f == last_time, , drop=FALSE] else NULL

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    is_composite = d_roi$is_composite[1],
    re_structure = fit_out$re_structure,
    n_obs = nrow(d_roi),
    interaction_F = if (!is.null(interaction_row)) interaction_row[["F value"]] else NA_real_,
    interaction_p = if (!is.null(interaction_row)) interaction_row[["Pr(>F)"]] else NA_real_,
    dance_effect_last_timepoint = if (!is.null(primary_row) && nrow(primary_row)) primary_row$estimate[1] else NA_real_,
    dance_effect_last_timepoint_p = if (!is.null(primary_row) && nrow(primary_row)) primary_row$p.value[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$interaction_p_fdr <- p.adjust(summary_df$interaction_p, method="fdr")
  summary_df$dance_effect_last_timepoint_p_fdr <- p.adjust(summary_df$dance_effect_last_timepoint_p, method="fdr")
  summary_df$significant_interaction <- summary_df$interaction_p_fdr < opt$alpha
  summary_df$significant_dance_effect <- summary_df$dance_effect_last_timepoint_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$interaction_p_fdr), ]
  write.csv(summary_df, file.path(opt$outdir, "amygdala_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (FDR-corrected across %d nuclei) to amygdala_summary.csv\n", nrow(summary_df))
} else {
  warning("No nucleus models were successfully fit; no summary written")
}
