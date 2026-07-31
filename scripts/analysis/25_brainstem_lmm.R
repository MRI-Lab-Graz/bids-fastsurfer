#!/usr/bin/env Rscript
#
# Brainstem substructure analysis: same confirmatory LMM design as
# 01_primary_lmm.R / 10_amygdala_lmm.R / 24_thalamic_lmm.R, applied to
# brainstem substructures (from FreeSurfer 8.2.0's
# `segment_subregions brainstem --long-base`).
#
# Unlike hippocampus/amygdala/thalamus, brainstem substructures are midline
# (single "midline" pseudo-hemisphere per FreeSurfer's own labeling, not
# lateralized lh/rh) -- so, unlike the sibling scripts, there is no
# `hemisphere` covariate here; there is exactly one row per subject/session/
# structure, not two.
#
# All 4 substructures (Medulla, Pons, Midbrain, SCP) plus Whole_brainstem
# are tested (FDR-corrected across them); no pre-registered a priori subset
# exists for this ROI set.
#
#   log(volume) ~ intervention * time_f + age_z + sex + etiv_z
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
              help="Path to brainstem_tidy.tsv (structure, hemisphere='midline', volume, etiv columns)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/25_brainstem_lmm",
              help="Output directory [default %default]"),
  make_option(c("--intervention-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'intervention' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("-m", "--moderators"), type="character", default=NULL,
              help="Optional path to a subject-level moderators TSV (one row per subject_id) -- if given, baseline CES-D is added as a standing covariate"),
  make_option(c("--cesd-col"), type="character", default="ads_score",
              help="Column name for the baseline depression score in --moderators [default %default]"),
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

intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

msg("Loading tidy brainstem data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
tidy$is_composite <- as.logical(tidy$is_composite)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(tidy, participants[, required_pcols], by="subject_id")
# eTIV varies slightly session-to-session (FreeSurfer re-estimation noise,
# not real anatomical change) -- replace the per-row (time-varying) eTIV
# from `tidy` with each subject's BASELINE (earliest session) eTIV as a
# fixed per-subject covariate.
baseline_ses_for_etiv <- sort(unique(tidy$session))[1]
etiv_lookup <- unique(tidy[tidy$session == baseline_ses_for_etiv, c("subject_id","etiv")])
agg$etiv <- NULL
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)
agg <- agg[agg$group %in% c(intervention_groups, control_group), , drop=FALSE]
agg$time_numeric <- as.numeric(factor(agg$session, levels=sort(unique(agg$session)))) - 1
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$intervention <- factor(ifelse(agg$group %in% intervention_groups, "intervention", "control"), levels=c("control","intervention"))
agg$sex <- factor(agg$sex)
# Rank-transformed (not raw) age -- see 01_primary_lmm.R for rationale: this
# cohort's age distribution has a sparse, unevenly-populated tail, and rank
# transformation caps its leverage without discarding subjects.
agg$age_z <- as.numeric(scale(rank(agg$age)))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

# Optional standing covariate: baseline depression score (e.g. CES-D/ADS),
# additive only (not crossed with intervention/time) -- moderation of the
# intervention effect by baseline depression is tested separately and more
# appropriately in 08_moderator_analysis.R, with FDR correction across ROIs.
if (!is.null(opt$moderators)) {
  if (!file.exists(opt$moderators)) stop(sprintf("moderators file not found: %s", opt$moderators))
  moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  cesd_col <- opt$`cesd-col`
  if (!cesd_col %in% names(moderators)) stop(sprintf("--cesd-col '%s' not found in --moderators", cesd_col))
  agg <- merge(agg, moderators[, c("subject_id", cesd_col)], by="subject_id", all.x=TRUE)
  names(agg)[names(agg) == cesd_col] <- "cesd_baseline"
  agg$cesd_z <- as.numeric(scale(agg$cesd_baseline))
}

n_subjects <- length(unique(agg$subject_id))
msg("Modelling data: %d rows, %d subjects, %d structures\n", nrow(agg), n_subjects, length(unique(agg$structure)))

base_formula <- log_volume ~ intervention * time_f + age_z + sex + etiv_z
if (!is.null(opt$moderators)) base_formula <- update(base_formula, . ~ . + cesd_z)

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
    warning(sprintf("Could not fit any model for structure '%s'", roi_name))
    return(NULL)
  }
  list(model=model, re_structure=re_structure)
}

summary_rows <- list()
for (roi_name in unique(agg$structure)) {
  msg("Fitting model for structure: %s\n", roi_name)
  d_roi <- agg[agg$structure == roi_name, , drop=FALSE]
  fit_out <- fit_one_roi(d_roi, roi_name)
  if (is.null(fit_out)) next
  model <- fit_out$model

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_model.txt")))
  cat("Structure:", roi_name, "-- Random-effects structure:", fit_out$re_structure, "\n\n")
  print(summary(model))
  sink()

  anova_tab <- tryCatch(anova(model), error=function(e) NULL)
  interaction_row <- if (!is.null(anova_tab) && "intervention:time_f" %in% rownames(anova_tab)) anova_tab["intervention:time_f", ] else NULL

  em <- tryCatch(emmeans(model, ~ intervention | time_f), error=function(e) NULL)
  intervention_contrasts <- if (!is.null(em)) as.data.frame(contrast(em, method="revpairwise")) else NULL
  if (!is.null(intervention_contrasts)) {
    write.csv(intervention_contrasts, file.path(opt$outdir, "models", paste0(roi_name, "_intervention_vs_control_by_time.csv")), row.names=FALSE)
  }
  last_time <- levels(agg$time_f)[length(levels(agg$time_f))]
  primary_row <- if (!is.null(intervention_contrasts)) intervention_contrasts[intervention_contrasts$time_f == last_time, , drop=FALSE] else NULL

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    is_composite = d_roi$is_composite[1],
    re_structure = fit_out$re_structure,
    n_obs = nrow(d_roi),
    interaction_F = if (!is.null(interaction_row)) interaction_row[["F value"]] else NA_real_,
    interaction_p = if (!is.null(interaction_row)) interaction_row[["Pr(>F)"]] else NA_real_,
    intervention_effect_last_timepoint = if (!is.null(primary_row) && nrow(primary_row)) primary_row$estimate[1] else NA_real_,
    intervention_effect_last_timepoint_p = if (!is.null(primary_row) && nrow(primary_row)) primary_row$p.value[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$interaction_p_fdr <- p.adjust(summary_df$interaction_p, method="fdr")
  summary_df$intervention_effect_last_timepoint_p_fdr <- p.adjust(summary_df$intervention_effect_last_timepoint_p, method="fdr")
  summary_df$significant_interaction <- summary_df$interaction_p_fdr < opt$alpha
  summary_df$significant_intervention_effect <- summary_df$intervention_effect_last_timepoint_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$interaction_p_fdr), ]
  write.csv(summary_df, file.path(opt$outdir, "brainstem_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (FDR-corrected across %d structures) to brainstem_summary.csv\n", nrow(summary_df))
} else {
  warning("No structure models were successfully fit; no summary written")
}
