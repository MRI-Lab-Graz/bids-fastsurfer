#!/usr/bin/env Rscript
#
# Primary confirmatory analysis: intervention effect on hippocampal subfield
# THICKNESS (hipsta shape/thickness analysis), longitudinal mixed-effects
# models, one per subfield.
#
#   thickness ~ intervention * time + hemisphere + age_z + sex
#               + (1 + time_numeric | subject_id)
#
# with fallback to (1 | subject_id) if the random-slope model is singular.
# This mirrors 27_cortical_thickness_lmm.R exactly (not 01_primary_lmm.R):
# like cortical thickness, hippocampal thickness is not scaled by head size,
# so eTIV is intentionally NOT included as a covariate, and thickness is
# modelled on its raw mm scale (no log-transform -- the "scale domination"
# rationale for log(volume) in 01_primary_lmm.R doesn't apply across
# subfields sharing the same ~1-3mm thickness scale). See
# scripts/extract_hipsta_thickness.py's docstring for why the tidy file's
# columns ("region", "value") were deliberately made to match
# aparc_a2009s_tidy_2tp_fs82.tsv's contract rather than
# hippo_subfields_tidy.tsv's ("subfield", "volume").
#
# Unlike 27, --group-column/--intervention-groups/--control-group default to
# this repo's actual configs/participants.tsv encoding for study 129
# (group_5: control/alone_2w/alone_4w/smallgroup_2w/smallgroup_4w) rather than
# a placeholder scheme -- override for other participants files/studies.
#
# For the question of whether the effect differs by social context (alone vs.
# group) or duration (2wk vs 4wk) specifically, see
# 30_hippo_thickness_factorial.R -- this script only tests the pooled
# any-intervention-vs-control contrast.
#
# Reports BOTH uncorrected and FDR-corrected (Benjamini-Hochberg, across all
# tested subfields) p-values in the same summary table, matching
# 27_cortical_thickness_lmm.R's convention.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_thickness_tidy.tsv (from extract_hipsta_thickness.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, <group-column>, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/29_hippo_thickness_lmm",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character", default="presubiculum,subiculum,CA1,CA2_CA3",
              help="Comma-separated subfield names, or 'all' for every labelled subfield (excludes 'unlabeled') [default %default]"),
  make_option(c("--group-column"), type="character", default="group_5",
              help="Column in --participants holding the group assignment [default %default]"),
  make_option(c("--intervention-groups"), type="character",
              default="alone_2w,alone_4w,smallgroup_2w,smallgroup_4w",
              help="Comma-separated group values pooled into the 'intervention' contrast [default %default]"),
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

intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)
group_col <- opt$`group-column`

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
msg("Loading tidy hippocampal thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$value <- suppressWarnings(as.numeric(tidy$value))

required_tidy_cols <- c("subject_id","session","hemisphere","region","value")
missing_cols <- setdiff(required_tidy_cols, names(tidy))
if (length(missing_cols)) stop(sprintf("tidy file missing required columns: %s", paste(missing_cols, collapse=", ")))

roi_set <- if (identical(opt$`roi-set`, "all")) {
  setdiff(sort(unique(tidy$region)), "unlabeled")
} else {
  trimws(strsplit(opt$`roi-set`, ",")[[1]])
}

msg("Loading participant metadata from %s...\n", opt$participants)
participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"

unknown_groups <- setdiff(unique(participants$group), c(intervention_groups, control_group))
if (length(unknown_groups)) {
  msg("Note: group values not in --intervention-groups/--control-group (subjects dropped): %s\n",
      paste(unknown_groups, collapse=", "))
}

roi_dat <- tidy[tidy$region %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file (e.g. presubiculum, subiculum, CA1, CA2_CA3)")

agg <- aggregate(value ~ subject_id + session + hemisphere + region, data=roi_dat, FUN=mean)
agg <- merge(agg, participants[, c("subject_id","group","age","sex")], by="subject_id")
agg <- agg[agg$group %in% c(intervention_groups, control_group), , drop=FALSE]

agg$time_numeric <- as.numeric(factor(agg$session, levels=sort(unique(agg$session)))) - 1
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$intervention <- factor(ifelse(agg$group %in% intervention_groups, "intervention", "control"), levels=c("control","intervention"))
agg$hemisphere <- factor(agg$hemisphere, levels=sort(unique(agg$hemisphere)))
agg$sex <- factor(agg$sex)
# Rank-transformed age (project-wide standard -- see 01_primary_lmm.R for
# the rationale: this cohort's age distribution is concentrated at 18-27
# with a sparse, unevenly-populated tail out to 38).
agg$age_z <- as.numeric(scale(rank(agg$age)))

n_subjects <- length(unique(agg$subject_id))
msg("Modelling data: %d rows, %d subjects, %d subfields\n", nrow(agg), n_subjects, length(unique(agg$region)))
if (n_subjects < 10) warning("Fewer than 10 subjects in the merged dataset -- check subject_id matching between --tidy and --participants")

base_formula <- value ~ intervention * time_f + hemisphere + age_z + sex

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
    warning(sprintf("Could not fit any model for subfield '%s'", roi_name))
    return(NULL)
  }
  list(model=model, re_structure=re_structure)
}

roi_results <- list()
summary_rows <- list()

for (roi_name in roi_set) {
  msg("Fitting model for subfield: %s\n", roi_name)
  d_roi <- agg[agg$region == roi_name, , drop=FALSE]
  if (!nrow(d_roi)) { msg("  (no data, skipping)\n"); next }

  fit_out <- fit_one_roi(d_roi, roi_name)
  if (is.null(fit_out)) next
  model <- fit_out$model
  roi_results[[roi_name]] <- fit_out

  sink(file.path(opt$outdir, "models", paste0(gsub("[^0-9A-Za-z_]", "_", roi_name), "_model.txt")))
  cat("Subfield:", roi_name, "\n")
  cat("Random-effects structure:", fit_out$re_structure, "\n\n")
  print(summary(model))
  sink()

  anova_tab <- tryCatch(anova(model), error=function(e) NULL)
  interaction_row <- if (!is.null(anova_tab) && "intervention:time_f" %in% rownames(anova_tab)) {
    anova_tab["intervention:time_f", ]
  } else NULL

  em <- tryCatch(emmeans(model, ~ intervention | time_f), error=function(e) NULL)
  intervention_contrasts <- if (!is.null(em)) as.data.frame(contrast(em, method="revpairwise")) else NULL
  if (!is.null(intervention_contrasts)) {
    write.csv(intervention_contrasts, file.path(opt$outdir, "models", paste0(gsub("[^0-9A-Za-z_]", "_", roi_name), "_intervention_vs_control_by_time.csv")), row.names=FALSE)
  }

  last_time <- levels(agg$time_f)[length(levels(agg$time_f))]
  primary_row <- if (!is.null(intervention_contrasts)) intervention_contrasts[intervention_contrasts$time_f == last_time, , drop=FALSE] else NULL

  summary_rows[[roi_name]] <- data.frame(
    subfield = roi_name,
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
  summary_df$significant_interaction_uncorrected <- summary_df$interaction_p < opt$alpha
  summary_df$significant_interaction_fdr <- summary_df$interaction_p_fdr < opt$alpha
  summary_df$significant_intervention_effect_uncorrected <- summary_df$intervention_effect_last_timepoint_p < opt$alpha
  summary_df$significant_intervention_effect_fdr <- summary_df$intervention_effect_last_timepoint_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$interaction_p), ]
  write.csv(summary_df, file.path(opt$outdir, "hippo_thickness_lmm_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (both uncorrected and FDR-corrected across %d subfields) to hippo_thickness_lmm_summary.csv\n", nrow(summary_df))
} else {
  warning("No subfield models were successfully fit; no summary written")
}

saveRDS(roi_results, file.path(opt$outdir, "all_roi_models.rds"))
msg("Done. Full model objects saved to all_roi_models.rds\n")
