#!/usr/bin/env Rscript
#
# Primary confirmatory analysis: intervention effect on cortical thickness,
# Destrieux atlas (aparc.a2009s), ALL 74 regions (per-region model pools both
# hemispheres via a `hemisphere` covariate, exactly as in 01_primary_lmm.R --
# one test per region, not per region-per-hemisphere).
#
#   thickness ~ intervention * time + hemisphere + age_z + sex
#               + (1 + time_numeric | subject_id)
#
# with fallback to (1 | subject_id) if the random-slope model is singular.
# Unlike subcortical volume, cortical thickness is NOT scaled by head size,
# so eTIV is intentionally NOT included as a covariate here (standard
# practice; see extract_freesurfer.py's docstring for the same point) and
# thickness is modelled on its raw mm scale (no log-transform -- the
# "scale domination" rationale for log(volume) in 01_primary_lmm.R doesn't
# apply here since all 74 regions share the same ~1-4mm thickness scale).
#
# Reports BOTH uncorrected and FDR-corrected (Benjamini-Hochberg, across all
# 74 regions) p-values in the same summary table -- no --fdr flag needed,
# matching 01_primary_lmm.R's convention.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to aparc_a2009s_tidy_2tp_fs82.tsv (from extract_freesurfer.py --source aparc.a2009s)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/27_cortical_thickness_lmm",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character", default="all",
              help="Comma-separated Destrieux region names, or 'all' for every region [default %default]"),
  make_option(c("--intervention-groups"), type="character", default="Single_2wk,Single_4wk,Group_2wk,Group_4wk",
              help="Comma-separated group values pooled into the 'intervention' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="Control",
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

msg("Loading tidy cortical thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$is_composite <- as.logical(tidy$is_composite)
tidy$value <- suppressWarnings(as.numeric(tidy$value))

required_tidy_cols <- c("subject_id","session","hemisphere","region","value")
missing_cols <- setdiff(required_tidy_cols, names(tidy))
if (length(missing_cols)) stop(sprintf("tidy file missing required columns: %s", paste(missing_cols, collapse=", ")))

roi_set <- if (identical(opt$`roi-set`, "all")) sort(unique(tidy$region)) else trimws(strsplit(opt$`roi-set`, ",")[[1]])

msg("Loading participant metadata from %s...\n", opt$participants)
participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

unknown_groups <- setdiff(unique(participants$group), c(intervention_groups, control_group))
if (length(unknown_groups)) {
  warning(sprintf("Participants have group values not in --intervention-groups/--control-group: %s (these subjects will be dropped)",
                   paste(unknown_groups, collapse=", ")))
}

roi_dat <- tidy[tidy$region %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check region names against the tidy file")

agg <- aggregate(value ~ subject_id + session + hemisphere + region, data=roi_dat, FUN=mean)
agg <- merge(agg, participants[, required_pcols], by="subject_id")
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
msg("Modelling data: %d rows, %d subjects, %d regions\n", nrow(agg), n_subjects, length(unique(agg$region)))
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
    warning(sprintf("Could not fit any model for region '%s'", roi_name))
    return(NULL)
  }
  list(model=model, re_structure=re_structure)
}

roi_results <- list()
summary_rows <- list()

for (roi_name in roi_set) {
  msg("Fitting model for region: %s\n", roi_name)
  d_roi <- agg[agg$region == roi_name, , drop=FALSE]
  if (!nrow(d_roi)) { msg("  (no data, skipping)\n"); next }

  fit_out <- fit_one_roi(d_roi, roi_name)
  if (is.null(fit_out)) next
  model <- fit_out$model
  roi_results[[roi_name]] <- fit_out

  sink(file.path(opt$outdir, "models", paste0(gsub("[^0-9A-Za-z_]", "_", roi_name), "_model.txt")))
  cat("Region:", roi_name, "\n")
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
    region = roi_name,
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
  write.csv(summary_df, file.path(opt$outdir, "cortical_thickness_lmm_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (both uncorrected and FDR-corrected across %d regions) to cortical_thickness_lmm_summary.csv\n", nrow(summary_df))
} else {
  warning("No region models were successfully fit; no summary written")
}

saveRDS(roi_results, file.path(opt$outdir, "all_roi_models.rds"))
msg("Done. Full model objects saved to all_roi_models.rds\n")
