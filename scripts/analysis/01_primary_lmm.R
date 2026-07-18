#!/usr/bin/env Rscript
#
# Primary confirmatory analysis: dance intervention effect on hippocampal
# subfield volumes, longitudinal mixed-effects models.
#
# For each ROI in a PRE-SPECIFIED set (whole hippocampus + DG/CA subfields,
# the a priori neurogenic-zone hypothesis), fits:
#
#   log(volume) ~ dance * time + hemisphere + age_z + sex + etiv_z
#                 + (1 + time_numeric | subject_id)
#
# with fallback to (1 | subject_id) if the random-slope model is singular.
# `dance` = pooled ballet+contemporary vs control (primary contrast, per
# study design); ballet-vs-contemporary is reported as a secondary contrast
# from the same fitted model. log(volume) is used (not raw volume) so the
# dance:time coefficient is a proportional change, comparable across
# subfields of very different absolute size -- see the "scale domination"
# lesson in the sibling 129run project's retrospective.
#
# This is the CONFIRMATORY path: the ROI set and model formula here must be
# frozen before running on real group labels. No AIC/model search here --
# see scripts/analysis/02_hierarchical_brms.R for the exploratory full-subfield
# complement with partial pooling instead of ROI pre-selection.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(broom.mixed)
})

# ------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------
option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/01_primary_lmm",
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
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))

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
tidy$is_composite <- as.logical(tidy$is_composite)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

required_tidy_cols <- c("subject_id","session","hemisphere","subfield","volume")
missing_cols <- setdiff(required_tidy_cols, names(tidy))
if (length(missing_cols)) stop(sprintf("tidy file missing required columns: %s", paste(missing_cols, collapse=", ")))

msg("Loading participant metadata from %s...\n", opt$participants)
participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

unknown_groups <- setdiff(unique(participants$group), c(dance_groups, control_group))
if (length(unknown_groups)) {
  warning(sprintf("Participants have group values not in --dance-groups/--control-group: %s (these subjects will be dropped)",
                   paste(unknown_groups, collapse=", ")))
}

# Restrict to ROI set and sum across head/body (subfield identity is already
# head/body-collapsed by the extractor's base_subfield_name(); rows differing
# only in region [head vs body] are summed here into one subfield-total per
# subject/session/hemisphere). Whole_hippocampus has a single row already.
roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=roi_dat, FUN=sum)
# eTIV is constant per subject; carry it through via a simple lookup
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

# Merge participant metadata
agg <- merge(agg, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(dance_groups, control_group), , drop=FALSE]

# Derived variables
agg$time_numeric <- as.numeric(factor(agg$session, levels=sort(unique(agg$session)))) - 1
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$dance <- factor(ifelse(agg$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
agg$style <- factor(agg$group, levels=c(control_group, dance_groups))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

n_subjects <- length(unique(agg$subject_id))
msg("Modelling data: %d rows, %d subjects, %d ROIs\n", nrow(agg), n_subjects, length(unique(agg$subfield)))
if (n_subjects < 10) warning("Fewer than 10 subjects in the merged dataset -- check subject_id matching between --tidy and --participants")

# ------------------------------------------------------------------------
# Per-ROI model fitting (pre-specified formula, no model search)
# ------------------------------------------------------------------------
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
    warning(sprintf("Could not fit any model for ROI '%s'", roi_name))
    return(NULL)
  }
  list(model=model, re_structure=re_structure)
}

roi_results <- list()
summary_rows <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting model for ROI: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]

  fit_out <- fit_one_roi(d_roi, roi_name)
  if (is.null(fit_out)) next
  model <- fit_out$model
  roi_results[[roi_name]] <- fit_out

  # Save full model summary
  sink(file.path(opt$outdir, "models", paste0(roi_name, "_model.txt")))
  cat("ROI:", roi_name, "\n")
  cat("Random-effects structure:", fit_out$re_structure, "\n\n")
  print(summary(model))
  sink()

  # Overall dance:time interaction (joint test across time levels)
  anova_tab <- tryCatch(anova(model), error=function(e) NULL)
  interaction_row <- if (!is.null(anova_tab) && "dance:time_f" %in% rownames(anova_tab)) {
    anova_tab["dance:time_f", ]
  } else NULL

  # emmeans: dance vs control contrast at each timepoint (primary readout is
  # the last timepoint, i.e. end of intervention)
  em <- tryCatch(emmeans(model, ~ dance | time_f), error=function(e) NULL)
  dance_contrasts <- if (!is.null(em)) as.data.frame(contrast(em, method="revpairwise")) else NULL
  if (!is.null(dance_contrasts)) {
    write.csv(dance_contrasts, file.path(opt$outdir, "models", paste0(roi_name, "_dance_vs_control_by_time.csv")), row.names=FALSE)
  }

  last_time <- levels(agg$time_f)[length(levels(agg$time_f))]
  primary_row <- if (!is.null(dance_contrasts)) dance_contrasts[dance_contrasts$time_f == last_time, , drop=FALSE] else NULL

  # Secondary contrast: ballet vs contemporary (style), same model refit with
  # 3-level group factor in place of the pooled dance contrast
  style_formula <- update(base_formula, paste(". ~ . -dance*time_f + style*time_f +",
                                               if (fit_out$re_structure == "random_slope") "(1 + time_numeric | subject_id)" else "(1 | subject_id)"))
  style_model <- tryCatch(lmerTest::lmer(style_formula, data=d_roi, REML=TRUE), error=function(e) NULL)
  style_contrasts <- NULL
  if (!is.null(style_model)) {
    em_style <- tryCatch(emmeans(style_model, ~ style | time_f), error=function(e) NULL)
    if (!is.null(em_style)) {
      style_contrasts <- as.data.frame(contrast(em_style, method="pairwise"))
      write.csv(style_contrasts, file.path(opt$outdir, "models", paste0(roi_name, "_style_contrasts_by_time.csv")), row.names=FALSE)
    }
  }

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    re_structure = fit_out$re_structure,
    n_obs = nrow(d_roi),
    interaction_F = if (!is.null(interaction_row)) interaction_row[["F value"]] else NA_real_,
    interaction_p = if (!is.null(interaction_row)) interaction_row[["Pr(>F)"]] else NA_real_,
    dance_effect_last_timepoint = if (!is.null(primary_row) && nrow(primary_row)) primary_row$estimate[1] else NA_real_,
    dance_effect_last_timepoint_p = if (!is.null(primary_row) && nrow(primary_row)) primary_row$p.value[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------
# Combined summary with FDR correction across the pre-specified ROI set
# ------------------------------------------------------------------------
summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$interaction_p_fdr <- p.adjust(summary_df$interaction_p, method="fdr")
  summary_df$dance_effect_last_timepoint_p_fdr <- p.adjust(summary_df$dance_effect_last_timepoint_p, method="fdr")
  summary_df$significant_interaction <- summary_df$interaction_p_fdr < opt$alpha
  summary_df$significant_dance_effect <- summary_df$dance_effect_last_timepoint_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$interaction_p_fdr), ]
  write.csv(summary_df, file.path(opt$outdir, "primary_lmm_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (FDR-corrected across %d ROIs) to primary_lmm_summary.csv\n", nrow(summary_df))
} else {
  warning("No ROI models were successfully fit; no summary written")
}

saveRDS(roi_results, file.path(opt$outdir, "all_roi_models.rds"))
msg("Done. Full model objects saved to all_roi_models.rds\n")
