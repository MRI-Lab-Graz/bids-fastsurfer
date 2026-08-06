# test_lmm.R -- generic confirmatory LMM battery, replacing
# 01_primary_lmm.R / 10_amygdala_lmm.R / 11_basal_ganglia_lmm.R /
# 24_thalamic_lmm.R / 25_brainstem_lmm.R / 27_cortical_thickness_lmm.R /
# 29_hippo_thickness_lmm.R -- those seven scripts differ from each other
# only in three things, which is exactly what a measure's `profile` encodes:
#   - log(volume) vs raw thickness            -> profile$transform
#   - eTIV covariate present or not           -> profile$etiv_covariate
#   - hemisphere pooled/covariate vs midline  -> profile$hemisphere
#
# For each ROI in the measure's roi_set (or every region present, if "all"):
#   y ~ intervention * time_f [+ hemisphere] + <covariates> [+ etiv_z]
#     + (1 + time_numeric | subject_id), falling back to (1 | subject_id) if
#   singular -- see core_model.R. Reports the overall intervention:time_f
#   interaction (joint test) and the intervention-vs-control contrast at the
#   final timepoint (primary readout), FDR-corrected across ROIs, plus the
#   secondary "style" contrast (individual group levels, not just pooled
#   intervention) -- exactly mirroring 01_primary_lmm.R's two-model design.

flex_run_lmm <- function(cfg, measure_name, spec, outdir) {
  dp <- flex_prepare_measure(cfg, measure_name)
  tidy <- dp$tidy
  profile <- dp$profile
  roi_column <- dp$roi_column

  agg <- flex_aggregate_roi(tidy, roi_column, profile$aggregate)

  if (isTRUE(profile$etiv_covariate)) {
    etiv_lookup <- flex_baseline_etiv(tidy, cfg$sessions$baseline)
    agg <- merge(agg, etiv_lookup, by = "subject_id", all.x = TRUE)
    agg$etiv_z <- as.numeric(scale(agg$etiv))
  }

  participants <- flex_load_participants(cfg)
  agg <- merge(agg, participants, by = "subject_id")
  agg <- flex_derive_contrast(agg, cfg)
  agg <- flex_derive_time(agg)
  agg <- flex_apply_covariates(agg, cfg)
  agg$y <- flex_apply_profile_transform(agg$value, profile)
  if (identical(profile$hemisphere, "pooled")) {
    agg$hemisphere <- factor(agg$hemisphere)
  }

  n_subjects <- length(unique(agg$subject_id))
  message(sprintf("[lmm/%s] %d rows, %d subjects, %d ROIs", measure_name, nrow(agg), n_subjects,
                   length(unique(agg[[roi_column]]))))
  if (n_subjects < 10) {
    warning("fewer than 10 subjects in the merged dataset -- check subject_id matching between the tidy file and participants")
  }

  covariate_terms <- flex_covariate_terms(cfg)
  rhs <- flex_build_rhs(profile, covariate_terms, factor_term = "intervention * time_f")
  style_rhs <- flex_build_rhs(profile, covariate_terms, factor_term = "style * time_f")

  roi_names <- unique(agg[[roi_column]])
  summary_rows <- list()
  roi_models <- list()

  for (roi_name in roi_names) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    fit_out <- flex_fit_lmm_with_fallback(d_roi, rhs)
    if (is.null(fit_out$model)) {
      warning(sprintf("could not fit any model for ROI '%s'", roi_name))
      next
    }
    model <- fit_out$model
    roi_models[[roi_name]] <- fit_out
    flex_write_model_dump(outdir, roi_name, model, fit_out$re_structure)

    anova_tab <- tryCatch(stats::anova(model), error = function(e) NULL)
    interaction_row <- if (!is.null(anova_tab) && "intervention:time_f" %in% rownames(anova_tab)) {
      anova_tab["intervention:time_f", ]
    } else NULL

    intervention_contrasts <- flex_emmeans_contrast(model, ~ intervention | time_f, method = "revpairwise")
    if (!is.null(intervention_contrasts)) {
      utils::write.csv(intervention_contrasts,
                        file.path(outdir, "models", paste0(roi_name, "_intervention_vs_control_by_time.csv")),
                        row.names = FALSE)
    }
    last_time <- levels(agg$time_f)[length(levels(agg$time_f))]
    primary_row <- if (!is.null(intervention_contrasts)) {
      intervention_contrasts[intervention_contrasts$time_f == last_time, , drop = FALSE]
    } else NULL

    style_re <- if (fit_out$re_structure == "random_slope") "(1 + time_numeric | subject_id)" else "(1 | subject_id)"
    style_model <- tryCatch(
      lmerTest::lmer(stats::as.formula(paste("y ~", style_rhs, "+", style_re)), data = d_roi, REML = TRUE),
      error = function(e) NULL
    )
    if (!is.null(style_model)) {
      style_contrasts <- flex_emmeans_contrast(style_model, ~ style | time_f, method = "pairwise")
      if (!is.null(style_contrasts)) {
        utils::write.csv(style_contrasts, file.path(outdir, "models", paste0(roi_name, "_style_contrasts_by_time.csv")),
                          row.names = FALSE)
      }
    }

    summary_rows[[roi_name]] <- data.frame(
      roi = roi_name,
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
  flex_write_summary(outdir, summary_df, p_cols = c("interaction_p", "intervention_effect_last_timepoint_p"),
                      alpha = flex_alpha(cfg), filename = "summary.csv")
  saveRDS(roi_models, file.path(outdir, "all_roi_models.rds"))
  flex_write_manifest(cfg, outdir, test_name = "lmm", extra = list(measure = measure_name, n_subjects = n_subjects, n_rois = length(roi_names)))

  invisible(summary_df)
}
