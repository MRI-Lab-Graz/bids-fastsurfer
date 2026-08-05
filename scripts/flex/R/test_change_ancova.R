# test_change_ancova.R -- generic baseline-corrected (ANCOVA-on-change-scores)
# battery, replacing 05_baseline_corrected_change.R (and, by the same
# change-score construction, generalizes to 15/18's pooled/proportion
# variants -- not yet ported, see docs/FLEX_PIPELINE.md).
#
# Outcome is the change score (baseline -> follow-up, on the profile's own
# transform scale -- log(volume) or raw thickness), with baseline value as
# an explicit ANCOVA covariate (Van Breukelen & van den Brand 2006): this
# makes "extra change beyond what a pre-existing baseline difference would
# predict" the direct test, rather than 01/lmm's implicit
# intervention:time_f interaction (whose emmeans readout still includes the
# baseline offset and is easy to misread).
#
# NOTE ON COVARIATE TRANSFORM: the original 05_baseline_corrected_change.R
# hardcodes raw-scale age_z (scale(age)), while 01/29/30/33 use rank_z --
# a pre-existing inconsistency across the 33-script battery, not a
# deliberate design choice. This engine applies design.covariates
# uniformly across every test (one declared transform per study), which is
# the more defensible default -- a parity check against the *unmodified*
# scripts/analysis/05_baseline_corrected_change.R will therefore only match
# to the extent its covariate transform happens to match the study's
# declared one.

flex_run_change_ancova <- function(cfg, measure_name, spec, outdir) {
  dp <- flex_prepare_measure(cfg, measure_name)
  tidy <- dp$tidy
  profile <- dp$profile
  roi_column <- dp$roi_column

  agg <- flex_aggregate_roi(tidy, roi_column)

  if (isTRUE(profile$etiv_covariate)) {
    etiv_lookup <- flex_baseline_etiv(tidy, cfg$sessions$baseline)
    agg <- merge(agg, etiv_lookup, by = "subject_id", all.x = TRUE)
    agg$etiv_z <- as.numeric(scale(agg$etiv))
  }

  participants <- flex_load_participants(cfg)
  agg <- merge(agg, participants, by = "subject_id")
  agg <- flex_derive_contrast(agg, cfg)
  agg <- flex_apply_covariates(agg, cfg)
  agg$y <- flex_apply_profile_transform(agg$value, profile)
  if (identical(profile$hemisphere, "pooled")) {
    agg$hemisphere <- factor(agg$hemisphere)
  }

  sessions <- sort(unique(agg$session))
  baseline_ses <- cfg$sessions$baseline %||% sessions[1]
  followup_sessions <- setdiff(sessions, baseline_ses)
  if (!length(followup_sessions)) stop("no follow-up sessions found beyond the baseline session")
  message(sprintf("[change_ancova/%s] baseline=%s, follow-up=%s", measure_name, baseline_ses,
                   paste(followup_sessions, collapse = ", ")))

  build_change_data <- function(d_roi) {
    base <- d_roi[d_roi$session == baseline_ses, c("subject_id", "hemisphere", "y")]
    names(base)[3] <- "y_baseline"
    fu <- d_roi[d_roi$session %in% followup_sessions, , drop = FALSE]
    merged <- merge(fu, base, by = c("subject_id", "hemisphere"))
    merged$y_change <- merged$y - merged$y_baseline
    merged$y_baseline_z <- as.numeric(scale(merged$y_baseline))
    merged$followup_f <- factor(merged$session, levels = followup_sessions)
    merged
  }

  covariate_terms <- flex_covariate_terms(cfg)
  # With only one follow-up session (the 2-timepoint design), followup_f is
  # constant and intervention:followup_f/followup_f would be aliased with
  # the intercept -- drop them rather than feeding lmer a rank-deficient
  # design (mirrors 05's own branch).
  factor_term <- if (length(followup_sessions) > 1) "intervention * followup_f" else "intervention"
  rhs <- flex_build_rhs(profile, c("y_baseline_z", covariate_terms), factor_term = factor_term)

  roi_names <- unique(agg[[roi_column]])
  summary_rows <- list()

  for (roi_name in roi_names) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    change_dat <- build_change_data(d_roi)

    model <- tryCatch(
      lmerTest::lmer(stats::as.formula(paste("y_change ~", rhs, "+ (1 | subject_id)")), data = change_dat, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(model)) {
      warning(sprintf("could not fit change-score model for ROI '%s'", roi_name))
      next
    }
    flex_write_model_dump(outdir, roi_name, model, "random_intercept_only")

    anova_tab <- tryCatch(stats::anova(model), error = function(e) NULL)
    intervention_row <- if (!is.null(anova_tab) && "intervention" %in% rownames(anova_tab)) anova_tab["intervention", ] else NULL
    interaction_row <- if (!is.null(anova_tab) && "intervention:followup_f" %in% rownames(anova_tab)) anova_tab["intervention:followup_f", ] else NULL

    intervention_contrasts <- flex_emmeans_contrast(model, ~ intervention | followup_f, method = "revpairwise")
    if (!is.null(intervention_contrasts)) {
      utils::write.csv(intervention_contrasts,
                        file.path(outdir, "models", paste0(roi_name, "_change_intervention_vs_control_by_followup.csv")),
                        row.names = FALSE)
    }

    summary_rows[[roi_name]] <- data.frame(
      roi = roi_name,
      n_obs = nrow(change_dat),
      n_subjects = length(unique(change_dat$subject_id)),
      intervention_main_F = if (!is.null(intervention_row)) intervention_row[["F value"]] else NA_real_,
      intervention_main_p = if (!is.null(intervention_row)) intervention_row[["Pr(>F)"]] else NA_real_,
      intervention_x_followup_F = if (!is.null(interaction_row)) interaction_row[["F value"]] else NA_real_,
      intervention_x_followup_p = if (!is.null(interaction_row)) interaction_row[["Pr(>F)"]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_rows)
  flex_write_summary(outdir, summary_df, p_cols = c("intervention_main_p", "intervention_x_followup_p"),
                      alpha = flex_alpha(cfg), filename = "summary.csv")
  flex_write_manifest(cfg, outdir, test_name = "change_ancova",
                       extra = list(measure = measure_name, baseline_session = baseline_ses,
                                    followup_sessions = followup_sessions))
  invisible(summary_df)
}
