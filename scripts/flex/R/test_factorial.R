# test_factorial.R -- generic within-intervention-arm factorial analysis,
# replacing 13_factorial_social_duration.R (and, by the same crossed-axis
# machinery, 30_hippo_thickness_factorial.R's thickness retarget -- not
# separately needed here, since profile already generalizes volume vs
# thickness the same way test_lmm.R does).
#
# Crosses every axis in spec$factorial_axes (default: every axis declared
# in design.factorial, e.g. "social" x "duration") on the baseline-corrected
# change score, restricted to rows classified under EVERY axis (control is
# excluded automatically, since design.factorial only ever assigns
# intervention-arm group values -- see 13's header). LOSO runs immediately
# after the main fit for every ROI, not gated behind a p-value screen (13's
# own design choice, unlike test_moderator.R's threshold-gated LOSO).

flex_run_factorial <- function(cfg, measure_name, spec, outdir) {
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

  axes <- unlist(spec$factorial_axes) %||% names(cfg$design$factorial)
  if (!length(axes)) stop("no factorial axes available: set design.factorial in the config or spec.factorial_axes")
  for (axis in axes) agg <- flex_derive_factorial(agg, cfg, axis)

  agg <- flex_apply_covariates(agg, cfg)
  agg$y <- flex_apply_profile_transform(agg$value, profile)
  has_hemisphere <- identical(profile$hemisphere, "pooled") && length(unique(agg$hemisphere)) > 1
  if (has_hemisphere) agg$hemisphere <- factor(agg$hemisphere, levels = sort(unique(agg$hemisphere)))

  message(sprintf("[factorial/%s] axes: %s -- %d subjects classified", measure_name,
                   paste(axes, collapse = " x "), length(unique(agg$subject_id))))

  sessions <- sort(unique(agg$session))
  baseline_ses <- cfg$sessions$baseline %||% sessions[1]
  followup_sessions <- setdiff(sessions, baseline_ses)

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
  factor_term <- paste(axes, collapse = " * ")
  rhs <- flex_build_rhs(profile, c("y_baseline_z", covariate_terms), factor_term = factor_term)
  if (!has_hemisphere) rhs <- sub("\\s*\\+\\s*hemisphere\\b", "", rhs)
  base_formula <- stats::as.formula(paste("y_change ~", rhs))

  fit_factorial_model <- function(d) {
    if (has_hemisphere) {
      tryCatch(lmerTest::lmer(stats::update(base_formula, . ~ . + (1 | subject_id)), data = d, REML = TRUE), error = function(e) NULL)
    } else {
      tryCatch(stats::lm(base_formula, data = d), error = function(e) NULL)
    }
  }

  fit_terms <- function(model) {
    at <- if (inherits(model, "lm")) tryCatch(as.data.frame(car::Anova(model, type = 3)), error = function(e) NULL)
          else tryCatch(stats::anova(model), error = function(e) NULL)
    if (is.null(at)) return(NULL)
    terms <- setdiff(rownames(at), c("(Intercept)", "Residuals"))
    stats::setNames(at[terms, "Pr(>F)"], terms)
  }

  roi_names <- unique(agg[[roi_column]])
  summary_rows <- list()
  change_by_roi <- list()

  for (roi_name in roi_names) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    change_dat <- build_change_data(d_roi)
    change_by_roi[[roi_name]] <- change_dat

    model <- fit_factorial_model(change_dat)
    if (is.null(model)) next
    flex_write_model_dump(outdir, roi_name, model, if (has_hemisphere) "random_intercept_only" else "lm_no_hemisphere")

    term_ps <- fit_terms(model)
    row <- data.frame(roi = roi_name, n_obs = nrow(change_dat), n_subjects = length(unique(change_dat$subject_id)),
                       stringsAsFactors = FALSE)
    for (term in names(term_ps)) row[[paste0(gsub(":", "_x_", term), "_p")]] <- term_ps[[term]]
    summary_rows[[roi_name]] <- row
  }

  summary_df <- do.call(rbind, summary_rows)
  p_cols <- if (!is.null(summary_df)) grep("_p$", names(summary_df), value = TRUE) else character(0)
  flex_write_summary(outdir, summary_df, p_cols = p_cols, alpha = flex_alpha(cfg), filename = "summary.csv")

  # LOSO robustness, run immediately for every ROI (not threshold-gated --
  # mirrors 13's design choice, unlike test_moderator.R).
  loso_rows <- list()
  for (roi_name in names(change_by_roi)) {
    change_dat <- change_by_roi[[roi_name]]
    full_model <- fit_factorial_model(change_dat)
    full_ps <- if (!is.null(full_model)) fit_terms(full_model) else NULL
    for (excl_subj in unique(change_dat$subject_id)) {
      d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop = FALSE]
      if (any(vapply(axes, function(a) length(unique(d_sub[[a]])) < 2, logical(1)))) next
      fit <- fit_factorial_model(d_sub)
      if (is.null(fit)) next
      sub_ps <- fit_terms(fit)
      row <- data.frame(roi = roi_name, excluded_subject = excl_subj, stringsAsFactors = FALSE)
      for (term in names(sub_ps)) {
        row[[paste0(gsub(":", "_x_", term), "_p")]] <- sub_ps[[term]]
        row[[paste0(gsub(":", "_x_", term), "_full_p")]] <- if (!is.null(full_ps)) full_ps[[term]] else NA_real_
      }
      loso_rows[[length(loso_rows) + 1]] <- row
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  if (!is.null(loso_df) && nrow(loso_df)) {
    utils::write.csv(loso_df, file.path(outdir, "loso_factorial.csv"), row.names = FALSE)
  }

  flex_write_manifest(cfg, outdir, test_name = "factorial", extra = list(measure = measure_name, axes = axes))
  invisible(summary_df)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
