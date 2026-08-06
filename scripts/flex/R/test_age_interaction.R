# test_age_interaction.R -- two modes for "does age's role change depending
# on something else" tests, replacing 19_age_x_social_duration.R /
# 34_hippo_thickness_age_x_social_duration.R (mode "pairwise") and
# 20_age_x_moderator_x_dance.R / 35_hippo_thickness_age_x_moderator_x_intervention.R
# (mode "threeway_intervention"). Both are profile-driven exactly like every
# other module here, so the same config works for a volume measure or a
# thickness measure without change -- 34/35 exist as separate scripts from
# 19/20 only because the original per-domain scripts don't share a formula-
# assembly layer; this module is that layer.
#
#   mode "pairwise" -- does age's relationship to change interact with each
#     of spec$factorial_axes (e.g. "social", "duration")? Restricted to
#     PAIRWISE interactions only (age x social, age x duration, social x
#     duration), no three-way age x social x duration term -- 19/34's own
#     power rationale, implemented via R's (a+b+c)^2 formula operator.
#     Reports each age x axis term separately, FDR-corrected within a
#     (axis, ROI) grid.
#
#   mode "threeway_intervention" -- does intervention x age's joint effect
#     depend on sex, or on a candidate moderator (design.moderators.
#     baseline_cols)? Tests the FULL three-way intervention x age x
#     third_var interaction, one model per third variable (not jointly),
#     exactly as 20/35 do. FDR-corrected across (third_var, ROI), with
#     targeted LOSO on hits below spec$loso_screen_threshold (default 0.10).
#
# Both modes use the study's declared age covariate (design.roles.age's
# entry in design.covariates, e.g. rank_z-transformed age_z) automatically
# -- there is no separate "age transform" setting here.

flex_age_term <- function(cfg) {
  age_col <- cfg$design$roles$age
  if (is.null(age_col)) stop("design.roles.age not set -- required for test 'age_interaction'")
  for (spec in cfg$design$covariates %||% list()) {
    if (identical(spec$col, age_col)) return(spec$as %||% spec$col)
  }
  age_col
}

flex_run_age_interaction_pairwise <- function(cfg, measure_name, spec, outdir) {
  dp <- flex_prepare_measure(cfg, measure_name)
  tidy <- dp$tidy
  profile <- dp$profile
  roi_column <- dp$roi_column
  age_term <- flex_age_term(cfg)

  agg <- flex_aggregate_roi(tidy, roi_column, profile$aggregate)
  if (isTRUE(profile$etiv_covariate)) {
    etiv_lookup <- flex_baseline_etiv(tidy, cfg$sessions$baseline)
    agg <- merge(agg, etiv_lookup, by = "subject_id", all.x = TRUE)
    agg$etiv_z <- as.numeric(scale(agg$etiv))
  }
  participants <- flex_load_participants(cfg)
  agg <- merge(agg, participants, by = "subject_id")

  axes <- unlist(spec$factorial_axes) %||% names(cfg$design$factorial)
  if (length(axes) < 2) stop("spec.factorial_axes needs at least 2 axes for mode 'pairwise' (e.g. c('social','duration'))")
  for (axis in axes) agg <- flex_derive_factorial(agg, cfg, axis)

  agg <- flex_apply_covariates(agg, cfg)
  agg$y <- flex_apply_profile_transform(agg$value, profile)
  has_hemisphere <- identical(profile$hemisphere, "pooled") && length(unique(agg$hemisphere)) > 1
  if (has_hemisphere) agg$hemisphere <- factor(agg$hemisphere, levels = sort(unique(agg$hemisphere)))

  message(sprintf("[age_interaction/pairwise/%s] axes: %s -- %d subjects classified",
                   measure_name, paste(axes, collapse = " x "), length(unique(agg$subject_id))))

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

  other_covariate_terms <- setdiff(flex_covariate_terms(cfg), age_term)
  build_formula <- function(change_dat) {
    followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
    hemisphere_term <- if (has_hemisphere) "hemisphere + " else ""
    etiv_term <- if (isTRUE(profile$etiv_covariate)) "etiv_z + " else ""
    other_term <- if (length(other_covariate_terms)) paste0(paste(other_covariate_terms, collapse = " + "), " + ") else ""
    pairwise_term <- paste0("(", paste(c(age_term, axes), collapse = " + "), ")^2")
    stats::as.formula(paste0("y_change ~ ", pairwise_term, " + y_baseline_z + ", hemisphere_term, followup_term,
                              etiv_term, other_term, "(1 | subject_id)"))
  }

  fit_term_p <- function(model, term) {
    at <- tryCatch(stats::anova(model), error = function(e) NULL)
    if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
  }
  age_axis_terms <- vapply(axes, function(a) paste(c(age_term, a), collapse = ":"), character(1))
  age_axis_regex <- setNames(vapply(age_axis_terms, function(t) {
    parts <- strsplit(t, ":")[[1]]
    paste(vapply(flex_permutations(parts), function(x) paste(x, collapse = ".*"), character(1)), collapse = "|")
  }, character(1)), axes)

  roi_names <- unique(agg[[roi_column]])
  summary_rows <- list()
  change_by_roi <- list()

  for (roi_name in roi_names) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    change_dat <- build_change_data(d_roi)
    change_by_roi[[roi_name]] <- change_dat

    model <- tryCatch(lmerTest::lmer(build_formula(change_dat), data = change_dat, REML = TRUE), error = function(e) NULL)
    if (is.null(model)) next
    flex_write_model_dump(outdir, roi_name, model, "random_intercept_only")

    at <- tryCatch(stats::anova(model), error = function(e) NULL)
    row <- data.frame(roi = roi_name, n_obs = nrow(change_dat), n_subjects = length(unique(change_dat$subject_id)),
                       stringsAsFactors = FALSE)
    for (axis in axes) {
      term_match <- if (!is.null(at)) rownames(at)[grepl(age_axis_regex[[axis]], rownames(at))] else character(0)
      row[[paste0("age_x_", axis, "_p")]] <- if (length(term_match)) at[term_match[1], "Pr(>F)"] else NA_real_
    }
    summary_rows[[roi_name]] <- row
  }

  summary_df <- do.call(rbind, summary_rows)
  p_cols <- if (!is.null(summary_df)) grep("^age_x_.*_p$", names(summary_df), value = TRUE) else character(0)
  flex_write_summary(outdir, summary_df, p_cols = p_cols, alpha = flex_alpha(cfg), filename = "summary.csv")

  loso_threshold <- spec$loso_screen_threshold %||% 0.10
  loso_rows <- list()
  if (!is.null(summary_df)) {
    for (axis in axes) {
      pcol <- paste0("age_x_", axis, "_p")
      hits <- summary_df[!is.na(summary_df[[pcol]]) & summary_df[[pcol]] < loso_threshold, "roi"]
      for (roi_name in hits) {
        change_dat <- change_by_roi[[roi_name]]
        full_model <- tryCatch(lmerTest::lmer(build_formula(change_dat), data = change_dat, REML = TRUE), error = function(e) NULL)
        full_p <- if (!is.null(full_model)) fit_term_p(full_model, grep(age_axis_regex[[axis]], rownames(stats::anova(full_model)), value = TRUE)[1]) else NA_real_
        for (excl_subj in unique(change_dat$subject_id)) {
          d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop = FALSE]
          if (any(vapply(axes, function(a) length(unique(d_sub[[a]])) < 2, logical(1)))) next
          fit <- tryCatch(lmerTest::lmer(build_formula(d_sub), data = d_sub, REML = TRUE), error = function(e) NULL)
          if (is.null(fit)) next
          at_sub <- tryCatch(stats::anova(fit), error = function(e) NULL)
          term_match <- if (!is.null(at_sub)) rownames(at_sub)[grepl(age_axis_regex[[axis]], rownames(at_sub))] else character(0)
          loso_p <- if (length(term_match)) at_sub[term_match[1], "Pr(>F)"] else NA_real_
          loso_rows[[length(loso_rows) + 1]] <- data.frame(roi = roi_name, axis = axis, excluded_subject = excl_subj,
                                                             loso_p = loso_p, full_p = full_p, stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(loso_rows)) utils::write.csv(do.call(rbind, loso_rows), file.path(outdir, "loso_hits.csv"), row.names = FALSE)

  flex_write_manifest(cfg, outdir, test_name = "age_interaction", extra = list(measure = measure_name, mode = "pairwise", axes = axes))
  invisible(summary_df)
}

flex_run_age_interaction_threeway <- function(cfg, measure_name, spec, outdir) {
  dp <- flex_prepare_measure(cfg, measure_name)
  tidy <- dp$tidy
  profile <- dp$profile
  roi_column <- dp$roi_column
  age_term <- flex_age_term(cfg)

  agg <- flex_aggregate_roi(tidy, roi_column, profile$aggregate)
  if (isTRUE(profile$etiv_covariate)) {
    etiv_lookup <- flex_baseline_etiv(tidy, cfg$sessions$baseline)
    agg <- merge(agg, etiv_lookup, by = "subject_id", all.x = TRUE)
    agg$etiv_z <- as.numeric(scale(agg$etiv))
  }

  participants <- flex_load_participants(cfg)
  moderators <- flex_load_moderators(cfg)
  meta <- participants
  extra_moderator_cols <- character(0)
  if (!is.null(moderators)) {
    meta <- merge(meta, moderators, by = "subject_id", all.x = TRUE)
    requested <- unlist(spec$moderators)
    baseline_cols <- unlist(cfg$design$moderators$baseline_cols)
    extra_moderator_cols <- if (!is.null(requested) && !identical(requested, "all")) intersect(requested, baseline_cols) else baseline_cols
  }

  agg <- merge(agg, meta, by = "subject_id")
  agg <- flex_derive_contrast(agg, cfg)
  agg <- flex_apply_covariates(agg, cfg)
  agg$y <- flex_apply_profile_transform(agg$value, profile)
  has_hemisphere <- identical(profile$hemisphere, "pooled") && length(unique(agg$hemisphere)) > 1
  if (has_hemisphere) agg$hemisphere <- factor(agg$hemisphere, levels = sort(unique(agg$hemisphere)))

  for (col in extra_moderator_cols) {
    agg[[paste0(col, "_z")]] <- as.numeric(scale(suppressWarnings(as.numeric(agg[[col]]))))
  }

  # sex (or whichever covariate is declared as a factor) is always a
  # candidate third variable, matching 20/35's "sex always tested" -- plus
  # every requested moderator column.
  third_var_specs <- list()
  sex_col <- cfg$design$roles$sex
  if (!is.null(sex_col)) {
    for (spec_c in cfg$design$covariates %||% list()) {
      if (identical(spec_c$col, sex_col) && identical(spec_c$type, "factor")) {
        third_var_specs[["sex"]] <- spec_c$as %||% spec_c$col
      }
    }
  }
  for (col in extra_moderator_cols) third_var_specs[[col]] <- paste0(col, "_z")
  if (!length(third_var_specs)) stop("no third variables available: set design.roles.sex as a factor covariate and/or design.moderators for mode 'threeway_intervention'")

  message(sprintf("[age_interaction/threeway/%s] testing %d third variable(s) across %d ROI(s)",
                   measure_name, length(third_var_specs), length(unique(agg[[roi_column]]))))

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

  fit_three_way <- function(change_dat, third_term) {
    followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
    hemisphere_term <- if (has_hemisphere) "hemisphere + " else ""
    etiv_term <- if (isTRUE(profile$etiv_covariate)) "etiv_z + " else ""
    f <- stats::as.formula(paste0("y_change ~ intervention * ", age_term, " * ", third_term,
                                   " + y_baseline_z + ", hemisphere_term, followup_term, etiv_term, "(1 | subject_id)"))
    fit <- tryCatch(lmerTest::lmer(f, data = change_dat, REML = TRUE), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    at <- tryCatch(stats::anova(fit), error = function(e) NULL)
    if (is.null(at)) return(NA_real_)
    term_regex <- paste(vapply(flex_permutations(c("intervention", age_term, third_term)),
                                function(x) paste(x, collapse = ".*"), character(1)), collapse = "|")
    match <- rownames(at)[grepl(term_regex, rownames(at))]
    if (length(match)) at[match[1], "Pr(>F)"] else NA_real_
  }

  all_results <- list()
  change_by_roi <- list()

  for (tv_name in names(third_var_specs)) {
    tv_term <- third_var_specs[[tv_name]]
    rows <- list()
    for (roi_name in unique(agg[[roi_column]])) {
      d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
      change_dat <- build_change_data(d_roi)
      change_dat <- change_dat[!is.na(change_dat[[tv_term]]), , drop = FALSE]
      change_by_roi[[paste(tv_name, roi_name, sep = "__")]] <- change_dat
      if (!nrow(change_dat) || length(unique(change_dat$intervention)) < 2) next
      p_val <- fit_three_way(change_dat, tv_term)
      rows[[roi_name]] <- data.frame(third_var = tv_name, roi = roi_name, n_obs = nrow(change_dat),
                                      interaction_p = p_val, stringsAsFactors = FALSE)
    }
    df <- do.call(rbind, rows)
    if (!is.null(df) && nrow(df)) {
      df$interaction_p_fdr <- stats::p.adjust(df$interaction_p, method = "fdr")
      df$significant <- df$interaction_p_fdr < flex_alpha(cfg)
      all_results[[tv_name]] <- df
    }
  }

  combined <- do.call(rbind, all_results)
  if (!is.null(combined)) {
    utils::write.csv(combined, file.path(outdir, "three_way_results.csv"), row.names = FALSE)
    utils::write.csv(combined, file.path(outdir, "summary.csv"), row.names = FALSE)
  }

  loso_threshold <- spec$loso_screen_threshold %||% 0.10
  screen_hits <- if (!is.null(combined)) combined[!is.na(combined$interaction_p) & combined$interaction_p < loso_threshold, ] else NULL
  if (!is.null(screen_hits) && nrow(screen_hits)) {
    message(sprintf("[age_interaction/threeway/%s] targeted LOSO for %d hit(s) with p < %.2f", measure_name, nrow(screen_hits), loso_threshold))
    loso_rows <- list()
    for (i in seq_len(nrow(screen_hits))) {
      tv_name <- screen_hits$third_var[i]; roi_name <- screen_hits$roi[i]
      tv_term <- third_var_specs[[tv_name]]
      change_dat <- change_by_roi[[paste(tv_name, roi_name, sep = "__")]]
      full_p <- fit_three_way(change_dat, tv_term)
      for (excl_subj in unique(change_dat$subject_id)) {
        d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop = FALSE]
        if (length(unique(d_sub$intervention)) < 2) next
        p_val <- fit_three_way(d_sub, tv_term)
        loso_rows[[length(loso_rows) + 1]] <- data.frame(third_var = tv_name, roi = roi_name, excluded_subject = excl_subj,
                                                           loso_p = p_val, full_p = full_p, stringsAsFactors = FALSE)
      }
    }
    utils::write.csv(do.call(rbind, loso_rows), file.path(outdir, "loso_hits.csv"), row.names = FALSE)
  }

  flex_write_manifest(cfg, outdir, test_name = "age_interaction",
                       extra = list(measure = measure_name, mode = "threeway_intervention", third_vars = names(third_var_specs)))
  invisible(combined)
}

flex_run_age_interaction <- function(cfg, measure_name, spec, outdir) {
  mode <- spec$mode %||% stop("spec.mode is required for test 'age_interaction' ('pairwise' or 'threeway_intervention')")
  if (identical(mode, "pairwise")) {
    flex_run_age_interaction_pairwise(cfg, measure_name, spec, outdir)
  } else if (identical(mode, "threeway_intervention")) {
    flex_run_age_interaction_threeway(cfg, measure_name, spec, outdir)
  } else {
    stop(sprintf("unknown age_interaction mode '%s' (known: pairwise, threeway_intervention)", mode))
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a
