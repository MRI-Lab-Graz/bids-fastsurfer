# test_moderator.R -- generic EXPLORATORY moderator analysis, replacing
# 08_moderator_analysis.R (and, by the same interaction-term machinery,
# generalizes to 19/20's age x factorial-axis / age x dance x third-var
# three-way variants and 33's thickness retarget -- not yet ported, see
# docs/FLEX_PIPELINE.md).
#
# Tests whether the baseline-corrected change score depends on each
# declared covariate (design.covariates, e.g. age/sex) AND every
# design.moderators.baseline_cols column, via an intervention x moderator
# interaction term, adjusting for the other covariates. Hypothesis-
# generating, not confirmatory -- see the header note in 08 this ports.
# LOSO robustness runs automatically, but only for hits below
# spec$loso_screen_threshold (default 0.10), mirroring 08's triage
# rationale: full LOSO on every moderator x ROI combination would mean
# thousands of refits for combinations already clearly null.
#
# One difference from 08: this module always adjusts for the study's full
# declared covariate set as the "moderator" factor (using the same column
# 08 would derive as sex_mf, but without 08's M/F-only value filter -- for
# a cohort whose sex column only ever takes M/F values this is identical;
# for one with other values, all levels are retained instead of dropped).

flex_run_moderator <- function(cfg, measure_name, spec, outdir) {
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
  moderators <- flex_load_moderators(cfg)
  meta <- participants
  extra_moderator_cols <- character(0)
  if (!is.null(moderators)) {
    meta <- merge(meta, moderators, by = "subject_id", all.x = TRUE)
    requested <- unlist(spec$moderators)
    baseline_cols <- unlist(cfg$design$moderators$baseline_cols)
    extra_moderator_cols <- if (!is.null(requested) && !identical(requested, "all")) {
      intersect(requested, baseline_cols)
    } else {
      baseline_cols
    }
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

  base_covariate_terms <- flex_covariate_terms(cfg)
  moderator_specs <- setNames(as.list(base_covariate_terms), base_covariate_terms)
  for (col in extra_moderator_cols) moderator_specs[[col]] <- paste0(col, "_z")

  message(sprintf("[moderator/%s] testing %d candidate moderator(s) across %d ROI(s)",
                   measure_name, length(moderator_specs), length(unique(agg[[roi_column]]))))

  fit_moderator_model <- function(change_dat, mod_term) {
    adjust_terms <- setdiff(base_covariate_terms, mod_term)
    adjust_str <- if (length(adjust_terms)) paste(adjust_terms, collapse = " + ") else NULL
    followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
    hemisphere_term <- if (has_hemisphere) "hemisphere + " else ""
    re_term <- if (has_hemisphere) " + (1 | subject_id)" else ""
    etiv_term <- if (isTRUE(profile$etiv_covariate)) " + etiv_z" else ""
    rhs_tail <- paste0(hemisphere_term, followup_term, adjust_str %||% "")
    rhs_tail <- sub("\\+\\s*$", "", rhs_tail)
    f_str <- paste0("y_change ~ intervention * ", mod_term, " + y_baseline_z + ", rhs_tail, etiv_term, re_term)
    interaction_term <- paste0("intervention:", mod_term)
    if (has_hemisphere) {
      fit <- tryCatch(lmerTest::lmer(stats::as.formula(f_str), data = change_dat, REML = TRUE), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      at <- tryCatch(stats::anova(fit), error = function(e) NULL)
    } else {
      fit <- tryCatch(stats::lm(stats::as.formula(f_str), data = change_dat), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      at <- tryCatch(as.data.frame(car::Anova(fit, type = 3)), error = function(e) NULL)
    }
    if (!is.null(at) && interaction_term %in% rownames(at)) at[interaction_term, "Pr(>F)"] else NA_real_
  }

  all_results <- list()
  roi_names <- unique(agg[[roi_column]])

  for (mod_name in names(moderator_specs)) {
    mod_term <- moderator_specs[[mod_name]]
    mod_rows <- list()
    for (roi_name in roi_names) {
      d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
      change_dat <- build_change_data(d_roi)
      change_dat <- change_dat[!is.na(change_dat[[mod_term]]), , drop = FALSE]
      if (!nrow(change_dat) || length(unique(change_dat$intervention)) < 2) next
      p_val <- fit_moderator_model(change_dat, mod_term)
      mod_rows[[roi_name]] <- data.frame(moderator = mod_name, roi = roi_name, n_obs = nrow(change_dat),
                                          interaction_p = p_val, stringsAsFactors = FALSE)
    }
    mod_df <- do.call(rbind, mod_rows)
    if (!is.null(mod_df) && nrow(mod_df)) {
      mod_df$interaction_p_fdr <- stats::p.adjust(mod_df$interaction_p, method = "fdr")
      mod_df$significant <- mod_df$interaction_p_fdr < flex_alpha(cfg)
      all_results[[mod_name]] <- mod_df
    }
  }

  combined <- do.call(rbind, all_results)
  if (!is.null(combined)) {
    utils::write.csv(combined, file.path(outdir, "moderator_results.csv"), row.names = FALSE)
    utils::write.csv(combined, file.path(outdir, "summary.csv"), row.names = FALSE)
  }

  loso_threshold <- spec$loso_screen_threshold %||% 0.10
  screen_hits <- if (!is.null(combined)) combined[!is.na(combined$interaction_p) & combined$interaction_p < loso_threshold, ] else NULL
  if (!is.null(screen_hits) && nrow(screen_hits)) {
    message(sprintf("[moderator/%s] running targeted LOSO for %d hit(s) with uncorrected p < %.2f",
                     measure_name, nrow(screen_hits), loso_threshold))
    loso_rows <- list()
    for (i in seq_len(nrow(screen_hits))) {
      mod_name <- screen_hits$moderator[i]; roi_name <- screen_hits$roi[i]
      mod_term <- moderator_specs[[mod_name]]
      d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
      change_dat <- build_change_data(d_roi)
      change_dat <- change_dat[!is.na(change_dat[[mod_term]]), , drop = FALSE]
      full_p <- fit_moderator_model(change_dat, mod_term)
      for (excl_subj in unique(change_dat$subject_id)) {
        d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop = FALSE]
        if (length(unique(d_sub$intervention)) < 2) next
        p_val <- fit_moderator_model(d_sub, mod_term)
        loso_rows[[length(loso_rows) + 1]] <- data.frame(moderator = mod_name, roi = roi_name, excluded_subject = excl_subj,
                                                           loso_p = p_val, full_p = full_p, stringsAsFactors = FALSE)
      }
    }
    loso_df <- do.call(rbind, loso_rows)
    utils::write.csv(loso_df, file.path(outdir, "loso_moderator_hits.csv"), row.names = FALSE)
  }

  flex_write_manifest(cfg, outdir, test_name = "moderator", extra = list(measure = measure_name, n_moderators = length(moderator_specs)))
  invisible(combined)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
