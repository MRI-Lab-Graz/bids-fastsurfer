# test_loso.R -- leave-one-subject-out robustness for the two primary null-
# result checks, replacing 07_loso_robustness.R: (1) the per-ROI baseline-
# corrected change-score model's intervention main effect (same
# construction as test_change_ancova.R), and (2) the MANOVA joint pattern
# test across ROIs (same construction as test_multivariate.R). Flags
# (ROI, excluded-subject) combinations whose exclusion flips significance,
# and the single most-influential exclusion per ROI.

flex_run_loso <- function(cfg, measure_name, spec, outdir) {
  dp <- flex_prepare_measure(cfg, measure_name)
  tidy <- dp$tidy
  profile <- dp$profile
  roi_column <- dp$roi_column
  alpha <- flex_alpha(cfg)

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
  if (identical(profile$hemisphere, "pooled")) agg$hemisphere <- factor(agg$hemisphere)

  sessions <- sort(unique(agg$session))
  baseline_ses <- cfg$sessions$baseline %||% sessions[1]
  final_ses <- sessions[length(sessions)]
  followup_sessions <- setdiff(sessions, baseline_ses)
  all_subjects <- unique(agg$subject_id)
  message(sprintf("[loso/%s] %d subjects available", measure_name, length(all_subjects)))

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
  factor_term <- if (length(followup_sessions) > 1) "intervention * followup_f" else "intervention"
  rhs <- flex_build_rhs(profile, c("y_baseline_z", covariate_terms), factor_term = factor_term)
  base_formula <- stats::as.formula(paste("y_change ~", rhs, "+ (1 | subject_id)"))

  # --- Part 1: per-ROI change-score model LOSO ---------------------------
  loso_rows <- list()
  for (roi_name in unique(agg[[roi_column]])) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    change_dat_full <- build_change_data(d_roi)

    full_fit <- tryCatch(lmerTest::lmer(base_formula, data = change_dat_full, REML = TRUE), error = function(e) NULL)
    full_p <- if (!is.null(full_fit)) {
      at <- tryCatch(stats::anova(full_fit), error = function(e) NULL)
      if (!is.null(at) && "intervention" %in% rownames(at)) at["intervention", "Pr(>F)"] else NA_real_
    } else NA_real_

    for (excl_subj in all_subjects) {
      d_sub <- change_dat_full[change_dat_full$subject_id != excl_subj, , drop = FALSE]
      if (length(unique(d_sub$intervention)) != 2 || nrow(d_sub) < 10) next
      fit <- tryCatch(lmerTest::lmer(base_formula, data = d_sub, REML = TRUE), error = function(e) NULL)
      p_val <- NA_real_; est <- NA_real_
      if (!is.null(fit)) {
        at <- tryCatch(stats::anova(fit), error = function(e) NULL)
        if (!is.null(at) && "intervention" %in% rownames(at)) p_val <- at["intervention", "Pr(>F)"]
        fe <- tryCatch(lme4::fixef(fit), error = function(e) NULL)
        if (!is.null(fe) && "interventionintervention" %in% names(fe)) est <- fe[["interventionintervention"]]
      }
      loso_rows[[length(loso_rows) + 1]] <- data.frame(
        roi = roi_name, excluded_subject = excl_subj,
        intervention_main_p = p_val, intervention_main_estimate = est, full_data_p = full_p,
        stringsAsFactors = FALSE
      )
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  utils::write.csv(loso_df, file.path(outdir, "loso_change_score_by_roi.csv"), row.names = FALSE)

  flip_rows <- loso_df[!is.na(loso_df$intervention_main_p) & !is.na(loso_df$full_data_p) &
                          (loso_df$intervention_main_p < alpha) != (loso_df$full_data_p < alpha), , drop = FALSE]
  utils::write.csv(flip_rows, file.path(outdir, "loso_significance_flips.csv"), row.names = FALSE)

  influence_summary <- do.call(rbind, lapply(split(loso_df, loso_df$roi), function(g) {
    g <- g[!is.na(g$intervention_main_p), , drop = FALSE]
    if (!nrow(g)) return(NULL)
    idx <- which.max(abs(g$intervention_main_p - g$full_data_p[1]))
    data.frame(roi = g$roi[1], full_data_p = g$full_data_p[1],
               loso_p_min = min(g$intervention_main_p), loso_p_max = max(g$intervention_main_p),
               most_influential_subject = g$excluded_subject[idx],
               most_influential_p_when_excluded = g$intervention_main_p[idx], stringsAsFactors = FALSE)
  }))
  utils::write.csv(influence_summary, file.path(outdir, "loso_influence_summary.csv"), row.names = FALSE)

  # --- Part 2: MANOVA joint pattern test LOSO -----------------------------
  # Deliberately summed across hemisphere too (unlike Part 1's agg) --
  # matches 07's group_cols = c(subject_id, session, subfield), no
  # hemisphere -- one change score per subject per ROI, as MANOVA needs.
  agg_sum_formula <- stats::as.formula(paste("value ~ subject_id + session +", roi_column))
  agg_sum <- stats::aggregate(agg_sum_formula, data = tidy, FUN = sum)
  agg_sum$y <- flex_apply_profile_transform(agg_sum$value, profile)
  base_vals <- agg_sum[agg_sum$session == baseline_ses, c("subject_id", roi_column, "y")]
  final_vals <- agg_sum[agg_sum$session == final_ses, c("subject_id", roi_column, "y")]
  names(base_vals)[3] <- "baseline"; names(final_vals)[3] <- "final"
  merged_long <- merge(base_vals, final_vals, by = c("subject_id", roi_column))
  merged_long$change <- merged_long$final - merged_long$baseline
  change_wide <- stats::reshape(merged_long[, c("subject_id", roi_column, "change")],
                                 idvar = "subject_id", timevar = roi_column, direction = "wide")
  change_cols <- setdiff(names(change_wide), "subject_id")
  names(change_wide)[names(change_wide) %in% change_cols] <- sub("^change\\.", "", change_cols)
  change_cols <- setdiff(names(change_wide), "subject_id")

  manova_dat <- merge(change_wide, participants, by = "subject_id")
  manova_dat <- flex_derive_contrast(manova_dat, cfg)
  manova_dat <- flex_apply_covariates(manova_dat, cfg)
  manova_dat <- manova_dat[stats::complete.cases(manova_dat[, c(change_cols, "intervention", covariate_terms), drop = FALSE]), , drop = FALSE]

  cov_rhs <- paste(c("intervention", covariate_terms), collapse = " + ")
  run_manova_p <- function(d) {
    if (length(unique(d$intervention)) < 2 || nrow(d) < length(change_cols) + 5) return(NA_real_)
    Y <- as.matrix(d[, change_cols])
    fit <- tryCatch(stats::manova(stats::as.formula(paste("Y ~", cov_rhs)), data = d), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    s <- tryCatch(as.data.frame(summary(fit, test = "Pillai")$stats), error = function(e) NULL)
    if (is.null(s) || !"intervention" %in% rownames(s)) return(NA_real_)
    s["intervention", "Pr(>F)"]
  }

  full_manova_p <- run_manova_p(manova_dat)
  manova_loso_rows <- list()
  for (excl_subj in unique(manova_dat$subject_id)) {
    d_sub <- manova_dat[manova_dat$subject_id != excl_subj, , drop = FALSE]
    manova_loso_rows[[excl_subj]] <- data.frame(excluded_subject = excl_subj,
                                                 manova_intervention_p = run_manova_p(d_sub),
                                                 full_data_p = full_manova_p, stringsAsFactors = FALSE)
  }
  manova_loso_df <- do.call(rbind, manova_loso_rows)
  utils::write.csv(manova_loso_df, file.path(outdir, "loso_manova.csv"), row.names = FALSE)

  manova_flips <- manova_loso_df[!is.na(manova_loso_df$manova_intervention_p) &
                                    (manova_loso_df$manova_intervention_p < alpha) != (manova_loso_df$full_data_p < alpha), , drop = FALSE]
  utils::write.csv(manova_flips, file.path(outdir, "loso_manova_flips.csv"), row.names = FALSE)

  summary_df <- data.frame(
    check = c("change_score_by_roi", "manova"),
    n_flips = c(nrow(flip_rows), nrow(manova_flips)),
    n_combinations = c(nrow(loso_df), nrow(manova_loso_df)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary_df, file.path(outdir, "summary.csv"), row.names = FALSE)
  flex_write_manifest(cfg, outdir, test_name = "loso", extra = list(measure = measure_name))
  invisible(summary_df)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
