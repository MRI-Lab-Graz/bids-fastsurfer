# test_rm_anova.R -- two engines for the "reviewer/replication-facing"
# family of tests:
#
#   engine "afex" (default) -- classical repeated-measures ANOVA (complete
#     cases only, listwise-deleted), replacing 03_rm_anova.R (between =
#     "intervention", the default) AND 22_hippotail_colleague_replication.R
#     (between = a named 3-level design.factorial axis, e.g. "group3" --
#     set spec$between to it). Raw-scale volume/thickness (afex::aov_ez's
#     own convention -- no log transform, unlike every LMM-based test in
#     this engine). Reports the highest-order interaction term (between x
#     every within factor).
#
#   engine "lmm_interaction" -- N-way factor-crossing LMM sweep with an
#     optional continuous covariate crossed in, replacing
#     23_hem_time_group_cesd.R (4-way: hemisphere x time x group3 x a
#     continuous moderator). spec$factors lists the categorical factors to
#     cross (each must be either "hemisphere"/"time_f" or a
#     design.factorial axis name); spec$continuous_moderator optionally
#     names a design.moderators column to cross in as a further continuous
#     term (raw z-scored, matching 23's cesd_z exactly -- not the study's
#     declared covariate transform). Reports the highest-order interaction
#     term, FDR-corrected across ROIs, with threshold-gated LOSO on hits
#     (spec$loso_screen_threshold, default 0.10).
#
# engine "jmv" (26_full_covariate_glm.R's jmv::anovaRM engine, chosen there
# specifically because it computes Type III SS for continuous x categorical
# interactions differently from car::Anova) is NOT implemented -- the jmv
# package isn't in this project's renv lockfile. Requesting it falls back
# to "afex" with a warning rather than failing the run.

# All orderings of a short vector (2-3 factors here) -- afex/anova term
# names aren't guaranteed to preserve input order, so the match regex needs
# every permutation, exactly as 03/22's own hardcoded alternation patterns do.
flex_permutations <- function(x) {
  if (length(x) <= 1) return(list(x))
  out <- list()
  for (i in seq_along(x)) {
    rest <- flex_permutations(x[-i])
    for (r in rest) out[[length(out) + 1]] <- c(x[i], r)
  }
  out
}

flex_run_rm_anova_afex <- function(cfg, measure_name, spec, outdir) {
  suppressPackageStartupMessages(library(afex))
  dp <- flex_prepare_measure(cfg, measure_name)
  tidy <- dp$tidy
  profile <- dp$profile
  roi_column <- dp$roi_column

  agg <- flex_aggregate_roi(tidy, roi_column)
  participants <- flex_load_participants(cfg)
  agg <- merge(agg, participants, by = "subject_id")

  between <- spec$between %||% "intervention"
  if (identical(between, "intervention")) {
    agg <- flex_derive_contrast(agg, cfg)
  } else {
    agg <- flex_derive_factorial(agg, cfg, between)
  }
  agg <- flex_derive_time(agg)
  agg$subject_id <- factor(agg$subject_id)
  has_hemisphere <- identical(profile$hemisphere, "pooled") && length(unique(agg$hemisphere)) > 1
  if (has_hemisphere) agg$hemisphere <- factor(agg$hemisphere)
  within <- if (has_hemisphere) c("time_f", "hemisphere") else "time_f"
  # Target interaction term: default is between x time_f only (03's 2-way
  # "does the group difference change over time" convention -- hemisphere
  # is fitted as a within factor but not part of the primary test unless
  # spec.interaction_factors says otherwise). 22's full 3-way replication
  # test sets spec.interaction_factors explicitly to c(between,"time_f","hemisphere").
  interaction_factors <- unlist(spec$interaction_factors) %||% c(between, "time_f")
  interaction_term <- paste(interaction_factors, collapse = ":")
  interaction_regex <- paste(vapply(flex_permutations(interaction_factors), function(x) paste(x, collapse = ".*"), character(1)), collapse = "|")

  roi_names <- unique(agg[[roi_column]])
  results_list <- list()
  for (roi_name in roi_names) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    if (!nrow(d_roi)) next
    n_before <- length(unique(d_roi$subject_id))
    fit <- tryCatch(
      afex::aov_ez(id = "subject_id", dv = "value", data = d_roi, within = within, between = between,
                    fun_aggregate = mean, anova_table = list(es = "pes")),
      error = function(e) { message(sprintf("aov_ez failed for %s: %s", roi_name, conditionMessage(e))); NULL }
    )
    if (is.null(fit)) next
    n_complete <- nrow(fit$data$long) / prod(vapply(within, function(w) nlevels(d_roi[[w]]), integer(1)))

    anova_table <- as.data.frame(fit$anova_table)
    anova_table$term <- rownames(anova_table)
    utils::write.csv(anova_table, file.path(outdir, paste0(roi_name, "_anova_table.csv")), row.names = FALSE)
    saveRDS(fit, file.path(outdir, paste0(roi_name, "_aov_fit.rds")))

    match_row <- anova_table[grepl(interaction_regex, anova_table$term), ]
    results_list[[roi_name]] <- data.frame(
      roi = roi_name, n_subjects_complete_cases = n_complete, n_subjects_available = n_before,
      interaction_p = if (nrow(match_row)) match_row[["Pr(>F)"]][1] else NA_real_,
      interaction_pes = if (nrow(match_row) && "pes" %in% names(match_row)) match_row[["pes"]][1] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, results_list)
  flex_write_summary(outdir, summary_df, p_cols = "interaction_p", alpha = flex_alpha(cfg), filename = "summary.csv")
  flex_write_manifest(cfg, outdir, test_name = "rm_anova", extra = list(measure = measure_name, engine = "afex", between = between, interaction_term = interaction_term))
  invisible(summary_df)
}

flex_run_rm_anova_lmm_interaction <- function(cfg, measure_name, spec, outdir) {
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
  agg <- flex_derive_time(agg)
  agg <- flex_apply_covariates(agg, cfg)
  agg$y <- flex_apply_profile_transform(agg$value, profile)

  factors <- unlist(spec$factors)
  if (!length(factors)) stop("spec.factors is required for engine 'lmm_interaction' (e.g. c('hemisphere','time_f','group3'))")
  for (f in factors) {
    if (f == "hemisphere") {
      agg$hemisphere <- factor(agg$hemisphere, levels = sort(unique(agg$hemisphere)))
    } else if (f == "time_f") {
      next  # already derived by flex_derive_time
    } else {
      agg <- flex_derive_factorial(agg, cfg, f)
    }
  }

  continuous_term <- NULL
  if (!is.null(spec$continuous_moderator)) {
    moderators <- flex_load_moderators(cfg)
    if (is.null(moderators)) stop("spec.continuous_moderator set but design.moderators is not configured")
    col <- spec$continuous_moderator
    if (!col %in% names(moderators)) stop(sprintf("spec.continuous_moderator '%s' not found in moderators file", col))
    agg <- merge(agg, moderators[, c("subject_id", col)], by = "subject_id", all.x = TRUE)
    continuous_term <- paste0(col, "_z")
    agg[[continuous_term]] <- as.numeric(scale(suppressWarnings(as.numeric(agg[[col]]))))
    agg <- agg[!is.na(agg[[continuous_term]]), , drop = FALSE]
  }

  all_factors <- c(factors, continuous_term)
  interaction_term <- paste(all_factors, collapse = ":")
  covariate_terms <- flex_covariate_terms(cfg)
  rhs <- flex_build_rhs(profile, covariate_terms, factor_term = paste(all_factors, collapse = " * "))
  base_formula <- stats::as.formula(paste("y ~", rhs, "+ (1 | subject_id)"))

  message(sprintf("[rm_anova/lmm_interaction/%s] testing %s across %d ROI(s), %d subjects",
                   measure_name, interaction_term, length(unique(agg[[roi_column]])), length(unique(agg$subject_id))))

  fit_term_p <- function(model, term) {
    at <- tryCatch(stats::anova(model), error = function(e) NULL)
    if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
  }

  roi_names <- unique(agg[[roi_column]])
  summary_rows <- list()
  data_by_roi <- list()

  for (roi_name in roi_names) {
    d_roi <- agg[agg[[roi_column]] == roi_name, , drop = FALSE]
    data_by_roi[[roi_name]] <- d_roi
    model <- tryCatch(lmerTest::lmer(base_formula, data = d_roi, REML = TRUE), error = function(e) NULL)
    if (is.null(model)) next
    flex_write_model_dump(outdir, roi_name, model, "random_intercept_only")
    p_val <- fit_term_p(model, interaction_term)
    summary_rows[[roi_name]] <- data.frame(roi = roi_name, n_obs = nrow(d_roi), n_subjects = length(unique(d_roi$subject_id)),
                                            interaction_p = p_val, stringsAsFactors = FALSE)
  }

  summary_df <- do.call(rbind, summary_rows)
  flex_write_summary(outdir, summary_df, p_cols = "interaction_p", alpha = flex_alpha(cfg), filename = "summary.csv")

  loso_threshold <- spec$loso_screen_threshold %||% 0.10
  screen <- if (!is.null(summary_df)) summary_df[!is.na(summary_df$interaction_p) & summary_df$interaction_p < loso_threshold, ] else NULL
  if (!is.null(screen) && nrow(screen)) {
    message(sprintf("[rm_anova/lmm_interaction/%s] targeted LOSO for %d ROI(s) with p < %.2f", measure_name, nrow(screen), loso_threshold))
    loso_all <- list()
    for (roi_name in screen$roi) {
      d_roi <- data_by_roi[[roi_name]]
      full_p <- screen$interaction_p[screen$roi == roi_name]
      loso_p <- sapply(unique(d_roi$subject_id), function(s) {
        d_sub <- d_roi[d_roi$subject_id != s, , drop = FALSE]
        if (any(vapply(factors, function(f) length(unique(d_sub[[f]])) < 2, logical(1)))) return(NA_real_)
        fit <- tryCatch(lmerTest::lmer(base_formula, data = d_sub, REML = TRUE), error = function(e) NULL)
        if (is.null(fit)) return(NA_real_)
        fit_term_p(fit, interaction_term)
      })
      loso_p <- loso_p[!is.na(loso_p)]
      loso_all[[roi_name]] <- data.frame(roi = roi_name, loso_p = loso_p, full_p = full_p)
    }
    utils::write.csv(do.call(rbind, loso_all), file.path(outdir, "loso_hits.csv"), row.names = FALSE)
  }

  flex_write_manifest(cfg, outdir, test_name = "rm_anova",
                       extra = list(measure = measure_name, engine = "lmm_interaction", factors = all_factors))
  invisible(summary_df)
}

flex_run_rm_anova <- function(cfg, measure_name, spec, outdir) {
  engine <- spec$engine %||% "afex"
  if (identical(engine, "jmv")) {
    if (!requireNamespace("jmv", quietly = TRUE)) {
      warning("engine 'jmv' requested but the jmv package is not installed -- falling back to 'afex'")
      engine <- "afex"
    }
  }
  if (identical(engine, "lmm_interaction")) {
    flex_run_rm_anova_lmm_interaction(cfg, measure_name, spec, outdir)
  } else if (identical(engine, "afex")) {
    flex_run_rm_anova_afex(cfg, measure_name, spec, outdir)
  } else {
    stop(sprintf("unknown rm_anova engine '%s' (known: afex, lmm_interaction)", engine))
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a
