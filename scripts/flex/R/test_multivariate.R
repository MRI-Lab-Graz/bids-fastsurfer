# test_multivariate.R -- generic multivariate pattern analysis (MANOVA +
# PCA on covariate-residualized change scores), replacing
# 06_manova_pca_pattern.R and, via the optional region_groups looping
# below, 28_cortical_cluster_multivariate.R (run once per a priori
# functional cluster rather than once across all ROIs jointly). Does NOT
# cover 31_hippo_thickness_cluster_multivariate.R's grid-point feature
# space (individual hipsta vertices, not named ROIs) -- see
# docs/FLEX_PIPELINE.md; use scripts/analysis/31_hippo_thickness_cluster_multivariate.R
# directly for that until test_pointwise.R lands.
#
# Same age-transform note as test_change_ancova.R: 06 hardcodes raw-scale
# age_z, this engine applies design.covariates' declared transform
# uniformly instead.

flex_run_one_multivariate <- function(cfg, tidy, profile, roi_column, covariate_terms, alpha, outdir, keep_hemisphere_separate = FALSE) {
  group_cols <- if (isTRUE(keep_hemisphere_separate)) c("subject_id", "session", "hemisphere", roi_column) else c("subject_id", "session", roi_column)
  f <- stats::as.formula(paste("value ~", paste(group_cols, collapse = " + ")))
  agg <- stats::aggregate(f, data = tidy, FUN = sum)

  sessions <- sort(unique(agg$session))
  baseline_ses <- cfg$sessions$baseline %||% sessions[1]
  final_ses <- sessions[length(sessions)]

  feature_col <- roi_column
  if (isTRUE(keep_hemisphere_separate)) {
    agg$feature <- paste(agg$hemisphere, agg[[roi_column]], sep = "_")
    feature_col <- "feature"
  }
  agg$y <- flex_apply_profile_transform(agg$value, profile)

  base_vals <- agg[agg$session == baseline_ses, c("subject_id", feature_col, "y")]
  final_vals <- agg[agg$session == final_ses, c("subject_id", feature_col, "y")]
  names(base_vals)[3] <- "baseline"
  names(final_vals)[3] <- "final"
  merged_long <- merge(base_vals, final_vals, by = c("subject_id", feature_col))
  merged_long$change <- merged_long$final - merged_long$baseline

  change_wide <- stats::reshape(merged_long[, c("subject_id", feature_col, "change")],
                                 idvar = "subject_id", timevar = feature_col, direction = "wide")
  change_cols <- setdiff(names(change_wide), "subject_id")
  names(change_wide)[names(change_wide) %in% change_cols] <- sub("^change\\.", "", change_cols)
  change_cols <- setdiff(names(change_wide), "subject_id")
  if (length(change_cols) < 2) {
    warning("fewer than 2 ROI features -- MANOVA/PCA needs at least 2; skipping")
    return(NULL)
  }

  participants <- flex_load_participants(cfg)
  dat <- merge(change_wide, participants, by = "subject_id")
  dat <- flex_derive_contrast(dat, cfg)
  dat <- flex_apply_covariates(dat, cfg)

  needed <- c(change_cols, "intervention", covariate_terms)
  n_before <- nrow(dat)
  dat <- dat[stats::complete.cases(dat[, needed, drop = FALSE]), , drop = FALSE]
  message(sprintf("Subjects with complete baseline+final data across %d ROIs: %d of %d", length(change_cols), nrow(dat), n_before))
  if (nrow(dat) < length(change_cols) + 10) {
    warning(sprintf("only %d complete subjects for %d ROI features -- MANOVA/PCA at this ratio should be treated cautiously",
                     nrow(dat), length(change_cols)))
  }

  utils::write.csv(dat[, c("subject_id", "intervention", change_cols)], file.path(outdir, "change_score_matrix.csv"), row.names = FALSE)

  Y <- as.matrix(dat[, change_cols])
  cov_rhs <- paste(c("intervention", covariate_terms), collapse = " + ")
  manova_fit <- stats::manova(stats::as.formula(paste("Y ~", cov_rhs)), data = dat)
  manova_summary <- as.data.frame(summary(manova_fit, test = "Pillai")$stats)
  manova_summary$term <- rownames(manova_summary)
  utils::write.csv(manova_summary, file.path(outdir, "manova_results.csv"), row.names = FALSE)
  intervention_manova_p <- manova_summary[manova_summary$term == "intervention", "Pr(>F)"]

  covariate_rhs_only <- if (length(covariate_terms)) paste(covariate_terms, collapse = " + ") else "1"
  resid_mat <- sapply(change_cols, function(col) {
    f <- stats::as.formula(paste0("`", col, "` ~ ", covariate_rhs_only))
    stats::residuals(stats::lm(f, data = dat))
  })
  colnames(resid_mat) <- change_cols

  pca_fit <- stats::prcomp(resid_mat, center = TRUE, scale. = TRUE)
  var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  n_keep <- max(1, sum(pca_fit$sdev^2 > 1))

  loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop = FALSE])
  loadings_df$roi <- rownames(loadings_df)
  utils::write.csv(loadings_df, file.path(outdir, "pca_loadings.csv"), row.names = FALSE)
  utils::write.csv(data.frame(component = paste0("PC", seq_along(var_explained)), variance_explained = var_explained),
                    file.path(outdir, "pca_variance_explained.csv"), row.names = FALSE)

  pc_scores <- as.data.frame(pca_fit$x[, seq_len(n_keep), drop = FALSE])
  pc_scores$intervention <- dat$intervention

  pc_test_rows <- list()
  for (i in seq_len(n_keep)) {
    pc_name <- paste0("PC", i)
    fit <- stats::lm(pc_scores[[pc_name]] ~ intervention, data = pc_scores)
    tt <- summary(fit)$coefficients
    pc_test_rows[[pc_name]] <- data.frame(
      component = pc_name,
      variance_explained = var_explained[i],
      intervention_estimate = if ("interventionintervention" %in% rownames(tt)) tt["interventionintervention", "Estimate"] else NA_real_,
      intervention_p = if ("interventionintervention" %in% rownames(tt)) tt["interventionintervention", "Pr(>|t|)"] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  pc_test_df <- do.call(rbind, pc_test_rows)
  flex_write_summary(outdir, pc_test_df, p_cols = "intervention_p", alpha = alpha, filename = "pca_component_group_tests.csv")

  list(manova_p = intervention_manova_p, n_features = length(change_cols), n_subjects = nrow(dat))
}

flex_run_multivariate <- function(cfg, measure_name, spec, outdir) {
  dp <- flex_prepare_measure(cfg, measure_name)
  covariate_terms <- flex_covariate_terms(cfg)
  alpha <- flex_alpha(cfg)
  keep_hemisphere_separate <- isTRUE(spec$keep_hemisphere_separate)

  region_groups_path <- dp$measure$region_groups
  if (is.null(region_groups_path)) {
    result <- flex_run_one_multivariate(cfg, dp$tidy, dp$profile, dp$roi_column, covariate_terms, alpha, outdir, keep_hemisphere_separate)
    manova_row <- data.frame(group = measure_name, manova_p = if (!is.null(result)) result$manova_p else NA_real_,
                              n_features = if (!is.null(result)) result$n_features else NA_integer_,
                              n_subjects = if (!is.null(result)) result$n_subjects else NA_integer_)
  } else {
    clusters <- utils::read.delim(flex_resolve_path(region_groups_path), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    cluster_col <- names(clusters)[2]
    region_col <- names(clusters)[1]
    manova_rows <- list()
    for (cl in unique(clusters[[cluster_col]])) {
      rois <- clusters[clusters[[cluster_col]] == cl, region_col]
      cl_tidy <- dp$tidy[dp$tidy[[dp$roi_column]] %in% rois, , drop = FALSE]
      if (!nrow(cl_tidy)) next
      cl_outdir <- file.path(outdir, cl)
      dir.create(cl_outdir, showWarnings = FALSE, recursive = TRUE)
      message(sprintf("[multivariate/%s] cluster: %s (%d ROIs)", measure_name, cl, length(unique(rois))))
      result <- tryCatch(
        flex_run_one_multivariate(cfg, cl_tidy, dp$profile, dp$roi_column, covariate_terms, alpha, cl_outdir, keep_hemisphere_separate),
        error = function(e) { warning(sprintf("cluster '%s' failed: %s", cl, conditionMessage(e))); NULL }
      )
      manova_rows[[cl]] <- data.frame(group = cl, manova_p = if (!is.null(result)) result$manova_p else NA_real_,
                                       n_features = if (!is.null(result)) result$n_features else NA_integer_,
                                       n_subjects = if (!is.null(result)) result$n_subjects else NA_integer_)
    }
    manova_row <- do.call(rbind, manova_rows)
  }

  flex_write_summary(outdir, manova_row, p_cols = "manova_p", alpha = alpha, filename = "summary.csv")
  flex_write_manifest(cfg, outdir, test_name = "multivariate", extra = list(measure = measure_name))
  invisible(manova_row)
}
