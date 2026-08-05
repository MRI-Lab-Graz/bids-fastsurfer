# core_report.R -- per-ROI model dumps, FDR-corrected summary CSVs, and a
# run manifest recording what produced a given results/flex/<study>/ tree.

flex_write_model_dump <- function(outdir, roi_name, model, re_structure) {
  models_dir <- file.path(outdir, "models")
  dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
  sink(file.path(models_dir, paste0(roi_name, "_model.txt")))
  cat("ROI:", roi_name, "\n")
  cat("Random-effects structure:", re_structure, "\n\n")
  print(summary(model))
  sink()
}

# Adds an FDR-corrected column (<col>_fdr) and a significant_<col> flag for
# every column named in p_cols, sorts by the first one, and writes the CSV.
# FDR scope (within this one summary_df's rows) matches output.fdr_within =
# "measure" -- callers write one summary_df per measure so this is
# automatically per-measure correction; a caller wanting "analysis"-wide
# correction instead should combine rows across measures before calling this.
flex_write_summary <- function(outdir, summary_df, p_cols, alpha = 0.05, filename = "summary.csv") {
  if (is.null(summary_df) || !nrow(summary_df)) {
    warning(sprintf("no rows to summarize for %s; no summary written", filename))
    return(invisible(NULL))
  }
  for (col in p_cols) {
    if (!col %in% names(summary_df)) next
    fdr_col <- paste0(col, "_fdr")
    sig_col <- paste0("significant_", col)
    summary_df[[fdr_col]] <- stats::p.adjust(summary_df[[col]], method = "fdr")
    summary_df[[sig_col]] <- summary_df[[fdr_col]] < alpha
  }
  if (length(p_cols) && paste0(p_cols[1], "_fdr") %in% names(summary_df)) {
    summary_df <- summary_df[order(summary_df[[paste0(p_cols[1], "_fdr")]]), ]
  }
  path <- file.path(outdir, filename)
  utils::write.csv(summary_df, path, row.names = FALSE)
  message(sprintf("Wrote %s (%d rows)", path, nrow(summary_df)))
  invisible(summary_df)
}

flex_write_manifest <- function(cfg, outdir, test_name, extra = list()) {
  manifest <- c(
    list(
      test = test_name,
      config_path = cfg$.config_path,
      config_sha1 = tools::md5sum(cfg$.config_path)[[1]],
      run_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      r_version = R.version.string
    ),
    extra
  )
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(manifest, file.path(outdir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
}
