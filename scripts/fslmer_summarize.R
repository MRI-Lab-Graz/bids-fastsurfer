#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(ggplot2)
})

option_list <- list(
  make_option("--results-dir", type = "character", default = NULL,
              help = "Directory containing fslmer outputs (used_config.json, individual_rois/, etc.)"),
  make_option("--effect", type = "character", default = "tp:sexM",
              help = "Coefficient name to summarise (default: 'tp:sexM')"),
  make_option("--alpha", type = "double", default = 0.05,
              help = "Significance threshold (default: 0.05)"),
  make_option("--trend-alpha", type = "double", default = 0.1,
              help = "Trend-level threshold (default: 0.10)"),
  make_option("--html", type = "character", default = NULL,
              help = "Optional HTML report path (default: <results-dir>/summary.html)"),
  make_option("--plot-dir", type = "character", default = NULL,
              help = "Directory to store plots (default: <results-dir>/plots)"),
  make_option("--verbose", action = "store_true", default = FALSE,
              help = "Print progress messages")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$`results-dir`)) stop("--results-dir is required")
results_dir <- normalizePath(opt$`results-dir`, mustWork = TRUE)
html_path <- if (is.null(opt$html)) file.path(results_dir, "summary.html") else opt$html
plot_dir <- if (is.null(opt$`plot-dir`)) file.path(results_dir, "plots") else opt$`plot-dir`
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

used_cfg_path <- file.path(results_dir, "used_config.json")
if (!file.exists(used_cfg_path)) stop(sprintf("used_config.json not found in %s", results_dir))
cfg <- jsonlite::fromJSON(used_cfg_path, simplifyVector = TRUE)

coerce_numeric_like <- function(x) {
  if (!is.character(x)) return(x)
  stripped <- trimws(tolower(x))
  stripped[stripped %in% c("", "na", "n/a", "nan", "null")] <- NA_character_
  nums <- suppressWarnings(as.numeric(stripped))
  if (any(!is.na(stripped) & is.na(nums))) return(x)
  nums
}

read_qdec <- function(path) {
  q <- tryCatch(read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(q)) q <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  if ("fsid-base" %in% names(q)) names(q)[names(q) == "fsid-base"] <- "fsid_base"
  if ("fsid.base" %in% names(q)) names(q)[names(q) == "fsid.base"] <- "fsid_base"
  for (nm in names(q)) q[[nm]] <- coerce_numeric_like(q[[nm]])
  q
}

read_aseg <- function(path, id_col_hint = NULL) {
  aseg <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  id_col <- if (!is.null(id_col_hint) && id_col_hint %in% names(aseg)) {
    id_col_hint
  } else if ("Measure.volume" %in% names(aseg)) {
    "Measure.volume"
  } else {
    "Measure:volume"
  }
  ids <- as.character(aseg[[id_col]])
  has_long <- grepl("\\.long\\.", ids)
  aseg$fsid <- ifelse(has_long, sub("^(.*)\\.long\\..*$", "\\1", ids, perl = TRUE), ids)
  aseg$fsid_base <- ifelse(has_long, sub("^.*\\.long\\.(.*)$", "\\1", ids, perl = TRUE), NA_character_)
  aseg
}

qdec <- read_qdec(cfg$qdec)
aseg <- read_aseg(cfg$aseg, cfg$id_col)
if (!all(c("fsid", "fsid_base") %in% names(qdec))) stop("qdec must contain fsid and fsid_base")
if (!all(c("fsid", "fsid_base") %in% names(aseg))) stop("aseg must contain derived fsid and fsid_base")

dat <- merge(qdec, aseg, by = c("fsid", "fsid_base"))
dat$sex <- factor(dat$sex)
if (!("tp" %in% names(dat))) stop("Merged data lacks 'tp' column")
dat$tp <- suppressWarnings(as.numeric(dat$tp))

coef_csv <- file.path(results_dir, "lme_coefficients.csv")
if (!file.exists(coef_csv)) stop("lme_coefficients.csv not found in results directory")
coef_df <- read.csv(coef_csv, stringsAsFactors = FALSE)
if (!("roi" %in% names(coef_df))) stop("lme_coefficients.csv missing 'roi' column")

effect_rows <- subset(coef_df, coef == opt$effect)
if (!nrow(effect_rows)) stop(sprintf("Effect '%s' not found in lme_coefficients.csv", opt$effect))

roi_to_filename <- function(roi) gsub("[^A-Za-z0-9_-]", "_", roi)
fit_dir <- file.path(results_dir, "individual_rois")
if (!dir.exists(fit_dir)) stop("individual_rois directory missing from results")

results <- list()
messages <- character(0)

for (i in seq_len(nrow(effect_rows))) {
  roi <- effect_rows$roi[i]
  beta <- effect_rows$beta[i]
  fit_path <- file.path(fit_dir, sprintf("%s_fit.rds", roi_to_filename(roi)))
  if (!file.exists(fit_path)) {
    messages <- c(messages, sprintf("[WARN] fit file missing for %s", roi))
    next
  }
  fit <- readRDS(fit_path)
  coef_names <- colnames(fit$X)
  if (is.null(coef_names) || !(opt$effect %in% coef_names)) {
    messages <- c(messages, sprintf("[WARN] Effect '%s' not present in design for %s", opt$effect, roi))
    next
  }
  idx <- match(opt$effect, coef_names)
  se <- sqrt(fit$CovBhat[idx, idx])
  z <- if (is.na(se) || se == 0) NA_real_ else beta / se
  p <- if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z))
  category <- if (is.na(p)) {
    "unknown"
  } else if (p < opt$alpha) {
    "significant"
  } else if (p < opt$`trend-alpha`) {
    "trend"
  } else {
    "non-significant"
  }
  # Some effects (e.g., derived design columns) are not ROI columns in merged data; handle gracefully
  if (roi %in% names(dat)) {
    roi_data <- dat[c("fsid", "fsid_base", "tp", "sex", roi)]
    names(roi_data)[names(roi_data) == roi] <- "value"
    roi_data <- roi_data[!is.na(roi_data$value) & !is.na(roi_data$tp), ]
    n_subj <- length(unique(roi_data$fsid_base))
  } else {
    roi_data <- data.frame()
    n_subj <- length(unique(dat$fsid_base))
  }

  plot_path <- NA_character_
  if (category == "significant" && nrow(roi_data) > 0) {
    plot_obj <- ggplot(roi_data, aes(x = tp, y = value, color = sex, group = interaction(fsid_base, sex))) +
      geom_line(alpha = 0.25) +
      geom_point(alpha = 0.4, size = 1.5) +
      stat_summary(aes(group = sex), fun = mean, geom = "line", linewidth = 1.1) +
      stat_summary(aes(group = sex), fun = mean, geom = "point", size = 2.2, shape = 21, fill = "white", color = NA) +
      labs(title = sprintf("%s — %s", roi, opt$effect),
           x = "Time point (tp)",
           y = roi,
           color = "Sex") +
      theme_minimal() +
      theme(legend.position = "bottom")
    plot_file <- file.path(plot_dir, sprintf("%s.png", roi_to_filename(roi)))
    ggsave(plot_file, plot = plot_obj, width = 6.5, height = 4.5, dpi = 150)
    plot_path <- normalizePath(plot_file)
  }

  results[[roi]] <- data.frame(
    roi = roi,
    beta = beta,
    se = se,
    z = z,
    p = p,
    category = category,
    n_subjects = n_subj,
    plot = plot_path,
    stringsAsFactors = FALSE
  )
}

if (length(messages) && opt$verbose) cat(paste(messages, collapse = "\n"), "\n", file = stderr())

summary_df <- if (length(results)) do.call(rbind, results) else data.frame()
if (!nrow(summary_df)) stop("No ROI summaries computed")
summary_df <- summary_df[order(summary_df$p), ]
effect_slug <- gsub("[^A-Za-z0-9]+", "_", opt$effect)
summary_csv <- file.path(results_dir, sprintf("summary_%s.csv", effect_slug))
write.csv(summary_df, summary_csv, row.names = FALSE)

rel_path <- function(path, root) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  path_norm <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (startsWith(path_norm, paste0(root_norm, "/"))) {
    substring(path_norm, nchar(root_norm) + 2L)
  } else {
    path_norm
  }
}

html_dir <- dirname(html_path)
html_lines <- c(
  "<!DOCTYPE html>",
  "<html>",
  "<head>",
  "  <meta charset='utf-8'>",
  sprintf("  <title>fslmer summary — %s</title>", opt$effect),
  "  <style>",
  "body { font-family: Arial, sans-serif; margin: 1.5rem; background: #f9fafc; color: #222; }",
  "h1 { margin-bottom: 0.2rem; }",
  "h2 { margin-top: 2rem; }",
  "table { border-collapse: collapse; margin-top: 0.5rem; width: 100%; max-width: 900px; background: #fff; }",
  "th, td { border: 1px solid #d8dee4; padding: 0.5rem 0.7rem; text-align: left; }",
  "th { background: #f0f3f8; }",
  "tr:nth-child(even) { background: #f8f9fc; }",
  ".plot { margin: 0.7rem 0 1.4rem 0; }",
  ".plot img { max-width: 680px; height: auto; border: 1px solid #d8dee4; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }",
  ".meta { font-size: 0.9rem; color: #555; margin-bottom: 1rem; }",
  "  </style>",
  "</head>",
  "<body>",
  sprintf("<h1>fslmer summary — effect %s</h1>", opt$effect),
  sprintf("<p class='meta'>Results directory: %s<br>Summary CSV: %s<br>Significant if p &lt; %.2f; trend if p &lt; %.2f.</p>",
          results_dir, basename(summary_csv), opt$alpha, opt$`trend-alpha`)
)

if (length(messages)) {
  html_lines <- c(html_lines, "<div class='meta'><strong>Warnings:</strong><ul>",
                  paste0("<li>", messages, "</li>"), "</ul></div>")
}

category_order <- c("significant", "trend", "non-significant", "unknown")
for (cat in category_order) {
  df_cat <- subset(summary_df, category == cat)
  if (!nrow(df_cat)) next
  html_lines <- c(html_lines, sprintf("<h2>%s (%d)</h2>", tools::toTitleCase(cat), nrow(df_cat)))
  html_lines <- c(html_lines, "<table>",
                  "<tr><th>ROI</th><th>Beta</th><th>SE</th><th>Z</th><th>p-value</th><th>Subjects</th></tr>")
  for (j in seq_len(nrow(df_cat))) {
    row <- df_cat[j, ]
    html_lines <- c(html_lines,
      sprintf("<tr><td>%s</td><td>%.4f</td><td>%s</td><td>%s</td><td>%s</td><td>%d</td></tr>",
              row$roi,
              row$beta,
              ifelse(is.na(row$se), "NA", sprintf("%.4f", row$se)),
              ifelse(is.na(row$z), "NA", sprintf("%.3f", row$z)),
              ifelse(is.na(row$p), "NA", sprintf("%.4g", row$p)),
              row$n_subjects))
  }
  html_lines <- c(html_lines, "</table>")
  if (cat == "significant") {
    for (j in seq_len(nrow(df_cat))) {
      row <- df_cat[j, ]
      if (is.na(row$plot) || !nzchar(row$plot)) next
      rel_plot <- rel_path(row$plot, html_dir)
      html_lines <- c(html_lines,
        sprintf("<div class='plot'><h3>%s</h3><img src='%s' alt='%s plot'></div>",
                row$roi, rel_plot, row$roi))
    }
  }
}

html_lines <- c(html_lines, "</body>", "</html>")
dir.create(dirname(html_path), recursive = TRUE, showWarnings = FALSE)
writeLines(html_lines, html_path)

if (opt$verbose) {
  cat(sprintf("Summary CSV: %s\n", summary_csv))
  cat(sprintf("HTML report: %s\n", html_path))
}
