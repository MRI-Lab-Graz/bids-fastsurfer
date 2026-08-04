#!/usr/bin/env Rscript
#
# QC for longitudinal hippocampal subfield THICKNESS (hipsta) prior to
# modelling. Consumes the subfield-level tidy long TSV produced by
# scripts/extract_hipsta_thickness.py (subject_id, session, hemisphere,
# region, value, median_value, sd_value, n_vertices). Mirrors
# scripts/qc_subfields.R's three checks, retargeted at thickness; unlike
# volume, thickness has no head/body split to sum first (the extractor
# already emits one row per subject/session/hemisphere/subfield).
#
# Checks:
#   1. Implausible thickness (<= --min-thickness or >= --max-thickness mm --
#      hipsta's grid can degenerate near cutting/parametrization boundaries)
#   2. Within-subject temporal instability per subfield/hemisphere
#      (coefficient of variation across sessions; flagged via MAD z-score of CV)
#   3. Session completeness per subject
#   4. Grid coverage (n_vertices) -- a subfield with an unusually low vertex
#      count for a given subject/session/hemisphere means its mean thickness
#      is an average over very few grid points and should be treated with
#      more caution than a well-covered one; this has no equivalent in
#      qc_subfields.R since FreeSurfer's native subfield volumes don't carry
#      a comparable "how much of the subfield did this measurement span" signal.
#
# Exclusions must be decided from this output and frozen BEFORE any group-wise
# model is fit (see scripts/analysis/29_hippo_thickness_lmm.R).

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to hippo_thickness_tidy.tsv"),
  make_option(c("-o", "--outdir"), type="character", default="qc_hipsta_thickness", help="Output directory [default %default]"),
  make_option(c("--min-thickness"), type="double", default=0.3, help="Thickness <= this (mm) is flagged as implausibly low [default %default]"),
  make_option(c("--max-thickness"), type="double", default=6.0, help="Thickness >= this (mm) is flagged as implausibly high [default %default]"),
  make_option(c("--cv-z-threshold"), type="double", default=3.0, help="MAD z-score threshold for flagging high temporal CV [default %default]"),
  make_option(c("--expected-sessions"), type="integer", default=2, help="Expected number of sessions per subject [default %default]"),
  make_option(c("--min-vertices"), type="integer", default=20, help="n_vertices <= this is flagged as poor grid coverage for that subfield [default %default]"),
  make_option(c("--include-unlabeled"), action="store_true", default=FALSE,
              help="Include the 'unlabeled' grid-edge category in the implausibility/CV checks"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required (path to hippo_thickness_tidy.tsv)")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

dat <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_cols <- c("subject_id","session","hemisphere","region","value","n_vertices")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols)) stop(sprintf("tidy file missing required columns: %s", paste(missing_cols, collapse=", ")))

dat$value <- suppressWarnings(as.numeric(dat$value))
dat$n_vertices <- suppressWarnings(as.integer(dat$n_vertices))

analysis_dat <- if (isTRUE(opt$`include-unlabeled`)) dat else dat[dat$region != "unlabeled", , drop=FALSE]

mad_z <- function(x) {
  x <- as.numeric(x)
  m <- median(x, na.rm=TRUE)
  s <- mad(x, center=m, constant=1.4826, na.rm=TRUE)
  if (isTRUE(is.na(s)) || s == 0) {
    mu <- mean(x, na.rm=TRUE); sdv <- sd(x, na.rm=TRUE)
    if (isTRUE(is.na(sdv)) || sdv == 0) return(rep(0, length(x)))
    return((x - mu) / sdv)
  }
  (x - m) / s
}

# ----------------------------------------------------------------------------
# 1. Implausible thickness
# ----------------------------------------------------------------------------
msg("Checking for implausible thickness (<= %.2f or >= %.2f mm)...\n", opt$`min-thickness`, opt$`max-thickness`)
implausible <- analysis_dat[!is.na(analysis_dat$value) &
                             (analysis_dat$value <= opt$`min-thickness` | analysis_dat$value >= opt$`max-thickness`), , drop=FALSE]
implausible <- implausible[order(implausible$subject_id, implausible$session, implausible$region), ]
write.csv(implausible, file.path(opt$outdir, "implausible_thickness.csv"), row.names=FALSE)
msg("  -> %d implausible rows written to implausible_thickness.csv\n", nrow(implausible))

# ----------------------------------------------------------------------------
# 2. Within-subject temporal instability (CV per subject x subfield x hemisphere)
# ----------------------------------------------------------------------------
msg("Checking within-subject temporal stability...\n")
key <- interaction(analysis_dat$subject_id, analysis_dat$hemisphere, analysis_dat$region, drop=TRUE)
cv_list <- lapply(split(analysis_dat, key), function(g) {
  g <- g[!is.na(g$value), , drop=FALSE]
  if (nrow(g) < 2) return(NULL)
  data.frame(
    subject_id = g$subject_id[1],
    hemisphere = g$hemisphere[1],
    region = g$region[1],
    n_sessions = length(unique(g$session)),
    mean_thickness = mean(g$value),
    sd_thickness = sd(g$value),
    cv_pct = 100 * sd(g$value) / mean(g$value),
    stringsAsFactors = FALSE
  )
})
cv_df <- do.call(rbind, cv_list[!vapply(cv_list, is.null, logical(1))])

if (!is.null(cv_df) && nrow(cv_df)) {
  # MAD z-score computed within subfield (different subfields have different
  # baseline CV -- e.g. thin CA1 is noisier by nature than presubiculum)
  cv_df$cv_z <- ave(cv_df$cv_pct, cv_df$region, FUN = mad_z)
  cv_df$flagged <- abs(cv_df$cv_z) > opt$`cv-z-threshold`

  write.csv(cv_df[order(-cv_df$cv_pct), ], file.path(opt$outdir, "temporal_cv.csv"), row.names=FALSE)
  flagged_cv <- cv_df[cv_df$flagged, , drop=FALSE]
  write.csv(flagged_cv[order(-flagged_cv$cv_pct), ], file.path(opt$outdir, "temporal_cv_flagged.csv"), row.names=FALSE)
  msg("  -> %d subject/subfield/hemisphere combos flagged (|z| > %.1f) of %d total, written to temporal_cv_flagged.csv\n",
      nrow(flagged_cv), opt$`cv-z-threshold`, nrow(cv_df))
} else {
  cv_df <- data.frame()
  msg("  -> no subjects with >= 2 sessions found; skipping CV check\n")
}

# ----------------------------------------------------------------------------
# 3. Session completeness per subject
# ----------------------------------------------------------------------------
msg("Checking session completeness...\n")
sessions_per_subject <- aggregate(session ~ subject_id, data=unique(dat[, c("subject_id","session")]), FUN=length)
names(sessions_per_subject) <- c("subject_id", "n_sessions_present")
sessions_per_subject$expected_sessions <- opt$`expected-sessions`
sessions_per_subject$complete <- sessions_per_subject$n_sessions_present >= opt$`expected-sessions`
incomplete <- sessions_per_subject[!sessions_per_subject$complete, , drop=FALSE]
incomplete <- incomplete[order(incomplete$n_sessions_present), ]
write.csv(sessions_per_subject[order(sessions_per_subject$subject_id), ], file.path(opt$outdir, "session_completeness.csv"), row.names=FALSE)
write.csv(incomplete, file.path(opt$outdir, "session_completeness_incomplete.csv"), row.names=FALSE)
msg("  -> %d of %d subjects have < %d sessions, written to session_completeness_incomplete.csv\n",
    nrow(incomplete), nrow(sessions_per_subject), opt$`expected-sessions`)

# ----------------------------------------------------------------------------
# 4. Grid coverage
# ----------------------------------------------------------------------------
msg("Checking grid coverage (n_vertices per subfield)...\n")
low_coverage <- analysis_dat[!is.na(analysis_dat$n_vertices) & analysis_dat$n_vertices <= opt$`min-vertices`, , drop=FALSE]
low_coverage <- low_coverage[order(low_coverage$n_vertices), ]
write.csv(low_coverage, file.path(opt$outdir, "low_grid_coverage.csv"), row.names=FALSE)
msg("  -> %d subject/session/subfield/hemisphere rows with <= %d grid vertices, written to low_grid_coverage.csv\n",
    nrow(low_coverage), opt$`min-vertices`)

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf(
  "QC hipsta thickness summary\n\nInput: %s\n\n1. Implausible thickness (<= %.2f or >= %.2f mm): %d rows\n2. Temporal instability (|CV z-score| > %.1f): %d subject/subfield/hemisphere combos\n3. Incomplete subjects (< %d sessions): %d of %d subjects\n4. Low grid coverage (<= %d vertices): %d rows\n\nSee implausible_thickness.csv, temporal_cv_flagged.csv, session_completeness_incomplete.csv,\nlow_grid_coverage.csv for details.\nExclusion decisions based on this QC must be frozen before fitting any group-wise model.\n",
  opt$tidy, opt$`min-thickness`, opt$`max-thickness`, nrow(implausible), opt$`cv-z-threshold`,
  if (!is.null(cv_df) && nrow(cv_df)) sum(cv_df$flagged) else 0L,
  opt$`expected-sessions`, nrow(incomplete), nrow(sessions_per_subject),
  opt$`min-vertices`, nrow(low_coverage)
), file=con)
msg("\nDone. See %s/summary.txt\n", opt$outdir)
