#!/usr/bin/env Rscript
#
# QC for longitudinal hippocampal subfield volumes prior to modelling.
# Consumes the tidy long TSV produced by scripts/extract_hippo_subfields.py
# (subject_id, session, hemisphere, subfield, region, is_composite, volume, etiv).
#
# Checks (mirrors the MAD-outlier style of qc_etiv_outliers.R, but scoped to
# subfield-specific issues that script doesn't cover):
#   1. Implausible volumes (<= --min-volume, e.g. failed/degenerate segmentation)
#   2. Within-subject temporal instability per subfield/hemisphere (coefficient
#      of variation across sessions; flagged via MAD z-score of CV)
#   3. Session completeness per subject (missing timepoints)
#
# Exclusions must be decided from this output and frozen BEFORE any group-wise
# model is fit (see scripts/analysis/01_primary_lmm.R).

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to hippo_subfields_tidy.tsv"),
  make_option(c("-o", "--outdir"), type="character", default="qc_subfields", help="Output directory [default %default]"),
  make_option(c("--min-volume"), type="double", default=1.0, help="Volumes <= this (mm^3) are flagged as implausible [default %default]"),
  make_option(c("--cv-z-threshold"), type="double", default=3.0, help="MAD z-score threshold for flagging high temporal CV [default %default]"),
  make_option(c("--expected-sessions"), type="integer", default=3, help="Expected number of sessions per subject [default %default]"),
  make_option(c("--include-composite"), action="store_true", default=FALSE, help="Include composite/summary measures (Whole_hippocampus etc.) in CV/implausibility checks"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required (path to hippo_subfields_tidy.tsv)")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

dat <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_cols <- c("subject_id","session","hemisphere","subfield","region","is_composite","volume")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols)) stop(sprintf("tidy file missing required columns: %s", paste(missing_cols, collapse=", ")))

dat$is_composite <- as.logical(dat$is_composite)
dat$volume <- suppressWarnings(as.numeric(dat$volume))

analysis_dat <- if (isTRUE(opt$`include-composite`)) dat else dat[!dat$is_composite, , drop=FALSE]

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
# 1. Implausible volumes
# ----------------------------------------------------------------------------
msg("Checking for implausible volumes (<= %.2f mm^3)...\n", opt$`min-volume`)
implausible <- analysis_dat[!is.na(analysis_dat$volume) & analysis_dat$volume <= opt$`min-volume`, , drop=FALSE]
implausible <- implausible[order(implausible$subject_id, implausible$session, implausible$subfield), ]
write.csv(implausible, file.path(opt$outdir, "implausible_volumes.csv"), row.names=FALSE)
msg("  -> %d implausible rows written to implausible_volumes.csv\n", nrow(implausible))

# ----------------------------------------------------------------------------
# 2. Within-subject temporal instability (CV per subject x subfield x hemisphere)
# ----------------------------------------------------------------------------
# Subfield identity (e.g. "CA1") spans separate head/body rows per session
# (region column). Those must be summed into one subfield-total per session
# BEFORE computing temporal CV -- otherwise head and body volumes (very
# different absolute scales) get mixed together as if they were repeated
# measurements of the same quantity across time, producing a meaningless,
# inflated CV.
msg("Checking within-subject temporal stability...\n")
session_totals <- aggregate(volume ~ subject_id + session + hemisphere + subfield,
                             data=analysis_dat, FUN=sum)

key <- interaction(session_totals$subject_id, session_totals$hemisphere, session_totals$subfield, drop=TRUE)
cv_list <- lapply(split(session_totals, key), function(g) {
  g <- g[!is.na(g$volume), , drop=FALSE]
  if (nrow(g) < 2) return(NULL)
  data.frame(
    subject_id = g$subject_id[1],
    hemisphere = g$hemisphere[1],
    subfield = g$subfield[1],
    n_sessions = length(unique(g$session)),
    mean_volume = mean(g$volume),
    sd_volume = sd(g$volume),
    cv_pct = 100 * sd(g$volume) / mean(g$volume),
    stringsAsFactors = FALSE
  )
})
cv_df <- do.call(rbind, cv_list[!vapply(cv_list, is.null, logical(1))])

if (!is.null(cv_df) && nrow(cv_df)) {
  # MAD z-score computed within subfield (different subfields have very
  # different baseline CV, e.g. tiny structures are noisier by nature)
  cv_df$cv_z <- ave(cv_df$cv_pct, cv_df$subfield, FUN = mad_z)
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
# Summary
# ----------------------------------------------------------------------------
con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf(
  "QC subfields summary\n\nInput: %s\n\n1. Implausible volumes (<= %.2f mm^3): %d rows\n2. Temporal instability (|CV z-score| > %.1f): %d subject/subfield/hemisphere combos\n3. Incomplete subjects (< %d sessions): %d of %d subjects\n\nSee implausible_volumes.csv, temporal_cv_flagged.csv, session_completeness_incomplete.csv for details.\nExclusion decisions based on this QC must be frozen before fitting any group-wise model.\n",
  opt$tidy, opt$`min-volume`, nrow(implausible), opt$`cv-z-threshold`,
  if (!is.null(cv_df) && nrow(cv_df)) sum(cv_df$flagged) else 0L,
  opt$`expected-sessions`, nrow(incomplete), nrow(sessions_per_subject)
), file=con)

msg("QC done. Outputs in %s\n", opt$outdir)
