#!/usr/bin/env Rscript
#
# Responder heterogeneity analysis: does a null GROUP-MEAN effect hide real,
# but individually heterogeneous, responses -- e.g. some participants
# showing an early rise that fades by the final session, others a later or
# more sustained change, cancelling out in a mean-based test?
#
# IMPORTANT DESIGN NOTE: this deliberately does NOT classify subjects into
# "early responder" / "late responder" groups from their own change scores
# and then test those groups -- that's circular (the classification and the
# test would use the same information, guaranteeing an inflated effect).
# Instead it asks two answerable, non-circular questions:
#
#   1. Variance heterogeneity: is the SPREAD of early change (ses1->ses2)
#      and late change (ses2->ses3) larger in dance than control (Levene's
#      test)? Larger variance with a null mean is exactly what you'd expect
#      if individuals respond in different directions/magnitudes/timing.
#
#   2. Early-late correlation: within each group, does early change
#      correlate with late change? A negative correlation in dance but not
#      control (early gain predicting later fade) would be a genuine
#      transient-responder signature; compared between groups via Fisher's
#      z-test for a difference in correlations.
#
# Both are legitimate hypothesis tests using the full sample, not a
# post-hoc subgroup comparison -- but still exploratory in the sense that
# it wasn't part of the original pre-specified analysis plan.

suppressPackageStartupMessages({
  library(optparse)
  library(car)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/12_responder_heterogeneity",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load data and build per-subject, bilateral-average early/late change
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

# subject x session x subfield (bilateral total: sum across hemisphere and
# across head/body region, consistent with the primary confirmatory ROI
# definition used throughout this pipeline)
agg <- aggregate(volume ~ subject_id + session + subfield, data=roi_dat, FUN=sum)
agg$log_volume <- log(agg$volume)

sessions <- sort(unique(agg$session))
if (length(sessions) < 3) stop("Need at least 3 sessions to compute early (ses1->ses2) and late (ses2->ses3) change")
ses1 <- sessions[1]; ses2 <- sessions[2]; ses3 <- sessions[3]

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

results_rows <- list()
wide_by_roi <- list()

for (roi_name in roi_set) {
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  wide <- reshape(d_roi[, c("subject_id","session","log_volume")], idvar="subject_id", timevar="session", direction="wide")
  names(wide) <- sub("^log_volume\\.", "", names(wide))
  needed <- c("subject_id", ses1, ses2, ses3)
  if (!all(needed %in% names(wide))) next
  wide <- wide[complete.cases(wide[, c(ses1, ses2, ses3)]), needed]
  names(wide) <- c("subject_id", "ses1", "ses2", "ses3")

  wide$early_change <- wide$ses2 - wide$ses1
  wide$late_change <- wide$ses3 - wide$ses2

  wide <- merge(wide, participants[, required_pcols], by="subject_id")
  wide <- wide[wide$group %in% c(dance_groups, control_group), , drop=FALSE]
  wide$dance <- factor(ifelse(wide$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))

  if (sum(wide$dance == "dance") < 5 || sum(wide$dance == "control") < 5) next
  wide_by_roi[[roi_name]] <- wide

  # --- 1. Variance heterogeneity (Levene's test, dance vs control) ---
  levene_early <- tryCatch(car::leveneTest(early_change ~ dance, data=wide), error=function(e) NULL)
  levene_late <- tryCatch(car::leveneTest(late_change ~ dance, data=wide), error=function(e) NULL)

  sd_dance_early <- sd(wide$early_change[wide$dance == "dance"])
  sd_control_early <- sd(wide$early_change[wide$dance == "control"])
  sd_dance_late <- sd(wide$late_change[wide$dance == "dance"])
  sd_control_late <- sd(wide$late_change[wide$dance == "control"])

  # --- 2. Early-late correlation, compared between groups (Fisher's z) ---
  cor_dance <- cor.test(wide$early_change[wide$dance == "dance"], wide$late_change[wide$dance == "dance"])
  cor_control <- cor.test(wide$early_change[wide$dance == "control"], wide$late_change[wide$dance == "control"])

  n_dance <- sum(wide$dance == "dance"); n_control <- sum(wide$dance == "control")
  z_dance <- atanh(cor_dance$estimate); z_control <- atanh(cor_control$estimate)
  se_diff <- sqrt(1/(n_dance - 3) + 1/(n_control - 3))
  z_stat <- (z_dance - z_control) / se_diff
  fisher_p <- 2 * pnorm(-abs(z_stat))

  results_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    n_dance = n_dance, n_control = n_control,
    sd_early_dance = sd_dance_early, sd_early_control = sd_control_early,
    levene_early_p = if (!is.null(levene_early)) levene_early[["Pr(>F)"]][1] else NA_real_,
    sd_late_dance = sd_dance_late, sd_late_control = sd_control_late,
    levene_late_p = if (!is.null(levene_late)) levene_late[["Pr(>F)"]][1] else NA_real_,
    cor_early_late_dance = as.numeric(cor_dance$estimate),
    cor_early_late_control = as.numeric(cor_control$estimate),
    fisher_z_diff_p = fisher_p,
    stringsAsFactors = FALSE
  )
}

results_df <- do.call(rbind, results_rows)
if (!is.null(results_df) && nrow(results_df)) {
  results_df$levene_early_p_fdr <- p.adjust(results_df$levene_early_p, method="fdr")
  results_df$levene_late_p_fdr <- p.adjust(results_df$levene_late_p, method="fdr")
  results_df$fisher_z_diff_p_fdr <- p.adjust(results_df$fisher_z_diff_p, method="fdr")
  write.csv(results_df, file.path(opt$outdir, "responder_heterogeneity_summary.csv"), row.names=FALSE)

  con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat(sprintf("Responder heterogeneity analysis (%s -> %s -> %s)\n\n", ses1, ses2, ses3), file=con)
  cat("1. Variance heterogeneity (Levene's test): is early/late change MORE VARIABLE in dance than control?\n", file=con)
  cat("   (larger variance with a null mean = consistent with heterogeneous individual responses)\n\n", file=con)
  for (i in seq_len(nrow(results_df))) {
    row <- results_df[i, ]
    cat(sprintf("   %s: early SD dance=%.4f vs control=%.4f (Levene p=%.3f, p_fdr=%.3f)%s | late SD dance=%.4f vs control=%.4f (Levene p=%.3f, p_fdr=%.3f)%s\n",
                row$roi, row$sd_early_dance, row$sd_early_control, row$levene_early_p, row$levene_early_p_fdr,
                if (row$levene_early_p_fdr < opt$alpha) " *" else "",
                row$sd_late_dance, row$sd_late_control, row$levene_late_p, row$levene_late_p_fdr,
                if (row$levene_late_p_fdr < opt$alpha) " *" else ""), file=con)
  }
  cat("\n2. Early-late correlation (transient-responder signature: negative correlation\n", file=con)
  cat("   in dance -- early gain predicting later fade -- differing from control)\n\n", file=con)
  for (i in seq_len(nrow(results_df))) {
    row <- results_df[i, ]
    cat(sprintf("   %s: cor(early,late) dance=%.3f vs control=%.3f (Fisher z-test p=%.3f, p_fdr=%.3f)%s\n",
                row$roi, row$cor_early_late_dance, row$cor_early_late_control, row$fisher_z_diff_p, row$fisher_z_diff_p_fdr,
                if (row$fisher_z_diff_p_fdr < opt$alpha) " *" else ""), file=con)
  }
  msg("Done. See %s/summary.txt\n", opt$outdir)
} else {
  warning("No ROI had sufficient complete-case data for this analysis")
}

# ------------------------------------------------------------------------
# LOSO robustness for the early-late correlation (Fisher's z) test -- the
# key question being "is this one or two subjects, or a real group-level
# pattern?" Leave-one-subject-out, per group (dance and control separately,
# since a real group-specific effect should not depend on any single
# subject in EITHER group), recompute the Fisher z-test p-value each time.
# ------------------------------------------------------------------------
if (!is.null(results_df) && nrow(results_df) && length(wide_by_roi)) {
  msg("\nRunning LOSO for the early-late correlation test...\n")
  loso_rows <- list()

  fisher_p_for <- function(w) {
    dance_vals <- w[w$dance == "dance", ]
    control_vals <- w[w$dance == "control", ]
    if (nrow(dance_vals) < 5 || nrow(control_vals) < 5) return(list(p = NA_real_, cor_dance = NA_real_, cor_control = NA_real_))
    cd <- suppressWarnings(cor.test(dance_vals$early_change, dance_vals$late_change))
    cc <- suppressWarnings(cor.test(control_vals$early_change, control_vals$late_change))
    nd <- nrow(dance_vals); nc <- nrow(control_vals)
    zd <- atanh(cd$estimate); zc <- atanh(cc$estimate)
    se <- sqrt(1/(nd - 3) + 1/(nc - 3))
    z_stat <- (zd - zc) / se
    list(p = 2 * pnorm(-abs(z_stat)), cor_dance = as.numeric(cd$estimate), cor_control = as.numeric(cc$estimate))
  }

  for (roi_name in names(wide_by_roi)) {
    w <- wide_by_roi[[roi_name]]
    full <- fisher_p_for(w)
    for (excl_subj in w$subject_id) {
      w_sub <- w[w$subject_id != excl_subj, , drop=FALSE]
      res <- fisher_p_for(w_sub)
      loso_rows[[length(loso_rows) + 1]] <- data.frame(
        roi = roi_name, excluded_subject = excl_subj,
        loso_p = res$p, full_p = full$p,
        stringsAsFactors = FALSE
      )
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  write.csv(loso_df, file.path(opt$outdir, "loso_fisher_z.csv"), row.names=FALSE)

  flips <- loso_df[!is.na(loso_df$loso_p) & !is.na(loso_df$full_p) &
                    (loso_df$loso_p < opt$alpha) != (loso_df$full_p < opt$alpha), , drop=FALSE]
  write.csv(flips, file.path(opt$outdir, "loso_fisher_z_flips.csv"), row.names=FALSE)

  con2 <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con2), add=TRUE)
  cat("LOSO robustness for the early-late correlation (Fisher's z) test\n\n", file=con2)
  for (roi_name in names(wide_by_roi)) {
    g <- loso_df[loso_df$roi == roi_name & !is.na(loso_df$loso_p), , drop=FALSE]
    if (!nrow(g)) next
    n_flips_roi <- sum((g$loso_p < opt$alpha) != (g$full_p[1] < opt$alpha))
    most_infl_idx <- which.max(abs(g$loso_p - g$full_p[1]))
    cat(sprintf("%s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance, most influential = %s (p becomes %.4f)\n",
                roi_name, g$full_p[1], min(g$loso_p), max(g$loso_p), n_flips_roi, nrow(g),
                g$excluded_subject[most_infl_idx], g$loso_p[most_infl_idx]), file=con2)
  }
  msg("LOSO done. See %s/loso_summary.txt\n", opt$outdir)
}
