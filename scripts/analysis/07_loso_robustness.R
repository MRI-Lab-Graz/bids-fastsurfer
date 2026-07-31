#!/usr/bin/env Rscript
#
# Leave-one-subject-out (LOSO) robustness check: refit the baseline-corrected
# change-score model (05) and the MANOVA joint pattern test (06) once per
# subject, excluding that subject each time. Reveals whether the (current
# null) result is fragile -- i.e. whether excluding any single subject would
# flip a conclusion, or whether a handful of subjects are pulling estimates
# around disproportionately (high leverage / possible data-quality issues).
#
# This does NOT replace QC (scripts/qc_subfields.R) -- it complements it.
# QC flags individual noisy measurements; LOSO flags subjects whose overall
# presence/absence changes the group-level conclusion, which QC alone
# cannot tell you (a flagged subject might still have negligible influence
# on the model, and an unflagged one could still be an influential outlier
# in the group contrast specifically).

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/07_loso_robustness",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--intervention-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'intervention' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value used as baseline/covariate [default: earliest session present]"),
  make_option(c("--final-session"), type="character", default=NULL,
              help="Session value used as endpoint for the MANOVA change score [default: latest session present]"),
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
intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load and prepare data (same construction as 05_baseline_corrected_change.R)
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=roi_dat, FUN=sum)
# eTIV varies slightly session-to-session (FreeSurfer re-estimation noise,
# not real anatomical change) -- use each subject's BASELINE (earliest
# session) eTIV as a fixed per-subject covariate, explicitly selected (not
# relying on incidental row order) so it can't silently pick up a different
# session's value.
baseline_ses_for_etiv <- sort(unique(tidy$session))[1]
etiv_lookup <- unique(tidy[tidy$session == baseline_ses_for_etiv, c("subject_id","etiv")])
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(agg, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(intervention_groups, control_group), , drop=FALSE]
agg$intervention <- factor(ifelse(agg$group %in% intervention_groups, "intervention", "control"), levels=c("control","intervention"))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
followup_sessions <- setdiff(sessions, baseline_ses)

build_change_data <- function(d_roi) {
  base <- d_roi[d_roi$session == baseline_ses, c("subject_id","hemisphere","log_volume")]
  names(base)[3] <- "log_baseline"
  fu <- d_roi[d_roi$session %in% followup_sessions, , drop=FALSE]
  merged <- merge(fu, base, by=c("subject_id","hemisphere"))
  merged$log_change <- merged$log_volume - merged$log_baseline
  merged$log_baseline_z <- as.numeric(scale(merged$log_baseline))
  merged$followup_f <- factor(merged$session, levels=followup_sessions)
  merged
}

all_subjects <- unique(agg$subject_id)
msg("Subjects available for LOSO: %d\n", length(all_subjects))

# ------------------------------------------------------------------------
# LOSO for the per-ROI baseline-corrected change-score model (05)
# intervention_main effect (the clean intervention test)
# ------------------------------------------------------------------------
# With only one follow-up session (the 2-timepoint design), followup_f is a
# constant and the intervention:followup_f / followup_f terms would be
# aliased with the intercept -- drop them rather than feeding lmer a
# rank-deficient design.
base_formula <- if (length(followup_sessions) > 1) {
  log_change ~ intervention * followup_f + log_baseline_z + hemisphere + age_z + sex + etiv_z + (1 | subject_id)
} else {
  log_change ~ intervention + log_baseline_z + hemisphere + age_z + sex + etiv_z + (1 | subject_id)
}

msg("Running LOSO for per-ROI change-score models (%d ROIs x %d subjects)...\n", length(roi_set), length(all_subjects))
loso_rows <- list()

for (roi_name in unique(agg$subfield)) {
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat_full <- build_change_data(d_roi)

  full_fit <- tryCatch(lmerTest::lmer(base_formula, data=change_dat_full, REML=TRUE), error=function(e) NULL)
  full_p <- if (!is.null(full_fit)) {
    at <- tryCatch(anova(full_fit), error=function(e) NULL)
    if (!is.null(at) && "intervention" %in% rownames(at)) at["intervention", "Pr(>F)"] else NA_real_
  } else NA_real_

  for (excl_subj in all_subjects) {
    d_sub <- change_dat_full[change_dat_full$subject_id != excl_subj, , drop=FALSE]
    if (!length(unique(d_sub$intervention)) == 2 || nrow(d_sub) < 10) next
    fit <- tryCatch(lmerTest::lmer(base_formula, data=d_sub, REML=TRUE), error=function(e) NULL)
    p_val <- NA_real_
    est <- NA_real_
    if (!is.null(fit)) {
      at <- tryCatch(anova(fit), error=function(e) NULL)
      if (!is.null(at) && "intervention" %in% rownames(at)) p_val <- at["intervention", "Pr(>F)"]
      fe <- tryCatch(lme4::fixef(fit), error=function(e) NULL)
      if (!is.null(fe) && "interventionintervention" %in% names(fe)) est <- fe[["interventionintervention"]]
    }
    loso_rows[[length(loso_rows) + 1]] <- data.frame(
      roi = roi_name, excluded_subject = excl_subj,
      intervention_main_p = p_val, intervention_main_estimate = est,
      full_data_p = full_p,
      stringsAsFactors = FALSE
    )
  }
}
loso_df <- do.call(rbind, loso_rows)
write.csv(loso_df, file.path(opt$outdir, "loso_change_score_by_roi.csv"), row.names=FALSE)

# Flag subjects whose exclusion flips significance, or causes the largest
# swing in the p-value / estimate relative to the full-data fit.
flip_rows <- loso_df[!is.na(loso_df$intervention_main_p) & !is.na(loso_df$full_data_p) &
                      (loso_df$intervention_main_p < opt$alpha) != (loso_df$full_data_p < opt$alpha), , drop=FALSE]
write.csv(flip_rows, file.path(opt$outdir, "loso_significance_flips.csv"), row.names=FALSE)
msg("LOSO change-score models: %d (ROI, excluded-subject) combinations flip significance status\n", nrow(flip_rows))

influence_summary <- do.call(rbind, lapply(split(loso_df, loso_df$roi), function(g) {
  g <- g[!is.na(g$intervention_main_p), , drop=FALSE]
  if (!nrow(g)) return(NULL)
  most_influential_idx <- which.max(abs(g$intervention_main_p - g$full_data_p[1]))
  data.frame(
    roi = g$roi[1],
    full_data_p = g$full_data_p[1],
    loso_p_min = min(g$intervention_main_p), loso_p_max = max(g$intervention_main_p),
    most_influential_subject = g$excluded_subject[most_influential_idx],
    most_influential_p_when_excluded = g$intervention_main_p[most_influential_idx],
    stringsAsFactors = FALSE
  )
}))
write.csv(influence_summary, file.path(opt$outdir, "loso_influence_summary.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# LOSO for the MANOVA joint pattern test (06)
# ------------------------------------------------------------------------
msg("Running LOSO for the MANOVA joint pattern test...\n")

group_cols <- c("subject_id","session","subfield")
agg_sum <- aggregate(volume ~ ., data=roi_dat[, c(group_cols, "volume")], FUN=sum)
base_vals <- agg_sum[agg_sum$session == baseline_ses, c("subject_id","subfield","volume")]
final_vals <- agg_sum[agg_sum$session == final_ses, c("subject_id","subfield","volume")]
names(base_vals)[3] <- "baseline"; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id","subfield"))
merged_long$change <- log(merged_long$final) - log(merged_long$baseline)
change_wide <- reshape(merged_long[, c("subject_id","subfield","change")], idvar="subject_id", timevar="subfield", direction="wide")
change_cols <- setdiff(names(change_wide), "subject_id")
names(change_wide)[names(change_wide) %in% change_cols] <- sub("^change\\.", "", change_cols)
change_cols <- setdiff(names(change_wide), "subject_id")

manova_dat <- merge(change_wide, participants[, required_pcols], by="subject_id")
manova_dat <- manova_dat[manova_dat$group %in% c(intervention_groups, control_group), , drop=FALSE]
manova_dat$intervention <- factor(ifelse(manova_dat$group %in% intervention_groups, "intervention", "control"), levels=c("control","intervention"))
manova_dat$sex <- factor(manova_dat$sex)
manova_dat$age_z <- as.numeric(scale(manova_dat$age))
manova_dat <- manova_dat[complete.cases(manova_dat[, c(change_cols,"intervention","age_z","sex")]), , drop=FALSE]

run_manova_p <- function(d) {
  if (length(unique(d$intervention)) < 2 || nrow(d) < length(change_cols) + 5) return(NA_real_)
  Y <- as.matrix(d[, change_cols])
  fit <- tryCatch(manova(Y ~ intervention + age_z + sex, data=d), error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- tryCatch(as.data.frame(summary(fit, test="Pillai")$stats), error=function(e) NULL)
  if (is.null(s) || !"intervention" %in% rownames(s)) return(NA_real_)
  s["intervention", "Pr(>F)"]
}

full_manova_p <- run_manova_p(manova_dat)
manova_loso_rows <- list()
for (excl_subj in unique(manova_dat$subject_id)) {
  d_sub <- manova_dat[manova_dat$subject_id != excl_subj, , drop=FALSE]
  manova_loso_rows[[excl_subj]] <- data.frame(
    excluded_subject = excl_subj,
    manova_intervention_p = run_manova_p(d_sub),
    full_data_p = full_manova_p,
    stringsAsFactors = FALSE
  )
}
manova_loso_df <- do.call(rbind, manova_loso_rows)
write.csv(manova_loso_df, file.path(opt$outdir, "loso_manova.csv"), row.names=FALSE)

manova_flips <- manova_loso_df[!is.na(manova_loso_df$manova_intervention_p) &
                                (manova_loso_df$manova_intervention_p < opt$alpha) != (manova_loso_df$full_data_p < opt$alpha), , drop=FALSE]
write.csv(manova_flips, file.path(opt$outdir, "loso_manova_flips.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------
con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf(
  "LOSO robustness summary\n\n1. Per-ROI baseline-corrected change-score model:\n   %d of %d (ROI, excluded-subject) combinations flip significance status (alpha=%.2f)\n\n",
  nrow(flip_rows), nrow(loso_df), opt$alpha
), file=con)
if (!is.null(influence_summary)) {
  for (i in seq_len(nrow(influence_summary))) {
    row <- influence_summary[i, ]
    cat(sprintf("   - %s: full p=%.3f, LOSO range [%.3f, %.3f], most influential exclusion = %s (p becomes %.3f)\n",
                row$roi, row$full_data_p, row$loso_p_min, row$loso_p_max,
                row$most_influential_subject, row$most_influential_p_when_excluded), file=con)
  }
} else {
  cat("   No ROI change-score model produced a usable intervention_main_p (see loso_change_score_by_roi.csv).\n", file=con)
}
cat(sprintf(
  "\n2. MANOVA joint pattern test:\n   Full-data p = %.4f\n   LOSO range: [%.4f, %.4f]\n   %d of %d exclusions flip significance status\n",
  full_manova_p, min(manova_loso_df$manova_intervention_p, na.rm=TRUE), max(manova_loso_df$manova_intervention_p, na.rm=TRUE),
  nrow(manova_flips), nrow(manova_loso_df)
), file=con)
if (nrow(manova_flips)) {
  cat("   Subjects whose exclusion flips the MANOVA conclusion:", paste(manova_flips$excluded_subject, collapse=", "), "\n", file=con)
}

msg("Done. See %s/summary.txt\n", opt$outdir)
