#!/usr/bin/env Rscript
#
# EXPLORATORY coupled-change analysis: does the CHANGE in a psychological
# measure (e.g. depression, anxiety, social support -- ses1 -> final
# session) track the CHANGE in brain volume over the same window? And does
# that coupling differ between the intervention and control groups?
#
# This is a different question from 08_moderator_analysis.R (which asks
# whether a BASELINE trait predicts the volume-change effect). Here the
# psychological measure is itself changing, and we ask whether its
# trajectory is coupled to the brain's trajectory -- directly relevant to
# a study framed around "does exercise-linked brain change track symptom
# change" (e.g. this study's "running away from depression" framing).
#
# Two tests per ROI x moderator combination:
#   1. Overall coupling: brain_change ~ moderator_change + age_z + sex + etiv_z
#      (pooled across all subjects, regardless of group)
#   2. Group-dependent coupling: brain_change ~ moderator_change * dance + ...
#      (does the coupling differ between intervention and control?)
#
# LOSO robustness is run automatically, but only for hits that are at
# least nominally interesting (uncorrected p < --loso-screen-threshold),
# for the same compute-efficiency reason as 08_moderator_analysis.R.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to a tidy volumes TSV (from extract_hippo_subfields.py / extract_amygdala_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-m", "--moderators"), type="character",
              help="Path to a LONGITUDINAL moderators TSV with columns: subject_id, session, plus moderator columns"),
  make_option(c("--moderator-cols"), type="character", required=TRUE,
              help="Comma-separated column names in --moderators to test (e.g. 'ads_score,fsozu_mean,stai_total')"),
  make_option(c("-o", "--outdir"), type="character", default="results/14_coupled_change_moderators",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified ROI names [default %default]"),
  make_option(c("--roi-column"), type="character", default="subfield",
              help="Column name identifying the ROI in the tidy file (subfield/nucleus/roi) [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the intervention contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value used as baseline [default: earliest session present]"),
  make_option(c("--final-session"), type="character", default=NULL,
              help="Session value used as endpoint [default: latest session present]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10,
              help="Only run LOSO follow-up on hits with uncorrected p below this [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (is.null(opt$moderators)) stop("--moderators is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))
if (!file.exists(opt$moderators)) stop(sprintf("moderators file not found: %s", opt$moderators))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)
roi_col <- opt$`roi-column`
moderator_cols <- trimws(strsplit(opt$`moderator-cols`, ",")[[1]])

# ------------------------------------------------------------------------
# Load brain data, compute per-subject per-ROI change score (baseline ->
# final session, bilateral-summed as elsewhere in this pipeline)
# ------------------------------------------------------------------------
msg("Loading tidy brain volume data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file", roi_col))
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check ROI names against the tidy file")

brain_agg <- aggregate(as.formula(paste("volume ~ subject_id + session +", roi_col)), data=roi_dat, FUN=sum)
names(brain_agg)[names(brain_agg) == roi_col] <- "roi"

sessions <- sort(unique(brain_agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Brain change window: %s -> %s\n", baseline_ses, final_ses)

base_vol <- brain_agg[brain_agg$session == baseline_ses, c("subject_id","roi","volume")]
final_vol <- brain_agg[brain_agg$session == final_ses, c("subject_id","roi","volume")]
names(base_vol)[3] <- "baseline_volume"; names(final_vol)[3] <- "final_volume"
brain_change <- merge(base_vol, final_vol, by=c("subject_id","roi"))
brain_change$brain_change <- log(brain_change$final_volume) - log(brain_change$baseline_volume)

# eTIV varies slightly session-to-session (FreeSurfer re-estimation noise,
# not real anatomical change) -- use each subject's BASELINE (earliest
# session) eTIV as a fixed per-subject covariate. (This was previously
# documented in the module docstring above but never actually wired into
# the model formulas below -- fixed here.)
etiv_lookup <- unique(tidy[tidy$session == baseline_ses, c("subject_id","etiv")])
brain_change <- merge(brain_change, etiv_lookup, by="subject_id", all.x=TRUE)

# ------------------------------------------------------------------------
# Load longitudinal moderators, compute per-subject change score for each
# candidate psychological measure over the SAME window
# ------------------------------------------------------------------------
msg("Loading longitudinal moderators from %s...\n", opt$moderators)
mod_long <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
missing_mod_cols <- setdiff(c("subject_id","session", moderator_cols), names(mod_long))
if (length(missing_mod_cols)) stop(sprintf("--moderators missing required columns: %s", paste(missing_mod_cols, collapse=", ")))

mod_base <- mod_long[mod_long$session == baseline_ses, c("subject_id", moderator_cols)]
mod_final <- mod_long[mod_long$session == final_ses, c("subject_id", moderator_cols)]
names(mod_base)[-1] <- paste0(moderator_cols, "_baseline")
names(mod_final)[-1] <- paste0(moderator_cols, "_final")
mod_change <- merge(mod_base, mod_final, by="subject_id")
for (col in moderator_cols) {
  b <- suppressWarnings(as.numeric(mod_change[[paste0(col, "_baseline")]]))
  f <- suppressWarnings(as.numeric(mod_change[[paste0(col, "_final")]]))
  mod_change[[paste0(col, "_change")]] <- f - b
}

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

dat <- merge(brain_change, participants[, required_pcols], by="subject_id")
dat <- merge(dat, mod_change[, c("subject_id", paste0(moderator_cols, "_change"))], by="subject_id")
dat <- dat[dat$group %in% c(dance_groups, control_group), , drop=FALSE]
dat$dance <- factor(ifelse(dat$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
dat$sex <- factor(dat$sex)
# Rank-transformed (not raw) age -- see 01_primary_lmm.R for rationale: this
# cohort's age distribution has a sparse, unevenly-populated tail, and rank
# transformation caps its leverage without discarding subjects.
dat$age_z <- as.numeric(scale(rank(dat$age)))
dat$etiv_z <- as.numeric(scale(dat$etiv))

msg("Subjects with brain-change + participant data: %d\n", length(unique(dat$subject_id)))

# ------------------------------------------------------------------------
# Per-ROI x moderator tests
# ------------------------------------------------------------------------
fit_overall <- function(d, mod_change_col) {
  f <- as.formula(paste("brain_change ~", mod_change_col, "+ age_z + sex + etiv_z"))
  fit <- tryCatch(lm(f, data=d), error=function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  s <- summary(fit)$coefficients
  if (!mod_change_col %in% rownames(s)) return(c(NA_real_, NA_real_))
  c(s[mod_change_col, "Estimate"], s[mod_change_col, "Pr(>|t|)"])
}

fit_interaction <- function(d, mod_change_col) {
  f <- as.formula(paste("brain_change ~", mod_change_col, "* dance + age_z + sex + etiv_z"))
  fit <- tryCatch(lm(f, data=d), error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  s <- summary(fit)$coefficients
  int_term <- paste0(mod_change_col, ":dancedance")
  if (!int_term %in% rownames(s)) return(NA_real_)
  s[int_term, "Pr(>|t|)"]
}

results_rows <- list()
data_by_roi_mod <- list()

for (roi_name in roi_set) {
  d_roi <- dat[dat$roi == roi_name, , drop=FALSE]
  if (!nrow(d_roi)) next
  for (mod_col in moderator_cols) {
    change_col <- paste0(mod_col, "_change")
    d_test <- d_roi[!is.na(d_roi[[change_col]]) & !is.na(d_roi$brain_change), , drop=FALSE]
    if (nrow(d_test) < 15) next

    key <- paste(roi_name, mod_col, sep="__")
    data_by_roi_mod[[key]] <- list(data=d_test, change_col=change_col)

    overall <- fit_overall(d_test, change_col)
    interaction_p <- fit_interaction(d_test, change_col)

    results_rows[[key]] <- data.frame(
      roi = roi_name, moderator = mod_col, n_obs = nrow(d_test),
      overall_coupling_estimate = overall[1], overall_coupling_p = overall[2],
      group_interaction_p = interaction_p,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, results_rows)
if (is.null(results_df) || !nrow(results_df)) {
  warning("No ROI x moderator combination had sufficient data")
  quit(status=0)
}

results_df$overall_coupling_p_fdr <- p.adjust(results_df$overall_coupling_p, method="fdr")
results_df$group_interaction_p_fdr <- p.adjust(results_df$group_interaction_p, method="fdr")
write.csv(results_df, file.path(opt$outdir, "coupled_change_summary.csv"), row.names=FALSE)
msg("Wrote coupled_change_summary.csv (%d ROI x moderator combinations)\n", nrow(results_df))

# ------------------------------------------------------------------------
# Targeted LOSO for nominally interesting hits (either test)
# ------------------------------------------------------------------------
screen <- results_df[(!is.na(results_df$overall_coupling_p) & results_df$overall_coupling_p < opt$`loso-screen-threshold`) |
                      (!is.na(results_df$group_interaction_p) & results_df$group_interaction_p < opt$`loso-screen-threshold`), ]
if (nrow(screen)) {
  msg("\nRunning targeted LOSO for %d ROI x moderator hits...\n", nrow(screen))
  loso_rows <- list()
  for (i in seq_len(nrow(screen))) {
    roi_name <- screen$roi[i]; mod_col <- screen$moderator[i]
    key <- paste(roi_name, mod_col, sep="__")
    d_test <- data_by_roi_mod[[key]]$data
    change_col <- data_by_roi_mod[[key]]$change_col
    full_overall <- fit_overall(d_test, change_col)
    full_interaction <- fit_interaction(d_test, change_col)

    for (excl_subj in unique(d_test$subject_id)) {
      d_sub <- d_test[d_test$subject_id != excl_subj, , drop=FALSE]
      ov <- fit_overall(d_sub, change_col)
      ip <- fit_interaction(d_sub, change_col)
      loso_rows[[length(loso_rows) + 1]] <- data.frame(
        roi = roi_name, moderator = mod_col, excluded_subject = excl_subj,
        overall_p = ov[2], full_overall_p = full_overall[2],
        interaction_p = ip, full_interaction_p = full_interaction,
        stringsAsFactors = FALSE
      )
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  write.csv(loso_df, file.path(opt$outdir, "loso_coupled_change.csv"), row.names=FALSE)

  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat("Targeted LOSO for coupled-change hits\n\n", file=con)
  for (i in seq_len(nrow(screen))) {
    roi_name <- screen$roi[i]; mod_col <- screen$moderator[i]
    g <- loso_df[loso_df$roi == roi_name & loso_df$moderator == mod_col, , drop=FALSE]
    if (!nrow(g)) next
    for (test in c("overall","interaction")) {
      p_col <- paste0(test, "_p"); full_col <- paste0("full_", test, "_p")
      gg <- g[!is.na(g[[p_col]]), ]
      if (!nrow(gg)) next
      n_flips <- sum((gg[[p_col]] < opt$alpha) != (gg[[full_col]][1] < opt$alpha))
      most_infl_idx <- which.max(abs(gg[[p_col]] - gg[[full_col]][1]))
      cat(sprintf("%s x %s [%s]: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d flip, most influential = %s (p becomes %.4f)\n",
                  roi_name, mod_col, test, gg[[full_col]][1], min(gg[[p_col]]), max(gg[[p_col]]), n_flips, nrow(gg),
                  gg$excluded_subject[most_infl_idx], gg[[p_col]][most_infl_idx]), file=con)
    }
  }
  msg("LOSO done. See %s/loso_summary.txt\n", opt$outdir)
} else {
  msg("\nNo hits below the LOSO screening threshold.\n")
}

msg("Done.\n")
