#!/usr/bin/env Rscript
#
# EXPLORATORY moderator analysis: does the (baseline-corrected) volume
# change score depend on age, sex, or additional candidate moderators
# (e.g. baseline sleep quality, stress, depression, social support --
# supplied generically via --extra-moderator-cols so this script works
# across studies/instruments without hardcoding column names)?
#
# IMPORTANT: this is run AFTER seeing a null main effect in the primary
# confirmatory analyses, so it is explicitly hypothesis-generating, not
# confirmatory -- testing several moderators x many ROIs after the fact
# carries real fishing risk. Report accordingly: treat any hit here as
# motivation for a dedicated follow-up, not a finding on its own.
#
# LOSO robustness is run automatically, but only for hits that are at
# least nominally interesting (uncorrected p < --loso-screen-threshold) --
# running full leave-one-out on every moderator x ROI combination
# regardless of its p-value would mean tens of thousands of model refits
# for combinations we already know are clearly null. This mirrors how a
# human analyst would triage: screen broadly and cheaply, then validate
# rigorously only the candidates worth a closer look.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(car)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to a tidy volumes TSV (from extract_hippo_subfields.py / extract_amygdala_subfields.py / extract_basal_ganglia.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-m", "--moderators"), type="character", default=NULL,
              help="Optional path to a subject-level moderators TSV (one row per subject_id)"),
  make_option(c("--extra-moderator-cols"), type="character", default=NULL,
              help="Comma-separated column names in --moderators to test as candidate moderators (e.g. 'psqi_global,pss_total' or 'ads_score,stai_total,fsozu_mean')"),
  make_option(c("-o", "--outdir"), type="character", default="results/08_moderator_analysis",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified ROI names [default %default]"),
  make_option(c("--roi-column"), type="character", default="subfield",
              help="Column name identifying the ROI in the tidy file (subfield/nucleus/roi) [default %default]"),
  make_option(c("--intervention-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the intervention contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value used as baseline/covariate [default: earliest session present]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10,
              help="Only run LOSO follow-up on moderator x ROI hits with uncorrected p below this [default %default]"),
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
roi_col <- opt$`roi-column`

# ------------------------------------------------------------------------
# Load and prepare data (change scores, as in 05_baseline_corrected_change.R)
# ------------------------------------------------------------------------
msg("Loading tidy data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file (columns: %s)", roi_col, paste(names(tidy), collapse=", ")))
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check ROI names against the tidy file")

agg <- aggregate(as.formula(paste("volume ~ subject_id + session + hemisphere +", roi_col)), data=roi_dat, FUN=sum)
names(agg)[names(agg) == roi_col] <- "subfield"
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

meta <- participants[, required_pcols]
extra_moderator_cols <- c()
if (!is.null(opt$moderators)) {
  if (!file.exists(opt$moderators)) stop(sprintf("moderators file not found: %s", opt$moderators))
  moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  meta <- merge(meta, moderators, by="subject_id", all.x=TRUE)
  if (!is.null(opt$`extra-moderator-cols`)) {
    extra_moderator_cols <- trimws(strsplit(opt$`extra-moderator-cols`, ",")[[1]])
    missing_mod_cols <- setdiff(extra_moderator_cols, names(moderators))
    if (length(missing_mod_cols)) stop(sprintf("--extra-moderator-cols not found in --moderators: %s", paste(missing_mod_cols, collapse=", ")))
  }
}

agg <- merge(agg, meta, by="subject_id")
agg <- agg[agg$group %in% c(intervention_groups, control_group), , drop=FALSE]
agg$intervention <- factor(ifelse(agg$group %in% intervention_groups, "intervention", "control"), levels=c("control","intervention"))
has_hemisphere <- length(unique(agg$hemisphere)) > 1
agg$hemisphere <- factor(agg$hemisphere, levels=sort(unique(agg$hemisphere)))
# Rank-transformed (not raw) age -- see 01_primary_lmm.R for rationale: this
# cohort's age distribution has a sparse, unevenly-populated tail, and rank
# transformation caps its leverage without discarding subjects. Used both as
# a covariate and (via moderator_specs below) as the "age" moderator itself.
agg$age_z <- as.numeric(scale(rank(agg$age)))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)
agg$sex_mf <- ifelse(agg$sex %in% c("M","F"), agg$sex, NA)
agg$sex_mf <- factor(agg$sex_mf, levels=c("M","F"))

for (col in extra_moderator_cols) {
  agg[[paste0(col, "_z")]] <- as.numeric(scale(suppressWarnings(as.numeric(agg[[col]]))))
}

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
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

# ------------------------------------------------------------------------
# Candidate moderators: base covariates always adjusted for; moderator adds
# an intervention x moderator interaction term on top.
# ------------------------------------------------------------------------
moderator_specs <- list(
  age = list(term = "age_z", label = "Age"),
  sex = list(term = "sex_mf", label = "Sex (M/F only)")
)
for (col in extra_moderator_cols) {
  moderator_specs[[col]] <- list(term = paste0(col, "_z"), label = col)
}

msg("Testing %d candidate moderators across %d ROIs (exploratory)...\n", length(moderator_specs), length(roi_set))

base_covariate_terms <- c("age_z", "sex_mf")

fit_moderator_model <- function(change_dat, mod_term) {
  adjust_terms <- setdiff(base_covariate_terms, mod_term)
  adjust_str <- paste(adjust_terms, collapse=" + ")
  # With only one follow-up session (the 2-timepoint design), followup_f is a
  # constant and would be aliased with the intercept -- drop it rather than
  # feeding lmer a rank-deficient design. Same for hemisphere on midline-only
  # structures (e.g. brainstem), which have a single hemisphere level -- and
  # without a hemisphere split, the change-score data has exactly one row per
  # subject, so a per-subject random intercept is unidentifiable (lme4 needs
  # fewer grouping-factor levels than observations); fall back to plain lm()
  # with Type III SS via car::Anova in that case.
  followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
  hemisphere_term <- if (has_hemisphere) "hemisphere + " else ""
  re_term <- if (has_hemisphere) " + (1 | subject_id)" else ""
  f_str <- paste0("log_change ~ intervention * ", mod_term, " + log_baseline_z + ", hemisphere_term, followup_term,
                   adjust_str, " + etiv_z", re_term)
  interaction_term <- paste0("intervention:", mod_term)
  if (has_hemisphere) {
    fit <- tryCatch(lmerTest::lmer(as.formula(f_str), data=change_dat, REML=TRUE), error=function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    at <- tryCatch(anova(fit), error=function(e) NULL)
  } else {
    fit <- tryCatch(lm(as.formula(f_str), data=change_dat), error=function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    at <- tryCatch(as.data.frame(car::Anova(fit, type=3)), error=function(e) NULL)
  }
  if (!is.null(at) && interaction_term %in% rownames(at)) at[interaction_term, "Pr(>F)"] else NA_real_
}

all_results <- list()
change_by_roi <- list()

for (mod_name in names(moderator_specs)) {
  mod_term <- moderator_specs[[mod_name]]$term
  msg("\n--- Moderator: %s (%s) ---\n", mod_name, moderator_specs[[mod_name]]$label)

  mod_rows <- list()
  for (roi_name in unique(agg$subfield)) {
    d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
    change_dat <- build_change_data(d_roi)
    change_dat <- change_dat[!is.na(change_dat[[mod_term]]), , drop=FALSE]
    change_by_roi[[roi_name]] <- change_dat  # cache for LOSO (moderator-independent columns only used per-call)
    if (!nrow(change_dat) || length(unique(change_dat$intervention)) < 2) next

    p_val <- fit_moderator_model(change_dat, mod_term)
    mod_rows[[roi_name]] <- data.frame(
      moderator = mod_name, roi = roi_name, n_obs = nrow(change_dat),
      interaction_p = p_val, stringsAsFactors = FALSE
    )
  }
  mod_df <- do.call(rbind, mod_rows)
  if (!is.null(mod_df) && nrow(mod_df)) {
    mod_df$interaction_p_fdr <- p.adjust(mod_df$interaction_p, method="fdr")
    mod_df$significant <- mod_df$interaction_p_fdr < opt$alpha
    all_results[[mod_name]] <- mod_df
    for (i in seq_len(nrow(mod_df))) {
      row <- mod_df[i, ]
      msg("  %s: p = %.3f, p_fdr = %.3f%s\n", row$roi, row$interaction_p, row$interaction_p_fdr,
          if (isTRUE(row$significant)) " *" else "")
    }
  }
}

combined <- do.call(rbind, all_results)
if (!is.null(combined)) {
  write.csv(combined, file.path(opt$outdir, "moderator_results.csv"), row.names=FALSE)
}

# ------------------------------------------------------------------------
# Targeted LOSO: only for hits below --loso-screen-threshold
# ------------------------------------------------------------------------
screen_hits <- if (!is.null(combined)) combined[!is.na(combined$interaction_p) & combined$interaction_p < opt$`loso-screen-threshold`, ] else NULL
if (!is.null(screen_hits) && nrow(screen_hits)) {
  msg("\nRunning targeted LOSO for %d moderator x ROI hits with uncorrected p < %.2f...\n", nrow(screen_hits), opt$`loso-screen-threshold`)
  loso_rows <- list()
  for (i in seq_len(nrow(screen_hits))) {
    mod_name <- screen_hits$moderator[i]; roi_name <- screen_hits$roi[i]
    mod_term <- moderator_specs[[mod_name]]$term
    d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
    change_dat <- build_change_data(d_roi)
    change_dat <- change_dat[!is.na(change_dat[[mod_term]]), , drop=FALSE]
    full_p <- fit_moderator_model(change_dat, mod_term)

    for (excl_subj in unique(change_dat$subject_id)) {
      d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop=FALSE]
      if (length(unique(d_sub$intervention)) < 2) next
      p_val <- fit_moderator_model(d_sub, mod_term)
      loso_rows[[length(loso_rows) + 1]] <- data.frame(
        moderator = mod_name, roi = roi_name, excluded_subject = excl_subj,
        loso_p = p_val, full_p = full_p, stringsAsFactors = FALSE
      )
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  write.csv(loso_df, file.path(opt$outdir, "loso_moderator_hits.csv"), row.names=FALSE)

  con_loso <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con_loso), add=TRUE)
  cat(sprintf("Targeted LOSO for moderator x ROI hits with uncorrected p < %.2f\n\n", opt$`loso-screen-threshold`), file=con_loso)
  for (i in seq_len(nrow(screen_hits))) {
    mod_name <- screen_hits$moderator[i]; roi_name <- screen_hits$roi[i]
    g <- loso_df[loso_df$moderator == mod_name & loso_df$roi == roi_name & !is.na(loso_df$loso_p), , drop=FALSE]
    if (!nrow(g)) next
    n_flips <- sum((g$loso_p < opt$alpha) != (g$full_p[1] < opt$alpha))
    most_infl_idx <- which.max(abs(g$loso_p - g$full_p[1]))
    cat(sprintf("%s x %s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance, most influential = %s (p becomes %.4f)\n",
                mod_name, roi_name, g$full_p[1], min(g$loso_p), max(g$loso_p), n_flips, nrow(g),
                g$excluded_subject[most_infl_idx], g$loso_p[most_infl_idx]), file=con_loso)
  }
  msg("LOSO done. See %s/loso_summary.txt\n", opt$outdir)
} else {
  msg("\nNo moderator x ROI combination fell below the LOSO screening threshold (p < %.2f) -- skipping LOSO.\n", opt$`loso-screen-threshold`)
}

con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat("EXPLORATORY moderator analysis (run after a null main effect -- treat any\n",
    "hit here as motivation for a dedicated follow-up study, not a finding)\n\n", file=con)
if (!is.null(combined)) {
  n_sig <- sum(combined$significant, na.rm=TRUE)
  cat(sprintf("Total tests: %d (FDR-corrected within each moderator across %d ROIs)\n", nrow(combined), length(roi_set)), file=con)
  cat(sprintf("Significant after FDR correction: %d\n\n", n_sig), file=con)
  if (n_sig > 0) {
    cat("Significant hits (check loso_summary.txt for robustness before trusting these):\n", file=con)
    sig_rows <- combined[which(combined$significant), ]
    for (i in seq_len(nrow(sig_rows))) {
      row <- sig_rows[i, ]
      cat(sprintf("  - %s moderator on %s: p_fdr = %.4f\n", row$moderator, row$roi, row$interaction_p_fdr), file=con)
    }
  } else {
    cat("No moderator x ROI combination survives FDR correction.\n", file=con)
  }
}
msg("\nDone. See %s/summary.txt\n", opt$outdir)
