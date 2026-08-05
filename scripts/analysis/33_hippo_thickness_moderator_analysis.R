#!/usr/bin/env Rscript
#
# EXPLORATORY moderator analysis for hippocampal subfield THICKNESS (hipsta),
# retargeting 08_moderator_analysis.R (volume) at thickness: does the
# baseline-corrected thickness change score depend on age, sex, or additional
# candidate psychological moderators (CES-D/ads_score, PANAS, STAI, FSozU --
# supplied generically via --extra-moderator-cols, as in 08)?
#
# Differences from 08 (mirroring 29/30's thickness convention, not 08's
# volume one):
#   - raw mm CHANGE score, not log(volume) change -- thickness isn't
#     eTIV-scaled and shares one scale across subfields, so the "scale
#     domination" rationale for log(volume) doesn't apply (see
#     29_hippo_thickness_lmm.R and extract_hipsta_thickness.py).
#   - no etiv_z covariate, for the same reason.
#   - age_z is RANK-transformed (as in 29/30/01), not raw scale (unlike 08,
#     which predates that convention) -- this cohort's age distribution is
#     concentrated at 18-27 with a sparse tail to 38.
#   - --group-column/--intervention-groups/--control-group default to this
#     repo's actual study 129 group_5 encoding (alone_2w/alone_4w/
#     smallgroup_2w/smallgroup_4w/control), not 08's placeholder
#     "ballet,contemporary" (that default was written for a different,
#     unrelated dance-intervention study template and doesn't match any
#     column in this project's configs/participants.tsv).
#   - --roi-column defaults to "region" and --roi-set to hipsta's four
#     available subfields (presubiculum/subiculum/CA1/CA2_CA3), matching
#     extract_hipsta_thickness.py's tidy contract.
#
# IMPORTANT: exploratory, not confirmatory -- same fishing-risk caveat as 08.
# LOSO runs automatically, but only for hits below --loso-screen-threshold.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(car)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_thickness_tidy.tsv (from extract_hipsta_thickness.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, <group-column>, age, sex"),
  make_option(c("-m", "--moderators"), type="character", default=NULL,
              help="Optional path to a subject-level moderators TSV (one row per subject_id)"),
  make_option(c("--extra-moderator-cols"), type="character", default=NULL,
              help="Comma-separated column names in --moderators to test as candidate moderators (e.g. 'ads_score,stai_total,fsozu_mean,panas_pa,panas_na')"),
  make_option(c("-o", "--outdir"), type="character", default="results/33_hippo_thickness_moderator_analysis",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character", default="presubiculum,subiculum,CA1,CA2_CA3",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--roi-column"), type="character", default="region",
              help="Column name identifying the subfield in the tidy file [default %default]"),
  make_option(c("--group-column"), type="character", default="group_5",
              help="Column in --participants holding the group assignment [default %default]"),
  make_option(c("--intervention-groups"), type="character",
              default="alone_2w,alone_4w,smallgroup_2w,smallgroup_4w",
              help="Comma-separated group values pooled into the intervention contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value used as baseline/covariate [default: earliest session present]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10,
              help="Only run LOSO follow-up on moderator x subfield hits with uncorrected p below this [default %default]"),
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
group_col <- opt$`group-column`
roi_col <- opt$`roi-column`

# ------------------------------------------------------------------------
# Load and prepare data (raw-mm change scores, not log)
# ------------------------------------------------------------------------
msg("Loading tidy thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$value <- suppressWarnings(as.numeric(tidy$value))

if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file (columns: %s)", roi_col, paste(names(tidy), collapse=", ")))
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(as.formula(paste("value ~ subject_id + session + hemisphere +", roi_col)), data=roi_dat, FUN=mean)
names(agg)[names(agg) == roi_col] <- "subfield"

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"

meta <- participants[, c("subject_id","group","age","sex")]
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
agg$age_z <- as.numeric(scale(rank(agg$age)))
agg$sex_mf <- ifelse(agg$sex %in% c("M","F"), agg$sex, NA)
agg$sex_mf <- factor(agg$sex_mf, levels=c("M","F"))

for (col in extra_moderator_cols) {
  agg[[paste0(col, "_z")]] <- as.numeric(scale(suppressWarnings(as.numeric(agg[[col]]))))
}

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
followup_sessions <- setdiff(sessions, baseline_ses)

build_change_data <- function(d_roi) {
  base <- d_roi[d_roi$session == baseline_ses, c("subject_id","hemisphere","value")]
  names(base)[3] <- "baseline"
  fu <- d_roi[d_roi$session %in% followup_sessions, , drop=FALSE]
  merged <- merge(fu, base, by=c("subject_id","hemisphere"))
  merged$change <- merged$value - merged$baseline
  merged$baseline_z <- as.numeric(scale(merged$baseline))
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

msg("Testing %d candidate moderators across %d subfields (exploratory)...\n", length(moderator_specs), length(roi_set))

base_covariate_terms <- c("age_z", "sex_mf")

fit_moderator_model <- function(change_dat, mod_term) {
  adjust_terms <- setdiff(base_covariate_terms, mod_term)
  adjust_str <- paste(adjust_terms, collapse=" + ")
  # With only one follow-up session, followup_f is constant and aliased with
  # the intercept -- drop it (mirrors 08). Random intercept requires >1
  # hemisphere row per subject; fall back to lm()+car::Anova otherwise.
  followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
  hemisphere_term <- if (has_hemisphere) "hemisphere + " else ""
  re_term <- if (has_hemisphere) " + (1 | subject_id)" else ""
  f_str <- paste0("change ~ intervention * ", mod_term, " + baseline_z + ", hemisphere_term, followup_term,
                   adjust_str, re_term)
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
    change_by_roi[[roi_name]] <- change_dat
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
  msg("\nRunning targeted LOSO for %d moderator x subfield hits with uncorrected p < %.2f...\n", nrow(screen_hits), opt$`loso-screen-threshold`)
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
  cat(sprintf("Targeted LOSO for moderator x subfield hits with uncorrected p < %.2f\n\n", opt$`loso-screen-threshold`), file=con_loso)
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
  msg("\nNo moderator x subfield combination fell below the LOSO screening threshold (p < %.2f) -- skipping LOSO.\n", opt$`loso-screen-threshold`)
}

con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat("EXPLORATORY moderator analysis for hippocampal subfield THICKNESS (run\n",
    "after a null main effect -- treat any hit here as motivation for a\n",
    "dedicated follow-up, not a finding)\n\n", file=con)
if (!is.null(combined)) {
  n_sig <- sum(combined$significant, na.rm=TRUE)
  cat(sprintf("Total tests: %d (FDR-corrected within each moderator across %d subfields)\n", nrow(combined), length(roi_set)), file=con)
  cat(sprintf("Significant after FDR correction: %d\n\n", n_sig), file=con)
  if (n_sig > 0) {
    cat("Significant hits (check loso_summary.txt for robustness before trusting these):\n", file=con)
    sig_rows <- combined[which(combined$significant), ]
    for (i in seq_len(nrow(sig_rows))) {
      row <- sig_rows[i, ]
      cat(sprintf("  - %s moderator on %s: p_fdr = %.4f\n", row$moderator, row$roi, row$interaction_p_fdr), file=con)
    }
  } else {
    cat("No moderator x subfield combination survives FDR correction.\n", file=con)
  }
}
msg("\nDone. See %s/summary.txt\n", opt$outdir)
