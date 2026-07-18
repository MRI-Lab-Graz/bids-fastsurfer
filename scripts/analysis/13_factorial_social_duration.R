#!/usr/bin/env Rscript
#
# Factorial analysis within the intervention arms only: does the effect
# differ by SOCIAL CONTEXT (training alone vs. in a small group) or by
# INTERVENTION DURATION (2 weeks vs. 4 weeks)? This is a different question
# than "any intervention vs control" (01/05/06/07) -- it asks whether
# pooling across social context and duration (as the earlier analyses did)
# was hiding a real difference between the study arms themselves.
#
# Design note: Control is excluded here -- social/duration are only defined
# within the 4 intervention arms, so this is a 2x2 factorial (social x
# duration) on the baseline-corrected change score (as in
# 05_baseline_corrected_change.R):
#
#   log_change ~ social * duration + log_baseline_z + hemisphere +
#                followup_f + age_z + sex + etiv_z + (1 | subject_id)
#
# LOSO robustness is run immediately after the main fit (not as a separate
# script) -- given how easily a single subject drove the earlier responder-
# heterogeneity result, any finding here should be checked before being
# reported.

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
  make_option(c("-o", "--outdir"), type="character", default="results/13_factorial_social_duration",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified ROI names [default %default]"),
  make_option(c("--roi-column"), type="character", default="subfield",
              help="Column name identifying the ROI in the tidy file (subfield/nucleus/roi) [default %default]"),
  make_option(c("--social-pattern-a"), type="character", default="Single",
              help="Regex matched against the group column for the first social-context level [default %default]"),
  make_option(c("--social-label-a"), type="character", default="Alone", help="Label for social pattern A [default %default]"),
  make_option(c("--social-pattern-b"), type="character", default="Group",
              help="Regex matched against the group column for the second social-context level [default %default]"),
  make_option(c("--social-label-b"), type="character", default="SmallGroup", help="Label for social pattern B [default %default]"),
  make_option(c("--duration-pattern-a"), type="character", default="2wk",
              help="Regex matched against the group column for the first duration level [default %default]"),
  make_option(c("--duration-label-a"), type="character", default="2wk", help="Label for duration pattern A [default %default]"),
  make_option(c("--duration-pattern-b"), type="character", default="4wk",
              help="Regex matched against the group column for the second duration level [default %default]"),
  make_option(c("--duration-label-b"), type="character", default="4wk", help="Label for duration pattern B [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value used as baseline/covariate [default: earliest session present]"),
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
dir.create(file.path(opt$outdir, "models"), showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_col <- opt$`roi-column`
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file (columns: %s)", roi_col, paste(names(tidy), collapse=", ")))

roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check ROI names against the tidy file")

agg <- aggregate(as.formula(paste("volume ~ subject_id + session + hemisphere +", roi_col)), data=roi_dat, FUN=sum)
names(agg)[names(agg) == roi_col] <- "subfield"
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

is_social_a <- grepl(opt$`social-pattern-a`, participants$group)
is_social_b <- grepl(opt$`social-pattern-b`, participants$group)
is_duration_a <- grepl(opt$`duration-pattern-a`, participants$group)
is_duration_b <- grepl(opt$`duration-pattern-b`, participants$group)

participants$social <- NA_character_
participants$social[is_social_a] <- opt$`social-label-a`
participants$social[is_social_b] <- opt$`social-label-b`
participants$duration <- NA_character_
participants$duration[is_duration_a] <- opt$`duration-label-a`
participants$duration[is_duration_b] <- opt$`duration-label-b`

intervention_participants <- participants[!is.na(participants$social) & !is.na(participants$duration), ]
msg("Intervention-arm subjects classified: %d (social: %s) / factorial cells:\n", nrow(intervention_participants),
    paste(table(intervention_participants$social), collapse=", "))
print(table(intervention_participants$social, intervention_participants$duration))

agg <- merge(agg, intervention_participants[, c(required_pcols, "social", "duration")], by="subject_id")
agg$social <- factor(agg$social, levels=c(opt$`social-label-a`, opt$`social-label-b`))
agg$duration <- factor(agg$duration, levels=c(opt$`duration-label-a`, opt$`duration-label-b`))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

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

base_formula <- log_change ~ social * duration + log_baseline_z + hemisphere + followup_f + age_z + sex + etiv_z

fit_term_p <- function(model, term) {
  at <- tryCatch(anova(model), error=function(e) NULL)
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}

summary_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting factorial model for ROI: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat <- build_change_data(d_roi)
  change_by_roi[[roi_name]] <- change_dat

  model <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE), error=function(e) NULL)
  if (is.null(model)) next

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_factorial_model.txt")))
  cat("ROI:", roi_name, "\n\n")
  print(summary(model))
  sink()

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    n_obs = nrow(change_dat),
    n_subjects = length(unique(change_dat$subject_id)),
    social_p = fit_term_p(model, "social"),
    duration_p = fit_term_p(model, "duration"),
    social_x_duration_p = fit_term_p(model, "social:duration"),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$social_p_fdr <- p.adjust(summary_df$social_p, method="fdr")
  summary_df$duration_p_fdr <- p.adjust(summary_df$duration_p, method="fdr")
  summary_df$social_x_duration_p_fdr <- p.adjust(summary_df$social_x_duration_p, method="fdr")
  write.csv(summary_df, file.path(opt$outdir, "factorial_summary.csv"), row.names=FALSE)
  msg("Wrote factorial_summary.csv (FDR-corrected across %d ROIs, separately per term)\n", nrow(summary_df))
} else {
  warning("No ROI models were successfully fit")
}

# ------------------------------------------------------------------------
# LOSO robustness (run immediately, not as a separate step)
# ------------------------------------------------------------------------
msg("\nRunning LOSO for the factorial model (%d ROIs)...\n", length(change_by_roi))
loso_rows <- list()
for (roi_name in names(change_by_roi)) {
  change_dat <- change_by_roi[[roi_name]]
  all_subjects <- unique(change_dat$subject_id)
  full_model <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE), error=function(e) NULL)
  full_social_p <- if (!is.null(full_model)) fit_term_p(full_model, "social") else NA_real_
  full_duration_p <- if (!is.null(full_model)) fit_term_p(full_model, "duration") else NA_real_
  full_interaction_p <- if (!is.null(full_model)) fit_term_p(full_model, "social:duration") else NA_real_

  for (excl_subj in all_subjects) {
    d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop=FALSE]
    if (length(unique(d_sub$social)) < 2 || length(unique(d_sub$duration)) < 2) next
    fit <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d_sub, REML=TRUE), error=function(e) NULL)
    if (is.null(fit)) next
    loso_rows[[length(loso_rows) + 1]] <- data.frame(
      roi = roi_name, excluded_subject = excl_subj,
      social_p = fit_term_p(fit, "social"), full_social_p = full_social_p,
      duration_p = fit_term_p(fit, "duration"), full_duration_p = full_duration_p,
      interaction_p = fit_term_p(fit, "social:duration"), full_interaction_p = full_interaction_p,
      stringsAsFactors = FALSE
    )
  }
}
loso_df <- do.call(rbind, loso_rows)
if (!is.null(loso_df) && nrow(loso_df)) {
  write.csv(loso_df, file.path(opt$outdir, "loso_factorial.csv"), row.names=FALSE)

  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat("LOSO robustness for the social x duration factorial model\n\n", file=con)
  for (roi_name in names(change_by_roi)) {
    g <- loso_df[loso_df$roi == roi_name, , drop=FALSE]
    if (!nrow(g)) next
    for (term in c("social","duration","interaction")) {
      p_col <- paste0(term, "_p"); full_col <- paste0("full_", term, "_p")
      gg <- g[!is.na(g[[p_col]]), ]
      if (!nrow(gg)) next
      n_flips <- sum((gg[[p_col]] < opt$alpha) != (gg[[full_col]][1] < opt$alpha))
      most_infl_idx <- which.max(abs(gg[[p_col]] - gg[[full_col]][1]))
      cat(sprintf("%s [%s]: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance, most influential = %s (p becomes %.4f)\n",
                  roi_name, term, gg[[full_col]][1], min(gg[[p_col]]), max(gg[[p_col]]), n_flips, nrow(gg),
                  gg$excluded_subject[most_infl_idx], gg[[p_col]][most_infl_idx]), file=con)
    }
  }
  msg("Done. See %s/factorial_summary.csv and %s/loso_summary.txt\n", opt$outdir, opt$outdir)
} else {
  warning("LOSO produced no results")
}
