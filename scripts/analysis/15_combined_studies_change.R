#!/usr/bin/env Rscript
#
# Combined-studies analysis: pools the dance intervention (study 134:
# ballet + contemporary) and the running intervention (study 129: all four
# arms) into a single "intervention" contrast against both studies' control
# groups, on the baseline-corrected change score (as in
# 05_baseline_corrected_change.R).
#
# Pooling across studies with different protocols/durations needs care, so
# this script tests THREE things per ROI, not one:
#   1. dance_main: does "any structured exercise" show more change than
#      control, pooling both studies (the user's actual question)?
#   2. dance_x_study: does the intervention effect differ between the two
#      studies -- i.e. is a pooled "hit" actually driven by just one study?
#      A significant interaction here means the pooled main effect should
#      NOT be interpreted as a single unified effect.
#   3. dance_x_followup: does the effect differ between follow-up sessions
#      (as in 05), now pooled across studies.
#
# `study` is included as an addititve covariate throughout (accounts for
# any systematic between-study offset in scanner/protocol/population),
# and LOSO robustness is run automatically on every term, given how many
# leads have already turned out to be single-subject artifacts in this
# project.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to combined tidy volumes TSV (with a 'study' column)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to combined participants TSV with columns: subject_id, study, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/15_combined_studies_change",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified ROI names [default %default]"),
  make_option(c("--roi-column"), type="character", default="subfield",
              help="Column name identifying the ROI in the tidy file (subfield/nucleus/roi) [default %default]"),
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
roi_col <- opt$`roi-column`

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
msg("Loading combined tidy data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file", roi_col))
if (!"study" %in% names(tidy)) stop("tidy file must have a 'study' column")

roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check ROI names against the tidy file")

agg <- aggregate(as.formula(paste("volume ~ subject_id + session + hemisphere + study +", roi_col)), data=roi_dat, FUN=sum)
names(agg)[names(agg) == roi_col] <- "subfield"
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","study","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(agg, participants[, c("subject_id","group","age","sex")], by="subject_id")
agg$dance <- factor(ifelse(agg$group == "intervention", "intervention", "control"), levels=c("control","intervention"))
# Levels are taken from the data rather than hardcoded to lh/rh so that
# midline structures (hemisphere == "midline", e.g. Brain-Stem, corpus
# callosum segments) work too: with a single observed level, hemisphere
# drops out of the model matrix instead of becoming all-NA and dropping
# every row.
agg$hemisphere <- factor(agg$hemisphere, levels=unique(agg$hemisphere))
agg$study <- factor(agg$study)
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

msg("Subjects: %d (%s)\n", length(unique(agg$subject_id)),
    paste(names(table(participants$study[participants$subject_id %in% agg$subject_id])),
          table(participants$study[participants$subject_id %in% agg$subject_id]), sep="=", collapse=", "))

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
followup_sessions <- setdiff(sessions, baseline_ses)
msg("Baseline: %s. Follow-ups: %s\n", baseline_ses, paste(followup_sessions, collapse=", "))

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

base_formula <- log_change ~ dance * study + dance * followup_f + log_baseline_z + hemisphere + age_z + sex + etiv_z

fit_term_p <- function(model, term) {
  at <- tryCatch(anova(model), error=function(e) NULL)
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}

# Midline structures (Brain-Stem, corpus callosum segments) have a single
# hemisphere level; R's contrast machinery errors on a 1-level factor term
# ("contrasts can be applied only to factors with 2 or more levels"), so
# `hemisphere` is dropped from the formula for those ROIs specifically.
roi_formula <- function(change_dat) {
  if (nlevels(droplevels(change_dat$hemisphere)) < 2) update(base_formula, . ~ . - hemisphere) else base_formula
}

summary_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting combined-studies model for ROI: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat <- build_change_data(d_roi)
  change_by_roi[[roi_name]] <- change_dat

  model <- tryCatch(lmerTest::lmer(update(roi_formula(change_dat), ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE), error=function(e) NULL)
  if (is.null(model)) next

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_combined_model.txt")))
  cat("ROI:", roi_name, "\n\n")
  print(summary(model))
  sink()

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    n_obs = nrow(change_dat),
    n_subjects = length(unique(change_dat$subject_id)),
    dance_main_p = fit_term_p(model, "dance"),
    dance_x_study_p = fit_term_p(model, "dance:study"),
    dance_x_followup_p = fit_term_p(model, "dance:followup_f"),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (is.null(summary_df) || !nrow(summary_df)) stop("No ROI models were successfully fit")

summary_df$dance_main_p_fdr <- p.adjust(summary_df$dance_main_p, method="fdr")
summary_df$dance_x_study_p_fdr <- p.adjust(summary_df$dance_x_study_p, method="fdr")
summary_df$dance_x_followup_p_fdr <- p.adjust(summary_df$dance_x_followup_p, method="fdr")
write.csv(summary_df, file.path(opt$outdir, "combined_summary.csv"), row.names=FALSE)
msg("\nWrote combined_summary.csv (FDR-corrected across %d ROIs, per term)\n", nrow(summary_df))
msg("Interpretation: 'dance_main' = pooled intervention effect (the main question). 'dance_x_study'\n")
msg("significant means the pooled effect is NOT consistent across studies -- don't trust dance_main\n")
msg("in that case without checking the two studies separately.\n")

# ------------------------------------------------------------------------
# LOSO robustness for every ROI x term (full, not just screened -- combined
# n is large enough and this is the user's primary question of interest)
# ------------------------------------------------------------------------
msg("\nRunning LOSO for the combined model (%d ROIs)...\n", length(change_by_roi))
loso_rows <- list()
for (roi_name in names(change_by_roi)) {
  change_dat <- change_by_roi[[roi_name]]
  all_subjects <- unique(change_dat$subject_id)
  roi_f <- roi_formula(change_dat)
  full_model <- tryCatch(lmerTest::lmer(update(roi_f, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE), error=function(e) NULL)
  full_dance_p <- if (!is.null(full_model)) fit_term_p(full_model, "dance") else NA_real_
  full_study_p <- if (!is.null(full_model)) fit_term_p(full_model, "dance:study") else NA_real_
  full_followup_p <- if (!is.null(full_model)) fit_term_p(full_model, "dance:followup_f") else NA_real_

  for (excl_subj in all_subjects) {
    d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop=FALSE]
    if (length(unique(d_sub$dance)) < 2 || length(unique(d_sub$study)) < 2) next
    fit <- tryCatch(lmerTest::lmer(update(roi_f, ". ~ . + (1 | subject_id)"), data=d_sub, REML=TRUE), error=function(e) NULL)
    if (is.null(fit)) next
    loso_rows[[length(loso_rows) + 1]] <- data.frame(
      roi = roi_name, excluded_subject = excl_subj,
      dance_p = fit_term_p(fit, "dance"), full_dance_p = full_dance_p,
      study_p = fit_term_p(fit, "dance:study"), full_study_p = full_study_p,
      followup_p = fit_term_p(fit, "dance:followup_f"), full_followup_p = full_followup_p,
      stringsAsFactors = FALSE
    )
  }
}
loso_df <- do.call(rbind, loso_rows)
write.csv(loso_df, file.path(opt$outdir, "loso_combined.csv"), row.names=FALSE)

con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat("LOSO robustness for the combined dance+running vs control model\n\n", file=con)
for (roi_name in names(change_by_roi)) {
  g <- loso_df[loso_df$roi == roi_name, , drop=FALSE]
  if (!nrow(g)) next
  for (term in c("dance","study","followup")) {
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
msg("Done. See %s/combined_summary.csv and %s/loso_summary.txt\n", opt$outdir, opt$outdir)
