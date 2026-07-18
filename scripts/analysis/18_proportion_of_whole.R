#!/usr/bin/env Rscript
#
# Proportion-of-whole-structure analysis: every prior script modelled each
# subfield/nucleus's ABSOLUTE volume (TIV-adjusted via the etiv_z
# covariate). This asks a different, complementary question: does a
# subfield change as a PROPORTION of its parent whole structure (e.g.
# CA3 / Whole_hippocampus)? This is more sensitive to a COMPOSITIONAL
# shift -- one subfield growing disproportionately relative to the rest of
# the same structure (e.g. dentate gyrus/CA-region neurogenesis) -- which
# a model of each subfield's absolute volume, with overall head-size and
# general hippocampal-size variation still present, can dilute.
#
# log(proportion) is used as the outcome (proportions are bounded (0,1);
# log keeps the same "proportional change" interpretation as log(volume)
# elsewhere in this pipeline, and change scores are computed exactly as
# before: log(final_proportion) - log(baseline_proportion)).
#
# Supports both the combined dance+running vs control test (as in
# 15_combined_studies_change.R) and a single-moderator check (e.g. does
# the age x CA3 finding from 08_moderator_analysis.R hold up when CA3 is
# expressed as a proportion of whole hippocampus instead of an absolute,
# TIV-adjusted volume?) via --mode.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to tidy volumes TSV (must include a 'study' column for --mode combined)"),
  make_option(c("-p", "--participants"), type="character", help="Path to participants TSV"),
  make_option(c("--whole-roi"), type="character", required=TRUE, help="ROI value denoting the whole/parent structure (e.g. 'Whole_hippocampus', 'Whole_amygdala')"),
  make_option(c("--roi-column"), type="character", default="subfield"),
  make_option(c("--roi-set"), type="character", required=TRUE, help="Comma-separated subfield/nucleus names to express as a proportion of --whole-roi"),
  make_option(c("--mode"), type="character", default="combined", help="'combined' (dance+running vs control, requires 'study' column) or 'moderator' (single continuous moderator x proportion interaction)"),
  make_option(c("--moderator-col"), type="character", default=NULL, help="For --mode moderator: column in --participants to test as a moderator (e.g. 'age')"),
  make_option(c("-o", "--outdir"), type="character", default="results/18_proportion_of_whole"),
  make_option(c("--dance-groups"), type="character", default=NULL, help="For --mode moderator on a single study: comma-separated intervention group values"),
  make_option(c("--control-group"), type="character", default=NULL),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy) || is.null(opt$participants)) stop("--tidy and --participants are required")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$outdir, "models"), showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-column`
roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])

# ------------------------------------------------------------------------
# Build subfield-as-proportion-of-whole per subject/session/hemisphere
# ------------------------------------------------------------------------
msg("Loading tidy data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

whole <- tidy[tidy[[roi_col]] == opt$`whole-roi`, c("subject_id","session","hemisphere","volume")]
names(whole)[4] <- "whole_volume"
if (!nrow(whole)) stop(sprintf("--whole-roi '%s' not found in tidy file", opt$`whole-roi`))

sub_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
sub_agg <- aggregate(as.formula(paste("volume ~ subject_id + session + hemisphere +", roi_col)), data=sub_dat, FUN=sum)
names(sub_agg)[names(sub_agg) == roi_col] <- "subfield"

merged <- merge(sub_agg, whole, by=c("subject_id","session","hemisphere"))
merged$proportion <- merged$volume / merged$whole_volume
merged$log_proportion <- log(merged$proportion)
msg("Proportion range: [%.4f, %.4f]\n", min(merged$proportion), max(merged$proportion))

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)

sessions <- sort(unique(merged$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
followup_sessions <- setdiff(sessions, baseline_ses)

build_change_data <- function(d_roi) {
  base <- d_roi[d_roi$session == baseline_ses, c("subject_id","hemisphere","log_proportion")]
  names(base)[3] <- "log_baseline_prop"
  fu <- d_roi[d_roi$session %in% followup_sessions, , drop=FALSE]
  m <- merge(fu, base, by=c("subject_id","hemisphere"))
  m$log_change <- m$log_proportion - m$log_baseline_prop
  m$log_baseline_prop_z <- as.numeric(scale(m$log_baseline_prop))
  m$followup_f <- factor(m$session, levels=followup_sessions)
  m
}

fit_term_p <- function(model, term) {
  at <- tryCatch(anova(model), error=function(e) NULL)
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}

if (opt$mode == "combined") {
  if (!"study" %in% names(tidy)) stop("--mode combined requires a 'study' column in the tidy file")
  study_lookup <- unique(tidy[, c("subject_id","study")])
  merged <- merge(merged, study_lookup, by="subject_id")
  merged <- merge(merged, participants[, c("subject_id","group","age","sex")], by="subject_id")
  merged$dance <- factor(ifelse(merged$group == "intervention", "intervention", "control"), levels=c("control","intervention"))
  merged$study <- factor(merged$study)
  merged$hemisphere <- factor(merged$hemisphere, levels=c("lh","rh"))
  merged$sex <- factor(merged$sex)
  merged$age_z <- as.numeric(scale(merged$age))

  base_formula <- log_change ~ dance * study + dance * followup_f + log_baseline_prop_z + hemisphere + age_z + sex

  summary_rows <- list()
  change_by_roi <- list()
  for (roi_name in unique(merged$subfield)) {
    msg("Fitting proportion-of-whole model for %s...\n", roi_name)
    d_roi <- merged[merged$subfield == roi_name, , drop=FALSE]
    change_dat <- build_change_data(d_roi)
    change_by_roi[[roi_name]] <- change_dat
    model <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE), error=function(e) NULL)
    if (is.null(model)) next
    sink(file.path(opt$outdir, "models", paste0(roi_name, "_proportion_model.txt")))
    cat("ROI:", roi_name, "(proportion of", opt$`whole-roi`, ")\n\n"); print(summary(model))
    sink()
    summary_rows[[roi_name]] <- data.frame(
      roi = roi_name, n_obs = nrow(change_dat), n_subjects = length(unique(change_dat$subject_id)),
      dance_main_p = fit_term_p(model, "dance"),
      dance_x_study_p = fit_term_p(model, "dance:study"),
      dance_x_followup_p = fit_term_p(model, "dance:followup_f"),
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  summary_df$dance_main_p_fdr <- p.adjust(summary_df$dance_main_p, method="fdr")
  write.csv(summary_df, file.path(opt$outdir, "proportion_combined_summary.csv"), row.names=FALSE)
  msg("\nWrote proportion_combined_summary.csv\n")

  msg("Running LOSO...\n")
  loso_rows <- list()
  for (roi_name in names(change_by_roi)) {
    change_dat <- change_by_roi[[roi_name]]
    full_model <- lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE)
    full_p <- fit_term_p(full_model, "dance")
    for (excl_subj in unique(change_dat$subject_id)) {
      d_sub <- change_dat[change_dat$subject_id != excl_subj, ]
      if (length(unique(d_sub$dance)) < 2 || length(unique(d_sub$study)) < 2) next
      fit <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d_sub, REML=TRUE), error=function(e) NULL)
      if (is.null(fit)) next
      loso_rows[[length(loso_rows)+1]] <- data.frame(roi=roi_name, excluded_subject=excl_subj, loso_p=fit_term_p(fit,"dance"), full_p=full_p, stringsAsFactors=FALSE)
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  write.csv(loso_df, file.path(opt$outdir, "loso_proportion.csv"), row.names=FALSE)
  con <- file(file.path(opt$outdir, "summary.txt"), open="wt"); on.exit(close(con), add=TRUE)
  cat("Proportion-of-whole-structure combined dance+running vs control\n\n", file=con)
  for (roi_name in names(change_by_roi)) {
    g <- loso_df[loso_df$roi==roi_name & !is.na(loso_df$loso_p),]
    if (!nrow(g)) next
    n_flips <- sum((g$loso_p<opt$alpha)!=(g$full_p[1]<opt$alpha))
    row <- summary_df[summary_df$roi==roi_name,]
    cat(sprintf("%s: p=%.4f, p_fdr=%.4f, LOSO range [%.4f,%.4f], %d/%d flips\n",
                roi_name, row$dance_main_p, row$dance_main_p_fdr, min(g$loso_p), max(g$loso_p), n_flips, nrow(g)), file=con)
  }
  msg("Done. See %s/summary.txt\n", opt$outdir)

} else if (opt$mode == "moderator") {
  if (is.null(opt$`moderator-col`)) stop("--mode moderator requires --moderator-col")
  if (is.null(opt$`dance-groups`) || is.null(opt$`control-group`)) stop("--mode moderator requires --dance-groups and --control-group")
  dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
  control_group <- trimws(opt$`control-group`)

  merged <- merge(merged, participants[, c("subject_id","group","age","sex", opt$`moderator-col`)], by="subject_id")
  merged <- merged[merged$group %in% c(dance_groups, control_group), ]
  merged$dance <- factor(ifelse(merged$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
  merged$hemisphere <- factor(merged$hemisphere, levels=c("lh","rh"))
  merged$sex <- factor(merged$sex)
  merged[[paste0(opt$`moderator-col`,"_z")]] <- as.numeric(scale(merged[[opt$`moderator-col`]]))
  mod_term <- paste0(opt$`moderator-col`, "_z")

  for (roi_name in roi_set) {
    d_roi <- merged[merged$subfield == roi_name, , drop=FALSE]
    change_dat <- build_change_data(d_roi)
    # `dance`, `sex`, and the moderator z-column are already present in
    # change_dat (carried through from d_roi's merge with `merged` above) --
    # no further merge needed here.

    f_str <- paste0("log_change ~ dance * ", mod_term, " + log_baseline_prop_z + hemisphere + followup_f + sex + (1 | subject_id)")
    model <- lmerTest::lmer(as.formula(f_str), data=change_dat, REML=TRUE)
    interaction_term <- paste0("dance:", mod_term)
    full_p <- fit_term_p(model, interaction_term)
    msg("%s x %s (proportion of %s): p = %.4f\n", opt$`moderator-col`, roi_name, opt$`whole-roi`, full_p)

    loso_p <- sapply(unique(change_dat$subject_id), function(s) {
      d_sub <- change_dat[change_dat$subject_id != s, ]
      fit <- tryCatch(lmerTest::lmer(as.formula(f_str), data=d_sub, REML=TRUE), error=function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      fit_term_p(fit, interaction_term)
    })
    loso_p <- loso_p[!is.na(loso_p)]
    n_flips <- sum((loso_p < opt$alpha) != (full_p < opt$alpha))
    con <- file(file.path(opt$outdir, paste0(roi_name, "_", opt$`moderator-col`, "_proportion_summary.txt")), open="wt")
    cat(sprintf("%s x %s, proportion of %s\nfull p = %.4f\nLOSO range [%.4f, %.4f]\n%d/%d exclusions flip significance\n",
                opt$`moderator-col`, roi_name, opt$`whole-roi`, full_p, min(loso_p), max(loso_p), n_flips, length(loso_p)), file=con)
    close(con)
    msg("  LOSO range [%.4f, %.4f], %d/%d flips\n", min(loso_p), max(loso_p), n_flips, length(loso_p))
  }
} else {
  stop("--mode must be 'combined' or 'moderator'")
}
