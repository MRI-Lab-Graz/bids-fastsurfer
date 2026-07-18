#!/usr/bin/env Rscript
#
# Baseline-corrected (ANCOVA-on-change-scores) analysis: does the dance
# intervention change hippocampal subfield volume MORE than control, once
# any pre-existing baseline group difference is explicitly removed?
#
# 01_primary_lmm.R's dance:time_f interaction term already implicitly tests
# this (it's the deviation from the baseline group gap), but its emmeans
# output reports the MARGINAL group difference at each timepoint, which
# still includes the baseline offset and is easy to misread as "the
# intervention worked" when it's actually a static pre-existing gap.
#
# This script makes the correction explicit: change score (log(follow-up) -
# log(baseline)) is the outcome, baseline volume is a covariate (ANCOVA),
# and `dance` now directly tests "extra change beyond what baseline
# differences would predict" (Van Breukelen & van den Brand 2006 -- ANCOVA
# on change scores is more powerful than raw repeated-measures analysis when
# there's baseline imbalance, exactly the situation found in this dataset).

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/05_baseline_corrected_change",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
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
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=roi_dat, FUN=sum)
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(agg, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(dance_groups, control_group), , drop=FALSE]
agg$dance <- factor(ifelse(agg$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
followup_sessions <- setdiff(sessions, baseline_ses)
if (!length(followup_sessions)) stop("No follow-up sessions found beyond the baseline session")
msg("Baseline session: %s. Follow-up sessions: %s\n", baseline_ses, paste(followup_sessions, collapse=", "))

# ------------------------------------------------------------------------
# Build per-ROI change-score datasets: one row per subject x hemisphere x
# follow-up session, with baseline log-volume carried as a covariate.
# ------------------------------------------------------------------------
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

base_formula <- log_change ~ dance * followup_f + log_baseline_z + hemisphere + age_z + sex + etiv_z

fit_one_roi <- function(d, roi_name) {
  fit_re <- function(re_formula) {
    full <- update(base_formula, paste(". ~ . +", re_formula))
    lmerTest::lmer(full, data=d, REML=TRUE)
  }
  model <- tryCatch(fit_re("(1 | subject_id)"), error=function(e) NULL)
  if (is.null(model)) {
    warning(sprintf("Could not fit change-score model for ROI '%s'", roi_name))
    return(NULL)
  }
  model
}

summary_rows <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting baseline-corrected change model for ROI: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat <- build_change_data(d_roi)

  model <- fit_one_roi(change_dat, roi_name)
  if (is.null(model)) next

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_change_model.txt")))
  cat("ROI:", roi_name, "\n")
  cat("Baseline session:", baseline_ses, "-- Follow-up sessions:", paste(followup_sessions, collapse=", "), "\n\n")
  print(summary(model))
  sink()

  anova_tab <- tryCatch(anova(model), error=function(e) NULL)
  dance_row <- if (!is.null(anova_tab) && "dance" %in% rownames(anova_tab)) anova_tab["dance", ] else NULL
  interaction_row <- if (!is.null(anova_tab) && "dance:followup_f" %in% rownames(anova_tab)) anova_tab["dance:followup_f", ] else NULL

  em <- tryCatch(emmeans(model, ~ dance | followup_f), error=function(e) NULL)
  dance_contrasts <- if (!is.null(em)) as.data.frame(contrast(em, method="revpairwise")) else NULL
  if (!is.null(dance_contrasts)) {
    write.csv(dance_contrasts, file.path(opt$outdir, "models", paste0(roi_name, "_change_dance_vs_control_by_followup.csv")), row.names=FALSE)
  }

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    n_obs = nrow(change_dat),
    n_subjects = length(unique(change_dat$subject_id)),
    dance_main_F = if (!is.null(dance_row)) dance_row[["F value"]] else NA_real_,
    dance_main_p = if (!is.null(dance_row)) dance_row[["Pr(>F)"]] else NA_real_,
    dance_x_followup_F = if (!is.null(interaction_row)) interaction_row[["F value"]] else NA_real_,
    dance_x_followup_p = if (!is.null(interaction_row)) interaction_row[["Pr(>F)"]] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$dance_main_p_fdr <- p.adjust(summary_df$dance_main_p, method="fdr")
  summary_df$dance_x_followup_p_fdr <- p.adjust(summary_df$dance_x_followup_p, method="fdr")
  summary_df$significant_dance_main <- summary_df$dance_main_p_fdr < opt$alpha
  summary_df$significant_dance_x_followup <- summary_df$dance_x_followup_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$dance_main_p_fdr), ]
  write.csv(summary_df, file.path(opt$outdir, "baseline_corrected_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (FDR-corrected across %d ROIs) to baseline_corrected_summary.csv\n", nrow(summary_df))
  msg("\nInterpretation:\n- 'dance_main' tests whether dance shows MORE change than control, baseline-corrected (this is the clean intervention-effect test).\n- 'dance_x_followup' tests whether that extra change differs between the two follow-up sessions (i.e. still accumulating vs already plateaued).\n")
} else {
  warning("No ROI change-score models were successfully fit; no summary written")
}
