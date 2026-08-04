#!/usr/bin/env Rscript
#
# Factorial analysis within the intervention arms only: does the hippocampal
# subfield THICKNESS effect differ by SOCIAL CONTEXT (training alone vs. in a
# small group) or by INTERVENTION DURATION (2 weeks vs. 4 weeks)? This is a
# different question than "any intervention vs control" (29_hippo_thickness_lmm.R)
# -- it asks whether pooling across social context and duration (as 29 does)
# is hiding a real difference between the study arms themselves. Mirrors
# 13_factorial_social_duration.R (volume) exactly, retargeted at thickness.
#
# Design note: Control is excluded here -- social/duration are only defined
# within the intervention arms, so this is a 2x2 factorial (social x
# duration) on the baseline-corrected change score:
#
#   change ~ social * duration + baseline_z + hemisphere +
#            followup_f + age_z + sex + (1 | subject_id)
#
# Unlike 13 (log_change, for volume), this uses the RAW change score -- see
# 27_cortical_thickness_lmm.R / 29_hippo_thickness_lmm.R for why thickness is
# not log-transformed or eTIV-corrected the way volume is.
#
# --group-column/--social-pattern-*/--duration-pattern-* default to this
# repo's actual configs/participants.tsv encoding for study 129 (group_5:
# control/alone_2w/alone_4w/smallgroup_2w/smallgroup_4w). Between this script
# (social x duration, decomposed) and 29 (any-intervention vs control,
# pooled), every grouping "angle" on the same 5-level design is covered by
# one or the other; there's no need for a third script testing e.g. "2wk vs
# 4wk" alone or "alone vs group" alone in isolation -- those are exactly the
# social/duration MAIN EFFECTS already reported here, tested jointly with
# their interaction rather than as separate uncorrected tests.
#
# LOSO robustness is run immediately after the main fit (not as a separate
# script), matching 13's convention -- any finding here should be checked
# before being reported.

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
  make_option(c("-o", "--outdir"), type="character", default="results/30_hippo_thickness_factorial",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character", default="presubiculum,subiculum,CA1,CA2_CA3",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--group-column"), type="character", default="group_5",
              help="Column in --participants holding the group assignment [default %default]"),
  make_option(c("--social-pattern-a"), type="character", default="^alone",
              help="Regex matched against the group column for the first social-context level [default %default]"),
  make_option(c("--social-label-a"), type="character", default="Alone", help="Label for social pattern A [default %default]"),
  make_option(c("--social-pattern-b"), type="character", default="^smallgroup",
              help="Regex matched against the group column for the second social-context level [default %default]"),
  make_option(c("--social-label-b"), type="character", default="SmallGroup", help="Label for social pattern B [default %default]"),
  make_option(c("--duration-pattern-a"), type="character", default="2w$",
              help="Regex matched against the group column for the first duration level [default %default]"),
  make_option(c("--duration-label-a"), type="character", default="2wk", help="Label for duration pattern A [default %default]"),
  make_option(c("--duration-pattern-b"), type="character", default="4w$",
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
group_col <- opt$`group-column`

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
msg("Loading tidy hippocampal thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$value <- suppressWarnings(as.numeric(tidy$value))

roi_dat <- tidy[tidy$region %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(value ~ subject_id + session + hemisphere + region, data=roi_dat, FUN=mean)
names(agg)[names(agg) == "region"] <- "subfield"

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"

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
msg("Intervention-arm subjects classified: %d / factorial cells:\n", nrow(intervention_participants))
print(table(intervention_participants$social, intervention_participants$duration))

agg <- merge(agg, intervention_participants[, c("subject_id","group","age","sex","social","duration")], by="subject_id")
agg$social <- factor(agg$social, levels=c(opt$`social-label-a`, opt$`social-label-b`))
agg$duration <- factor(agg$duration, levels=c(opt$`duration-label-a`, opt$`duration-label-b`))
has_hemisphere <- length(unique(agg$hemisphere)) > 1
agg$hemisphere <- factor(agg$hemisphere, levels=sort(unique(agg$hemisphere)))
agg$sex <- factor(agg$sex)
# Rank-transformed age -- see 01_primary_lmm.R for rationale.
agg$age_z <- as.numeric(scale(rank(agg$age)))

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

# With only one follow-up session, followup_f is a constant and would be
# aliased with the intercept -- drop it rather than feeding lmer a
# rank-deficient design (same handling as 13_factorial_social_duration.R).
followup_term <- if (length(followup_sessions) > 1) "followup_f + " else ""
hemisphere_term <- if (has_hemisphere) "hemisphere + " else ""
base_formula <- as.formula(paste("change ~ social * duration + baseline_z +", hemisphere_term, followup_term, "age_z + sex"))

fit_term_p <- function(model, term) {
  if (inherits(model, "lm")) {
    at <- tryCatch(as.data.frame(car::Anova(model, type=3)), error=function(e) NULL)
  } else {
    at <- tryCatch(anova(model), error=function(e) NULL)
  }
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}

# Without a hemisphere split, a per-subject random intercept would be
# unidentifiable (one row per subject) -- fall back to plain lm() with Type
# III SS in that case, same as 13.
fit_factorial_model <- function(d) {
  if (has_hemisphere) {
    tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d, REML=TRUE), error=function(e) NULL)
  } else {
    tryCatch(lm(base_formula, data=d), error=function(e) NULL)
  }
}

summary_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting factorial model for subfield: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat <- build_change_data(d_roi)
  change_by_roi[[roi_name]] <- change_dat

  model <- fit_factorial_model(change_dat)
  if (is.null(model)) next

  sink(file.path(opt$outdir, "models", paste0(gsub("[^0-9A-Za-z_]", "_", roi_name), "_factorial_model.txt")))
  cat("Subfield:", roi_name, "\n\n")
  print(summary(model))
  sink()

  summary_rows[[roi_name]] <- data.frame(
    subfield = roi_name,
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
  msg("Wrote factorial_summary.csv (FDR-corrected across %d subfields, separately per term)\n", nrow(summary_df))
} else {
  warning("No subfield models were successfully fit")
}

# ------------------------------------------------------------------------
# LOSO robustness (run immediately, not as a separate step)
# ------------------------------------------------------------------------
msg("\nRunning LOSO for the factorial model (%d subfields)...\n", length(change_by_roi))
loso_rows <- list()
for (roi_name in names(change_by_roi)) {
  change_dat <- change_by_roi[[roi_name]]
  all_subjects <- unique(change_dat$subject_id)
  full_model <- fit_factorial_model(change_dat)
  full_social_p <- if (!is.null(full_model)) fit_term_p(full_model, "social") else NA_real_
  full_duration_p <- if (!is.null(full_model)) fit_term_p(full_model, "duration") else NA_real_
  full_interaction_p <- if (!is.null(full_model)) fit_term_p(full_model, "social:duration") else NA_real_

  for (excl_subj in all_subjects) {
    d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop=FALSE]
    if (length(unique(d_sub$social)) < 2 || length(unique(d_sub$duration)) < 2) next
    fit <- fit_factorial_model(d_sub)
    if (is.null(fit)) next
    loso_rows[[length(loso_rows) + 1]] <- data.frame(
      subfield = roi_name, excluded_subject = excl_subj,
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
  cat("LOSO robustness for the social x duration factorial model (hippocampal thickness)\n\n", file=con)
  for (roi_name in names(change_by_roi)) {
    g <- loso_df[loso_df$subfield == roi_name, , drop=FALSE]
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
