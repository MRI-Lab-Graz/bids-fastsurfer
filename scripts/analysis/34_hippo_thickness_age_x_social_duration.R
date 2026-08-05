#!/usr/bin/env Rscript
#
# Does age's role in hippocampal subfield THICKNESS change depend on social
# context (alone vs. small group) or duration (2wk vs 4wk)? Retargets
# 19_age_x_social_duration.R (volume) at thickness -- same two-way design
# (age x social, age x duration as separate terms, not a three-way, per 19's
# power rationale), restricted to the intervention arms only (social/duration
# are undefined for control, as in 13/30).
#
#   change ~ age_z * social + age_z * duration + social * duration +
#            baseline_z + hemisphere + followup_f + sex + (1 | subject_id)
#
# Differences from 19 (mirroring 29/30/33's thickness convention): raw mm
# change (not log_volume change), no etiv_z, --group-column/--social-pattern-*
# /--duration-pattern-* default to this repo's actual group_5 encoding
# (alone_2w/alone_4w/smallgroup_2w/smallgroup_4w), and age_z is RANK-
# transformed as in 33 (not 19's raw scale -- 19 predates that convention).
#
# LOSO runs immediately on any hit, per this project's standing practice.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_thickness_tidy.tsv"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, <group-column>, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/34_hippo_thickness_age_x_social_duration",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character", default="presubiculum,subiculum,CA1,CA2_CA3",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--roi-column"), type="character", default="region",
              help="Column name identifying the subfield in the tidy file [default %default]"),
  make_option(c("--group-column"), type="character", default="group_5",
              help="Column in --participants holding the group assignment [default %default]"),
  make_option(c("--social-pattern-a"), type="character", default="^alone_"),
  make_option(c("--social-label-a"), type="character", default="Alone"),
  make_option(c("--social-pattern-b"), type="character", default="^smallgroup_"),
  make_option(c("--social-label-b"), type="character", default="SmallGroup"),
  make_option(c("--duration-pattern-a"), type="character", default="2w$"),
  make_option(c("--duration-label-a"), type="character", default="2wk"),
  make_option(c("--duration-pattern-b"), type="character", default="4w$"),
  make_option(c("--duration-label-b"), type="character", default="4wk"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10,
              help="Only run LOSO for hits with uncorrected p below this [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE)
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
group_col <- opt$`group-column`

msg("Loading tidy thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$value <- suppressWarnings(as.numeric(tidy$value))
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file", roi_col))

roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

agg <- aggregate(as.formula(paste("value ~ subject_id + session + hemisphere +", roi_col)), data=roi_dat, FUN=mean)
names(agg)[names(agg) == roi_col] <- "subfield"

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
msg("Intervention-arm subjects classified: %d\n", nrow(intervention_participants))
print(table(intervention_participants$social, intervention_participants$duration))

agg <- merge(agg, intervention_participants[, c("subject_id","group","age","sex", "social", "duration")], by="subject_id")
agg$social <- factor(agg$social, levels=c(opt$`social-label-a`, opt$`social-label-b`))
agg$duration <- factor(agg$duration, levels=c(opt$`duration-label-a`, opt$`duration-label-b`))
agg$hemisphere <- factor(agg$hemisphere, levels=sort(unique(agg$hemisphere)))
agg$sex <- factor(agg$sex)
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

# With only one follow-up session (this study's 2-timepoint design),
# followup_f is a constant and would be aliased with the intercept --
# drop it rather than feeding lmer a rank-deficient design (same guard as
# 08/33). Volume's 19_age_x_social_duration.R didn't need this because that
# tidy file has 3 sessions; ours has 2.
build_formula <- function(change_dat) {
  followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
  as.formula(paste0("change ~ age_z * social + age_z * duration + social * duration + ",
                     "baseline_z + hemisphere + ", followup_term, "sex + (1 | subject_id)"))
}

fit_term_p <- function(model, term) {
  at <- tryCatch(anova(model), error=function(e) NULL)
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}

summary_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting age x social/duration model for subfield: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat <- build_change_data(d_roi)
  change_by_roi[[roi_name]] <- change_dat

  model <- tryCatch(lmerTest::lmer(build_formula(change_dat), data=change_dat, REML=TRUE), error=function(e) NULL)
  if (is.null(model)) next

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_model.txt")))
  cat("Subfield:", roi_name, "\n\n")
  print(summary(model))
  sink()

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name,
    n_obs = nrow(change_dat),
    n_subjects = length(unique(change_dat$subject_id)),
    age_x_social_p = fit_term_p(model, "age_z:social"),
    age_x_duration_p = fit_term_p(model, "age_z:duration"),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (is.null(summary_df) || !nrow(summary_df)) stop("No subfield models were successfully fit")

summary_df$age_x_social_p_fdr <- p.adjust(summary_df$age_x_social_p, method="fdr")
summary_df$age_x_duration_p_fdr <- p.adjust(summary_df$age_x_duration_p, method="fdr")
write.csv(summary_df, file.path(opt$outdir, "summary.csv"), row.names=FALSE)
msg("\nWrote summary.csv (FDR-corrected across %d subfields, per term)\n", nrow(summary_df))

# ------------------------------------------------------------------------
# Targeted LOSO for any nominal hit
# ------------------------------------------------------------------------
screen_hits <- list()
for (i in seq_len(nrow(summary_df))) {
  if (!is.na(summary_df$age_x_social_p[i]) && summary_df$age_x_social_p[i] < opt$`loso-screen-threshold`) {
    screen_hits[[length(screen_hits)+1]] <- list(roi=summary_df$roi[i], term="age_z:social")
  }
  if (!is.na(summary_df$age_x_duration_p[i]) && summary_df$age_x_duration_p[i] < opt$`loso-screen-threshold`) {
    screen_hits[[length(screen_hits)+1]] <- list(roi=summary_df$roi[i], term="age_z:duration")
  }
}

if (length(screen_hits)) {
  msg("\nRunning targeted LOSO for %d (subfield, term) hits with uncorrected p < %.2f...\n", length(screen_hits), opt$`loso-screen-threshold`)
  loso_rows <- list()
  for (hit in screen_hits) {
    change_dat <- change_by_roi[[hit$roi]]
    full_model <- tryCatch(lmerTest::lmer(build_formula(change_dat), data=change_dat, REML=TRUE), error=function(e) NULL)
    full_p <- if (!is.null(full_model)) fit_term_p(full_model, hit$term) else NA_real_

    for (excl_subj in unique(change_dat$subject_id)) {
      d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop=FALSE]
      if (length(unique(d_sub$social)) < 2 || length(unique(d_sub$duration)) < 2) next
      fit <- tryCatch(lmerTest::lmer(build_formula(d_sub), data=d_sub, REML=TRUE), error=function(e) NULL)
      if (is.null(fit)) next
      loso_rows[[length(loso_rows)+1]] <- data.frame(
        roi=hit$roi, term=hit$term, excluded_subject=excl_subj,
        loso_p=fit_term_p(fit, hit$term), full_p=full_p, stringsAsFactors=FALSE
      )
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  write.csv(loso_df, file.path(opt$outdir, "loso_hits.csv"), row.names=FALSE)

  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat(sprintf("Targeted LOSO for (subfield, term) hits with uncorrected p < %.2f\n\n", opt$`loso-screen-threshold`), file=con)
  for (hit in screen_hits) {
    g <- loso_df[loso_df$roi == hit$roi & loso_df$term == hit$term & !is.na(loso_df$loso_p), ]
    if (!nrow(g)) next
    n_flips <- sum((g$loso_p < opt$alpha) != (g$full_p[1] < opt$alpha))
    most_infl_idx <- which.max(abs(g$loso_p - g$full_p[1]))
    cat(sprintf("%s x %s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance, most influential = %s (p becomes %.4f)\n",
                hit$roi, hit$term, g$full_p[1], min(g$loso_p), max(g$loso_p), n_flips, nrow(g),
                g$excluded_subject[most_infl_idx], g$loso_p[most_infl_idx]), file=con)
  }
  msg("Done. See %s/summary.csv and %s/loso_summary.txt\n", opt$outdir, opt$outdir)
} else {
  msg("\nNo (subfield, term) combination had uncorrected p < %.2f; nothing to LOSO-check.\n", opt$`loso-screen-threshold`)
}
