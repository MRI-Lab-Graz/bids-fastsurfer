#!/usr/bin/env Rscript
#
# Does age's relationship to hippocampal subfield THICKNESS change depend on
# sex, or on a psychological trait/state (CES-D/ads_score depression, PANAS
# affect, STAI anxiety, FSozU social support)? Retargets
# 20_age_x_moderator_x_dance.R (volume) at thickness: tests a three-way
# intervention x age_z x third_var interaction for each candidate third
# variable, across the pre-specified subfield set.
#
# Sex is always tested; each survey moderator is tested as a SEPARATE
# three-way model (not jointly), same power rationale as 20.
#
# Differences from 20 (mirroring 29/30/33/34's thickness convention): raw mm
# change (not log_volume change), no etiv_z, --group-column/--intervention-
# groups/--control-group default to this repo's actual group_5 encoding, and
# age_z is RANK-transformed as in 33/34 (not 20's raw scale).
#
# LOSO runs only on hits below --loso-screen-threshold, mirroring 20/33's
# triage logic -- exploratory (post-hoc), treat any hit as motivation for a
# dedicated follow-up, not a finding on its own, until LOSO-checked.

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
  make_option(c("-m", "--moderators"), type="character", default=NULL,
              help="Optional path to a subject-level moderators TSV"),
  make_option(c("--extra-moderator-cols"), type="character", default=NULL,
              help="Comma-separated column names in --moderators to test as third variables (e.g. 'ads_score,panas_pa,panas_na,fsozu_mean,stai_total')"),
  make_option(c("-o", "--outdir"), type="character", default="results/35_hippo_thickness_age_x_moderator_x_intervention",
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
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10),
  make_option(c("--quiet"), action="store_true", default=FALSE)
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
# Third variables: sex always tested; extras from --extra-moderator-cols
# ------------------------------------------------------------------------
third_var_specs <- list(sex = list(term = "sex_mf", label = "Sex (M/F only)"))
for (col in extra_moderator_cols) {
  third_var_specs[[col]] <- list(term = paste0(col, "_z"), label = col)
}

fit_three_way <- function(change_dat, third_term) {
  # As in 33/34: with only one follow-up session, followup_f is constant and
  # aliased with the intercept -- drop it rather than feeding lmer a
  # rank-deficient design.
  followup_term <- if (nlevels(change_dat$followup_f) > 1) "followup_f + " else ""
  f_str <- paste0("change ~ intervention * age_z * ", third_term,
                   " + baseline_z + hemisphere + ", followup_term, "(1 | subject_id)")
  fit <- tryCatch(lmerTest::lmer(as.formula(f_str), data=change_dat, REML=TRUE), error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  at <- tryCatch(anova(fit), error=function(e) NULL)
  interaction_term <- paste0("intervention:age_z:", third_term)
  if (!is.null(at) && interaction_term %in% rownames(at)) at[interaction_term, "Pr(>F)"] else NA_real_
}

msg("Testing %d third variables x %d subfields for intervention:age_z:third_var (exploratory)...\n", length(third_var_specs), length(roi_set))

all_results <- list()
change_by_roi <- list()

for (tv_name in names(third_var_specs)) {
  tv_term <- third_var_specs[[tv_name]]$term
  msg("\n--- Third variable: %s (%s) ---\n", tv_name, third_var_specs[[tv_name]]$label)

  rows <- list()
  for (roi_name in unique(agg$subfield)) {
    d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
    change_dat <- build_change_data(d_roi)
    change_dat <- change_dat[!is.na(change_dat[[tv_term]]), , drop=FALSE]
    change_by_roi[[paste(tv_name, roi_name, sep="__")]] <- change_dat
    if (!nrow(change_dat) || length(unique(change_dat$intervention)) < 2) next

    p_val <- fit_three_way(change_dat, tv_term)
    rows[[roi_name]] <- data.frame(
      third_var = tv_name, roi = roi_name, n_obs = nrow(change_dat),
      interaction_p = p_val, stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  if (!is.null(df) && nrow(df)) {
    df$interaction_p_fdr <- p.adjust(df$interaction_p, method="fdr")
    df$significant <- df$interaction_p_fdr < opt$alpha
    all_results[[tv_name]] <- df
    for (i in seq_len(nrow(df))) {
      row <- df[i, ]
      msg("  %s: p = %.3f, p_fdr = %.3f%s\n", row$roi, row$interaction_p, row$interaction_p_fdr,
          if (row$significant) " *" else "")
    }
  }
}

combined <- do.call(rbind, all_results)
if (!is.null(combined)) write.csv(combined, file.path(opt$outdir, "three_way_results.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# Targeted LOSO for hits below --loso-screen-threshold
# ------------------------------------------------------------------------
screen_hits <- if (!is.null(combined)) combined[!is.na(combined$interaction_p) & combined$interaction_p < opt$`loso-screen-threshold`, ] else NULL
if (!is.null(screen_hits) && nrow(screen_hits)) {
  msg("\nRunning targeted LOSO for %d (third_var, subfield) hits with uncorrected p < %.2f...\n", nrow(screen_hits), opt$`loso-screen-threshold`)
  loso_rows <- list()
  for (i in seq_len(nrow(screen_hits))) {
    tv_name <- screen_hits$third_var[i]; roi_name <- screen_hits$roi[i]
    tv_term <- third_var_specs[[tv_name]]$term
    change_dat <- change_by_roi[[paste(tv_name, roi_name, sep="__")]]
    full_p <- fit_three_way(change_dat, tv_term)

    for (excl_subj in unique(change_dat$subject_id)) {
      d_sub <- change_dat[change_dat$subject_id != excl_subj, , drop=FALSE]
      if (length(unique(d_sub$intervention)) < 2) next
      p_val <- fit_three_way(d_sub, tv_term)
      loso_rows[[length(loso_rows) + 1]] <- data.frame(
        third_var = tv_name, roi = roi_name, excluded_subject = excl_subj,
        loso_p = p_val, full_p = full_p, stringsAsFactors = FALSE
      )
    }
  }
  loso_df <- do.call(rbind, loso_rows)
  write.csv(loso_df, file.path(opt$outdir, "loso_hits.csv"), row.names=FALSE)

  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat(sprintf("Targeted LOSO for (third_var, subfield) hits with uncorrected p < %.2f\n\n", opt$`loso-screen-threshold`), file=con)
  for (i in seq_len(nrow(screen_hits))) {
    tv_name <- screen_hits$third_var[i]; roi_name <- screen_hits$roi[i]
    g <- loso_df[loso_df$third_var == tv_name & loso_df$roi == roi_name & !is.na(loso_df$loso_p), ]
    if (!nrow(g)) next
    n_flips <- sum((g$loso_p < opt$alpha) != (g$full_p[1] < opt$alpha))
    most_infl_idx <- which.max(abs(g$loso_p - g$full_p[1]))
    cat(sprintf("%s x %s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance, most influential = %s (p becomes %.4f)\n",
                tv_name, roi_name, g$full_p[1], min(g$loso_p), max(g$loso_p), n_flips, nrow(g),
                g$excluded_subject[most_infl_idx], g$loso_p[most_infl_idx]), file=con)
  }
  msg("Done. See %s/three_way_results.csv and %s/loso_summary.txt\n", opt$outdir, opt$outdir)
} else {
  msg("\nNo (third_var, subfield) combination had uncorrected p < %.2f; nothing to LOSO-check.\n", opt$`loso-screen-threshold`)
}
