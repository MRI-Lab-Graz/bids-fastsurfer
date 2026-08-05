#!/usr/bin/env Rscript
#
# HEMISPHERE x TIME x GROUP x CES-D (baseline depression, continuous) for
# hippocampal subfield THICKNESS. Retargets 23_hem_time_group_cesd.R (volume)
# at thickness -- same rationale: an LMM handles the continuous CES-D x
# categorical GROUP x categorical TIME x categorical HEMISPHERE four-way
# interaction without classical rm-ANOVA's sum-of-squares ambiguity for
# covariate-by-within-factor terms, and uses all available data rather than
# requiring complete cases across every session.
#
#   value ~ hemisphere * time_f * group3 * cesd_z + age_z + sex
#           + (1 | subject_id)
#
# GROUP3 is the 3-level Alone/SmallGroup/Control factor built from this
# repo's real group_5 encoding (alone_2w+alone_4w pooled vs.
# smallgroup_2w+smallgroup_4w pooled vs. control) -- 23's own default
# (--solo-pattern/--group-pattern regexes against "Single"/"Group" values)
# was written for a different, unrelated study template and doesn't match
# any column in configs/participants.tsv, so the defaults here are
# study-129-specific instead (mirroring 34/36's --group-column approach).
#
# Differences from 23 (mirroring 29/30/33-35's thickness convention): raw mm
# THICKNESS (not log_volume), no etiv_z, age_z RANK-transformed (not raw
# scale). The 4-way interaction term is the primary test, FDR-corrected
# across every subfield tested, with targeted LOSO on any hit below the
# screen threshold -- exploratory multi-comparison sweep, same standing
# safeguards as every other exploratory script in this project.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to hippo_thickness_tidy.tsv"),
  make_option(c("-p", "--participants"), type="character", help="Path to participants TSV (subject_id, <group-column>, age, sex)"),
  make_option(c("-m", "--moderators"), type="character", help="Path to a moderators TSV with a CES-D/ADS column"),
  make_option(c("--cesd-col"), type="character", default="ads_score", help="Column name for the depression score [default %default]"),
  make_option(c("-o", "--outdir"), type="character", default="results/36_hippo_thickness_hem_time_group_cesd"),
  make_option(c("--roi-set"), type="character", default="presubiculum,subiculum,CA1,CA2_CA3"),
  make_option(c("--roi-column"), type="character", default="region"),
  make_option(c("--group-column"), type="character", default="group_5"),
  make_option(c("--alone-pattern"), type="character", default="^alone_"),
  make_option(c("--smallgroup-pattern"), type="character", default="^smallgroup_"),
  make_option(c("--control-value"), type="character", default="control"),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--loso-screen-threshold"), type="double", default=0.10),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (is.null(opt$moderators)) stop("--moderators is required")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$outdir, "models"), showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_col <- opt$`roi-column`
group_col <- opt$`group-column`

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$value <- suppressWarnings(as.numeric(tidy$value))
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file", roi_col))
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set")

agg <- aggregate(as.formula(paste("value ~ subject_id + session + hemisphere +", roi_col)), data=roi_dat, FUN=mean)
names(agg)[names(agg) == roi_col] <- "subfield"

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"

participants$group3 <- ifelse(grepl(opt$`alone-pattern`, participants$group), "Solo",
                        ifelse(grepl(opt$`smallgroup-pattern`, participants$group), "Group",
                        ifelse(participants$group == opt$`control-value`, "Control", NA_character_)))
n_unclassified <- sum(is.na(participants$group3))
if (n_unclassified) msg("Note: %d subjects didn't classify into Solo/Group/Control (dropped): %s\n",
                         n_unclassified, paste(unique(participants$group[is.na(participants$group3)]), collapse=", "))

moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!opt$`cesd-col` %in% names(moderators)) stop(sprintf("--cesd-col '%s' not found in --moderators", opt$`cesd-col`))
moderators <- moderators[, c("subject_id", opt$`cesd-col`)]
names(moderators)[2] <- "cesd_raw"
moderators$cesd_raw <- suppressWarnings(as.numeric(as.character(moderators$cesd_raw)))

meta <- merge(participants[, c("subject_id","group","age","sex","group3")], moderators, by="subject_id", all.x=TRUE)
agg <- merge(agg, meta, by="subject_id")
agg <- agg[!is.na(agg$cesd_raw) & !is.na(agg$group3), , drop=FALSE]
agg$hemisphere <- factor(agg$hemisphere, levels=sort(unique(agg$hemisphere)))
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$group3 <- factor(agg$group3, levels=c("Control","Solo","Group"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(rank(agg$age)))
agg$cesd_z <- as.numeric(scale(agg$cesd_raw))

msg("n subjects with CES-D/ADS available: %d\n", length(unique(agg$subject_id)))
msg("group3 sizes: %s\n", paste(names(table(unique(agg[,c("subject_id","group3")])$group3)),
                                  table(unique(agg[,c("subject_id","group3")])$group3), sep="=", collapse=", "))

base_formula <- value ~ hemisphere * time_f * group3 * cesd_z + age_z + sex

fit_term_p <- function(model, term) {
  at <- tryCatch(anova(model), error=function(e) NULL)
  if (!is.null(at) && term %in% rownames(at)) at[term, "Pr(>F)"] else NA_real_
}
fourway_term <- "hemisphere:time_f:group3:cesd_z"

summary_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  msg("Fitting HEM x TIME x GROUP x CES-D model for %s...\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_by_roi[[roi_name]] <- d_roi

  model <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d_roi, REML=TRUE), error=function(e) NULL)
  if (is.null(model)) { msg("  -> failed to fit\n"); next }

  sink(file.path(opt$outdir, "models", paste0(roi_name, "_model.txt")))
  cat("Subfield:", roi_name, "\n\n")
  print(summary(model))
  sink()

  p4 <- fit_term_p(model, fourway_term)
  msg("  4-way p = %s\n", if (is.na(p4)) "NA (term not estimable)" else sprintf("%.4f", p4))

  summary_rows[[roi_name]] <- data.frame(
    roi = roi_name, n_obs = nrow(d_roi), n_subjects = length(unique(d_roi$subject_id)),
    fourway_p = p4, stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
if (is.null(summary_df) || !nrow(summary_df)) stop("No subfield models were successfully fit")
summary_df$fourway_p_fdr <- p.adjust(summary_df$fourway_p, method="fdr")
write.csv(summary_df, file.path(opt$outdir, "summary.csv"), row.names=FALSE)
msg("\nWrote summary.csv (FDR-corrected across %d subfields)\n", nrow(summary_df))

# ------------------------------------------------------------------------
# Targeted LOSO for hits below screen threshold
# ------------------------------------------------------------------------
screen <- summary_df[!is.na(summary_df$fourway_p) & summary_df$fourway_p < opt$`loso-screen-threshold`, ]
if (nrow(screen)) {
  msg("\nRunning targeted LOSO for %d subfield(s) with p < %.2f...\n", nrow(screen), opt$`loso-screen-threshold`)
  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat(sprintf("Targeted LOSO for HEM x TIME x GROUP x CES-D, p < %.2f\n\n", opt$`loso-screen-threshold`), file=con)
  loso_all <- list()
  for (roi_name in screen$roi) {
    d_roi <- change_by_roi[[roi_name]]
    full_p <- screen$fourway_p[screen$roi == roi_name]
    subjects <- unique(d_roi$subject_id)
    loso_p <- sapply(subjects, function(s) {
      d_sub <- d_roi[d_roi$subject_id != s, , drop=FALSE]
      if (length(unique(d_sub$group3)) < 3) return(NA_real_)
      fit <- tryCatch(lmerTest::lmer(update(base_formula, ". ~ . + (1 | subject_id)"), data=d_sub, REML=TRUE), error=function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      fit_term_p(fit, fourway_term)
    })
    loso_p <- loso_p[!is.na(loso_p)]
    n_flips <- sum((loso_p < opt$alpha) != (full_p < opt$alpha))
    cat(sprintf("%s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance\n",
                roi_name, full_p, min(loso_p), max(loso_p), n_flips, length(loso_p)), file=con)
    loso_all[[roi_name]] <- data.frame(roi=roi_name, loso_p=loso_p)
  }
  write.csv(do.call(rbind, loso_all), file.path(opt$outdir, "loso_hits.csv"), row.names=FALSE)
  msg("Done. See %s/loso_summary.txt\n", opt$outdir)
} else {
  msg("\nNo subfield reached p < %.2f; nothing to LOSO-check.\n", opt$`loso-screen-threshold`)
}
