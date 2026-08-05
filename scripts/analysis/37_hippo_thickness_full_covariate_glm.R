#!/usr/bin/env Rscript
#
# Full-covariate repeated-measures GLM for hippocampal subfield THICKNESS,
# matching the colleague's own design (26_full_covariate_glm.R, volume):
# HEM x TIME x GROUP interaction (plus TIME x GROUP without the hemisphere
# breakdown), with GENDER and baseline CES-D as additional between-subjects
# factors/covariate, ALL interactions included (fully factorial, Type III
# SS), via jmv::anovaRM (not car::Anova -- see 26's docstring for why:
# car::Anova was verified NOT to reproduce the colleague's reported
# statistics on their own embedded volumes, jmv::anovaRM does).
#
# Differences from 26:
#   - THICKNESS (raw mm, no log-transform), not log(volume) -- matching
#     29/30/33-36's convention.
#   - GROUP3 (Solo/Group/Control) is built from this repo's real group_5
#     encoding (alone_2w+alone_4w / smallgroup_2w+smallgroup_4w / control)
#     via --group-column/--alone-pattern/--smallgroup-pattern/--control-value,
#     not 26's "Single"/"Group"/"Control" regexes against a literal "group"
#     column (written for a different, unrelated study template).
#   - --roi-column defaults to "region" (hipsta has no head/body split the
#     way FreeSurfer's aseg subfield stats do -- one mean thickness per
#     whole subfield per hemisphere/session), not 26's "subfield_raw".
#
# Correction: OFF by default (--fdr to enable), same rationale as 26 -- this
# script exists to mirror the colleague's own (uncorrected) approach as
# closely as possible for direct comparison, as well as support a more
# rigorous FDR-corrected re-analysis, without silently picking one for you.

suppressPackageStartupMessages({
  library(optparse)
  library(jmv)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_thickness_tidy.tsv (subject_id, session, hemisphere, region, value)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, <group-column>, age, sex"),
  make_option(c("-m", "--moderators"), type="character",
              help="Path to a moderators TSV with a baseline CES-D/ADS column"),
  make_option(c("--cesd-col"), type="character", default="ads_score",
              help="Column name for the baseline depression score in --moderators [default %default]"),
  make_option(c("-o", "--outdir"), type="character", default="results/37_hippo_thickness_full_covariate_glm",
              help="Output directory [default %default]"),
  make_option(c("--roi-column"), type="character", default="region",
              help="Column identifying the subfield in the tidy file [default %default]"),
  make_option(c("--roi-set"), type="character", default=NULL,
              help="Comma-separated subfield names to restrict to [default: all present in the tidy file]"),
  make_option(c("--group-column"), type="character", default="group_5",
              help="Column in --participants holding the group assignment [default %default]"),
  make_option(c("--alone-pattern"), type="character", default="^alone_",
              help="Regex matched against --group-column for the 'Solo' GROUP level [default %default]"),
  make_option(c("--smallgroup-pattern"), type="character", default="^smallgroup_",
              help="Regex matched against --group-column for the 'Group' GROUP level [default %default]"),
  make_option(c("--control-value"), type="character", default="control",
              help="Group value treated as the third (control) GROUP level [default %default]"),
  make_option(c("--fdr"), action="store_true", default=FALSE,
              help="Apply FDR correction across subfields (colleague's own approach used no correction) [default %default]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (is.null(opt$moderators)) stop("--moderators is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))
if (!file.exists(opt$moderators)) stop(sprintf("moderators file not found: %s", opt$moderators))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-column`
cesd_col <- opt$`cesd-col`
group_col <- opt$`group-column`

# ------------------------------------------------------------------------
# Load and prepare data
# ------------------------------------------------------------------------
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file (columns: %s)", roi_col, paste(names(tidy), collapse=", ")))
tidy$value <- suppressWarnings(as.numeric(tidy$value))

sessions <- sort(unique(tidy$session))
if (length(sessions) != 2) stop(sprintf("This design requires exactly 2 sessions (TIME pre/post); found: %s", paste(sessions, collapse=", ")))

has_hemisphere <- length(unique(tidy$hemisphere)) > 1
msg("Design: %s, sessions %s\n", if (has_hemisphere) "HEM x TIME x GROUP (+ TIME x GROUP)" else "TIME x GROUP only (midline structures, no hemisphere factor)",
    paste(sessions, collapse=" vs "))

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"

participants$group3 <- ifelse(grepl(opt$`alone-pattern`, participants$group), "Solo",
                        ifelse(grepl(opt$`smallgroup-pattern`, participants$group), "Group",
                        ifelse(participants$group == opt$`control-value`, "Control", NA_character_)))
if (any(is.na(participants$group3))) {
  warning(sprintf("Some subjects didn't classify into Solo/Group/Control -- check --alone-pattern/--smallgroup-pattern/--control-value: %s",
                   paste(unique(participants$group[is.na(participants$group3)]), collapse=", ")))
}

moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!cesd_col %in% names(moderators)) stop(sprintf("--cesd-col '%s' not found in --moderators", cesd_col))
meta <- merge(participants, moderators[, c("subject_id", cesd_col)], by="subject_id", all.x=TRUE)
names(meta)[names(meta) == cesd_col] <- "cesd_baseline"

roi_set <- if (!is.null(opt$`roi-set`)) trimws(strsplit(opt$`roi-set`, ",")[[1]]) else sort(unique(tidy[[roi_col]]))

# ------------------------------------------------------------------------
# Fully-factorial between-subjects term set: GROUP x SEX x baseline-CESD,
# all interactions -- matches 26/the colleague's design.
# ------------------------------------------------------------------------
bsTerms <- list("group3","sex","cesd_baseline",
                c("group3","sex"), c("group3","cesd_baseline"), c("sex","cesd_baseline"),
                c("group3","sex","cesd_baseline"))

fit_one_roi <- function(roi) {
  w <- tidy[tidy[[roi_col]] == roi, c("subject_id","session","hemisphere","value")]
  if (has_hemisphere) {
    w$key <- paste0(ifelse(w$hemisphere == unique(w$hemisphere)[1], "A", "B"), match(w$session, sessions))
  } else {
    w$key <- paste0("T", match(w$session, sessions))
  }
  piv <- reshape(w[, c("subject_id","key","value")], idvar="subject_id", timevar="key", direction="wide")
  names(piv) <- sub("^value\\.", "", names(piv))
  m <- merge(piv, meta[, c("subject_id","group3","sex","cesd_baseline")], by="subject_id")
  m <- m[complete.cases(m), , drop=FALSE]
  n <- nrow(m)
  if (n < 20) return(list(n=n, hem_time_group_F=NA, hem_time_group_p=NA, hem_time_group_pes=NA,
                           time_group_F=NA, time_group_p=NA, time_group_pes=NA))

  wide <- data.frame(group3=factor(m$group3), sex=factor(m$sex), cesd_baseline=m$cesd_baseline)
  if (has_hemisphere) {
    wide$A1 <- m$A1; wide$A2 <- m$A2; wide$B1 <- m$B1; wide$B2 <- m$B2
    res <- tryCatch(anovaRM(
      data = wide,
      rm = list(list(label="hemisphere", levels=c("A","B")), list(label="time", levels=c("1","2"))),
      rmCells = list(list(measure="A1", cell=c("A","1")), list(measure="A2", cell=c("A","2")),
                     list(measure="B1", cell=c("B","1")), list(measure="B2", cell=c("B","2"))),
      bs = c("group3","sex"), cov = c("cesd_baseline"), bsTerms = bsTerms,
      ss = "3", effectSize = c("partEta")
    ), error=function(e) NULL)
  } else {
    wide$T1 <- m$T1; wide$T2 <- m$T2
    res <- tryCatch(anovaRM(
      data = wide,
      rm = list(list(label="time", levels=c("1","2"))),
      rmCells = list(list(measure="T1", cell="1"), list(measure="T2", cell="2")),
      bs = c("group3","sex"), cov = c("cesd_baseline"), bsTerms = bsTerms,
      ss = "3", effectSize = c("partEta")
    ), error=function(e) NULL)
  }
  if (is.null(res)) return(list(n=n, hem_time_group_F=NA, hem_time_group_p=NA, hem_time_group_pes=NA,
                                 time_group_F=NA, time_group_p=NA, time_group_pes=NA))

  rmt <- as.data.frame(res$rmTable$asDF)
  names(rmt) <- sub("\\[none\\]$", "", names(rmt))
  get_row <- function(term) rmt[rmt$name == term, ]
  tg <- get_row("time:group3")
  if (has_hemisphere) {
    htg <- get_row("hemisphere:time:group3")
    list(n=n,
         hem_time_group_F = if (nrow(htg)) htg$F else NA, hem_time_group_p = if (nrow(htg)) htg$p else NA,
         hem_time_group_pes = if (nrow(htg)) htg$partEta else NA,
         time_group_F = if (nrow(tg)) tg$F else NA, time_group_p = if (nrow(tg)) tg$p else NA,
         time_group_pes = if (nrow(tg)) tg$partEta else NA)
  } else {
    list(n=n, hem_time_group_F=NA, hem_time_group_p=NA, hem_time_group_pes=NA,
         time_group_F = if (nrow(tg)) tg$F else NA, time_group_p = if (nrow(tg)) tg$p else NA,
         time_group_pes = if (nrow(tg)) tg$partEta else NA)
  }
}

summary_rows <- list()
for (roi in roi_set) {
  msg("Fitting subfield: %s\n", roi)
  r <- fit_one_roi(roi)
  summary_rows[[roi]] <- data.frame(roi=roi, n=r$n,
                                     hem_time_group_F=r$hem_time_group_F, hem_time_group_p=r$hem_time_group_p, hem_time_group_pes=r$hem_time_group_pes,
                                     time_group_F=r$time_group_F, time_group_p=r$time_group_p, time_group_pes=r$time_group_pes,
                                     stringsAsFactors=FALSE)
}
summary_df <- do.call(rbind, summary_rows)

if (isTRUE(opt$fdr)) {
  summary_df$hem_time_group_p_fdr <- p.adjust(summary_df$hem_time_group_p, method="fdr")
  summary_df$time_group_p_fdr <- p.adjust(summary_df$time_group_p, method="fdr")
  summary_df$significant_hem_time_group <- summary_df$hem_time_group_p_fdr < opt$alpha
  summary_df$significant_time_group <- summary_df$time_group_p_fdr < opt$alpha
} else {
  summary_df$significant_hem_time_group <- summary_df$hem_time_group_p < opt$alpha
  summary_df$significant_time_group <- summary_df$time_group_p < opt$alpha
  summary_df$trend_hem_time_group <- summary_df$hem_time_group_p < 0.10
  summary_df$trend_time_group <- summary_df$time_group_p < 0.10
}

summary_df <- summary_df[order(summary_df$time_group_p), ]
write.csv(summary_df, file.path(opt$outdir, "full_covariate_glm_summary.csv"), row.names=FALSE)
msg("\nWrote %s across %d subfields to full_covariate_glm_summary.csv (FDR %s)\n",
    if (has_hemisphere) "HEM x TIME x GROUP and TIME x GROUP" else "TIME x GROUP",
    nrow(summary_df), if (isTRUE(opt$fdr)) "applied" else "NOT applied (matches colleague's own approach)")
