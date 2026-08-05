#!/usr/bin/env Rscript
#
# Targeted, pre-specified confirmatory test: does baseline TSD-Z (creativity,
# Test zum schoepferischen Denken - zeichnerisch) relate to structural brain
# change, for the two domains (brainstem, hippocampal thickness) flagged by
# 38_unsupervised_change_clustering.R's EXPLORATORY post-hoc cluster
# profiling? This script exists specifically to follow up that lead properly
# rather than treat it as confirmed or dismiss it as noise.
#
# Why ANCOVA (final ~ baseline + moderator), not a raw change score
# (final - baseline ~ moderator): 38's own follow-up check showed the
# baseline-TSD-Z and change-TSD-Z cluster associations in both domains have
# the signature of regression-to-the-mean in TSD-Z itself (the higher-
# baseline-TSD-Z cluster mechanically declines more) -- using baseline as a
# COVARIATE rather than differencing it out is the standard fix (Lord's
# paradox / ANCOVA-vs-change-score literature): it tests whether baseline
# creativity predicts the follow-up brain measure over and above where that
# subject started, without manufacturing a spurious baseline-change
# correlation the way a raw difference score does.
#
# Two complementary tests, both controlling for age_z (rank-transformed),
# sex, and intervention group (to net out the creativity effect from group
# assignment, not test group itself):
#   1. Per-ROI ANCOVA: final ~ baseline + moderator_z + age_z + sex + intervention,
#      one test per ROI in the domain, FDR-corrected across those ROIs.
#   2. Composite: the SAME residualized-change PC1 that drove the original
#      unsupervised cluster split (recomputed here, not re-used from 38, so
#      this script is self-contained) tested as
#      PC1 ~ moderator_z + age_z + sex + intervention -- the single test that
#      most directly asks "does creativity relate to the pattern that was
#      actually flagged", with no per-ROI multiplicity.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to a domain tidy TSV"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-m", "--moderators"), type="character",
              default="/data/local/129_PK01/moderators_baseline.tsv",
              help="Path to baseline moderators TSV [default %default]"),
  make_option(c("--moderator-col"), type="character", default="tsdz_total",
              help="Column in --moderators to test [default %default]"),
  make_option(c("--roi-col"), type="character", default="subfield"),
  make_option(c("--value-col"), type="character", default="volume"),
  make_option(c("--roi-set"), type="character", default="all"),
  make_option(c("--hemisphere-agg"), type="character", default="sum", help="'sum' or 'mean' [default %default]"),
  make_option(c("--log-transform"), action="store_true", default=FALSE,
              help="Model log(value) instead of raw value (volume domains) [default FALSE]"),
  make_option(c("--intervention-groups"), type="character", default="Single_2wk,Single_4wk,Group_2wk,Group_4wk"),
  make_option(c("--control-group"), type="character", default="Control"),
  make_option(c("-o", "--outdir"), type="character", default="results/39_creativity_brain_ancova"),
  make_option(c("--label"), type="character", default="domain"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--final-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-col`; value_col <- opt$`value-col`
intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)
mod_col <- opt$`moderator-col`

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy[[value_col]] <- suppressWarnings(as.numeric(tidy[[value_col]]))
roi_set <- if (identical(opt$`roi-set`, "all")) sort(unique(tidy[[roi_col]])) else trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]

agg_fun <- if (identical(opt$`hemisphere-agg`, "mean")) mean else sum
agg <- aggregate(as.formula(paste0(value_col, " ~ subject_id + session + `", roi_col, "`")), data=roi_dat, FUN=agg_fun)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("[%s] baseline=%s, final=%s, moderator=%s\n", opt$label, baseline_ses, final_ses, mod_col)

base_vals <- agg[agg$session == baseline_ses, c("subject_id", roi_col, value_col)]; names(base_vals)[3] <- "baseline"
final_vals <- agg[agg$session == final_ses, c("subject_id", roi_col, value_col)]; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id", roi_col))
if (isTRUE(opt$`log-transform`)) {
  merged_long$baseline <- log(merged_long$baseline)
  merged_long$final <- log(merged_long$final)
}
merged_long$change <- merged_long$final - merged_long$baseline

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
participants$intervention <- factor(ifelse(participants$group %in% intervention_groups, "intervention",
                                     ifelse(participants$group == control_group, "control", NA)), levels=c("control","intervention"))
participants$sex <- factor(participants$sex)
participants$age_z <- as.numeric(scale(rank(participants$age)))

moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!mod_col %in% names(moderators)) stop(sprintf("--moderator-col '%s' not found in --moderators", mod_col))
moderators <- moderators[, c("subject_id", mod_col)]
names(moderators)[2] <- "moderator_raw"

meta <- merge(participants[, c("subject_id","intervention","age_z","sex")], moderators, by="subject_id")
meta <- meta[!is.na(meta$intervention) & !is.na(meta$moderator_raw), , drop=FALSE]
meta$moderator_z <- as.numeric(scale(meta$moderator_raw))

# ------------------------------------------------------------------------
# 1. Per-ROI ANCOVA: final ~ baseline + moderator_z + age_z + sex + intervention
# ------------------------------------------------------------------------
roi_rows <- list()
for (roi_name in roi_set) {
  d <- merged_long[merged_long[[roi_col]] == roi_name, , drop=FALSE]
  d <- merge(d, meta, by="subject_id")
  d <- d[complete.cases(d[, c("baseline","final","moderator_z","age_z","sex","intervention")]), , drop=FALSE]
  if (nrow(d) < 15) { msg("[%s] %s: only %d complete subjects, skipping\n", opt$label, roi_name, nrow(d)); next }
  fit <- lm(final ~ baseline + moderator_z + age_z + sex + intervention, data=d)
  tt <- summary(fit)$coefficients
  roi_rows[[roi_name]] <- data.frame(
    roi = roi_name, n = nrow(d),
    moderator_beta = tt["moderator_z", "Estimate"],
    moderator_se = tt["moderator_z", "Std. Error"],
    moderator_t = tt["moderator_z", "t value"],
    moderator_p = tt["moderator_z", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}
roi_df <- do.call(rbind, roi_rows)
if (!is.null(roi_df) && nrow(roi_df)) {
  roi_df$moderator_p_fdr <- p.adjust(roi_df$moderator_p, method="fdr")
  roi_df$significant_uncorrected <- roi_df$moderator_p < opt$alpha
  roi_df$significant_fdr <- roi_df$moderator_p_fdr < opt$alpha
  roi_df <- roi_df[order(roi_df$moderator_p), ]
  write.csv(roi_df, file.path(opt$outdir, paste0(opt$label, "_per_roi_ancova.csv")), row.names=FALSE)
}

# ------------------------------------------------------------------------
# 2. Composite: residualized-change PC1 (same construction as 38) ~ moderator_z
# ------------------------------------------------------------------------
change_wide <- reshape(merged_long[, c("subject_id", roi_col, "change")], idvar="subject_id", timevar=roi_col, direction="wide")
names(change_wide) <- sub("^change\\.", "", names(change_wide))
change_cols <- setdiff(names(change_wide), "subject_id")

dat <- merge(change_wide, meta, by="subject_id")
dat <- dat[complete.cases(dat[, c(change_cols, "moderator_z","age_z","sex","intervention")]), , drop=FALSE]

resid_mat <- sapply(change_cols, function(col) residuals(lm(as.formula(paste0("`", col, "` ~ age_z + sex")), data=dat)))
colnames(resid_mat) <- change_cols
pca_fit <- prcomp(resid_mat, center=TRUE, scale.=TRUE)
dat$PC1 <- pca_fit$x[, 1]
var1 <- pca_fit$sdev[1]^2 / sum(pca_fit$sdev^2)

fit_pc1 <- lm(PC1 ~ moderator_z + age_z + sex + intervention, data=dat)
tt_pc1 <- summary(fit_pc1)$coefficients

con <- file(file.path(opt$outdir, paste0(opt$label, "_summary.txt")), open="wt")
cat(sprintf(
  "Targeted confirmatory test: baseline %s vs. %s structural change\nBaseline -> final: %s -> %s. n(ANCOVA composite) = %d\n\n1. COMPOSITE test (primary): residualized-change PC1 (%.1f%% of variance) ~ %s_z + age_z + sex + intervention\n   %s beta = %.4f, SE = %.4f, t = %.2f, p = %.4f\n\n2. Per-ROI ANCOVA (final ~ baseline + %s_z + age_z + sex + intervention), FDR across %d ROIs:\n",
  mod_col, opt$label, baseline_ses, final_ses, nrow(dat),
  100*var1, mod_col,
  mod_col, tt_pc1["moderator_z","Estimate"], tt_pc1["moderator_z","Std. Error"], tt_pc1["moderator_z","t value"], tt_pc1["moderator_z","Pr(>|t|)"],
  mod_col, if (!is.null(roi_df)) nrow(roi_df) else 0
), file=con)
if (!is.null(roi_df) && nrow(roi_df)) {
  for (i in seq_len(nrow(roi_df))) {
    row <- roi_df[i, ]
    cat(sprintf("   %-30s beta=%.4f p=%.4f p_fdr=%.4f%s\n", row$roi, row$moderator_beta, row$moderator_p, row$moderator_p_fdr,
                if (row$significant_fdr) " *" else if (row$significant_uncorrected) " (nominal)" else ""), file=con)
  }
}
close(con)

msg("[%s] Composite PC1 ~ %s_z: beta=%.4f, p=%.4f\n", opt$label, mod_col, tt_pc1["moderator_z","Estimate"], tt_pc1["moderator_z","Pr(>|t|)"])
msg("[%s] Done. See %s/%s_summary.txt\n", opt$label, opt$outdir, opt$label)
