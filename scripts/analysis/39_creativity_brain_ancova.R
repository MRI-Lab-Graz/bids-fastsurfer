#!/usr/bin/env Rscript
#
# Targeted ANCOVA test, generalized over its primary variable of interest:
#
#   --primary-variable moderator (original use): does a baseline
#     psychological moderator (default: TSD-Z creativity) relate to
#     structural brain change, net of intervention group? Built to follow
#     up the TSD-Z lead flagged by 38_unsupervised_change_clustering.R's
#     exploratory post-hoc cluster profiling.
#
#   --primary-variable intervention: does INTERVENTION GROUP itself relate
#     to structural brain change? This re-tests the project's original
#     confirmatory question (already covered by 01/10/24/25/27's
#     LMM-on-change-scores design) but with the ANCOVA design instead
#     (final ~ baseline + intervention + age + sex) as a cross-check using
#     a different, arguably more robust-to-regression-to-the-mean
#     statistical framework -- baseline enters as an explicit covariate
#     rather than being differenced out or handled via a time*group
#     interaction term.
#
# Why ANCOVA (final ~ baseline + X), not a raw change score
# (final - baseline ~ X), in either mode: differencing out baseline
# manufactures a spurious baseline-change correlation for ANY predictor X
# whenever the outcome has imperfect test-retest reliability (Lord's
# paradox) -- using baseline as a covariate instead of differencing it out
# is the standard fix and was needed here because 38's own follow-up check
# showed exactly this regression-to-the-mean signature for TSD-Z.
#
# Two complementary tests, both controlling for age_z (rank-transformed)
# and sex, plus (in --primary-variable moderator mode) intervention group
# as a nuisance covariate to net the moderator effect out of group
# assignment:
#   1. Per-ROI ANCOVA: final ~ baseline + <primary> + age_z + sex [+ intervention],
#      one test per ROI in the domain, FDR-corrected across those ROIs.
#   2. Composite: residualized-change PC1 (covariates: age_z + sex) tested as
#      PC1 ~ <primary> + age_z + sex [+ intervention] -- the single test that
#      most directly asks "does <primary> relate to the dominant pattern of
#      change", with no per-ROI multiplicity.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to a domain tidy TSV"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("--primary-variable"), type="character", default="moderator",
              help="'moderator' (test --moderator-col, net of intervention) or 'intervention' (test intervention group itself) [default %default]"),
  make_option(c("-m", "--moderators"), type="character",
              default="/data/local/129_PK01/moderators_baseline.tsv",
              help="Path to baseline moderators TSV (only used when --primary-variable=moderator) [default %default]"),
  make_option(c("--moderator-col"), type="character", default="tsdz_total",
              help="Column in --moderators to test (only used when --primary-variable=moderator) [default %default]"),
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
if (!opt$`primary-variable` %in% c("moderator","intervention")) stop("--primary-variable must be 'moderator' or 'intervention'")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-col`; value_col <- opt$`value-col`
intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)
mode_moderator <- identical(opt$`primary-variable`, "moderator")
mod_col <- opt$`moderator-col`
# the term name as it appears in lm()'s coefficient table, and the label used
# in messages/output files
primary_term <- if (mode_moderator) "moderator_z" else "interventionintervention"
primary_label <- if (mode_moderator) mod_col else "intervention"

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy[[value_col]] <- suppressWarnings(as.numeric(tidy[[value_col]]))
roi_set <- if (identical(opt$`roi-set`, "all")) sort(unique(tidy[[roi_col]])) else trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]

agg_fun <- if (identical(opt$`hemisphere-agg`, "mean")) mean else sum
agg <- aggregate(as.formula(paste0(value_col, " ~ subject_id + session + `", roi_col, "`")), data=roi_dat, FUN=agg_fun)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("[%s] baseline=%s, final=%s, primary variable=%s\n", opt$label, baseline_ses, final_ses, primary_label)

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

if (mode_moderator) {
  moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if (!mod_col %in% names(moderators)) stop(sprintf("--moderator-col '%s' not found in --moderators", mod_col))
  moderators <- moderators[, c("subject_id", mod_col)]
  names(moderators)[2] <- "moderator_raw"
  meta <- merge(participants[, c("subject_id","intervention","age_z","sex")], moderators, by="subject_id")
  meta <- meta[!is.na(meta$intervention) & !is.na(meta$moderator_raw), , drop=FALSE]
  meta$moderator_z <- as.numeric(scale(meta$moderator_raw))
  roi_formula_rhs <- "baseline + moderator_z + age_z + sex + intervention"
  pc1_formula_rhs <- "moderator_z + age_z + sex + intervention"
} else {
  meta <- participants[!is.na(participants$intervention), c("subject_id","intervention","age_z","sex")]
  roi_formula_rhs <- "baseline + intervention + age_z + sex"
  pc1_formula_rhs <- "intervention + age_z + sex"
}

# ------------------------------------------------------------------------
# 1. Per-ROI ANCOVA
# ------------------------------------------------------------------------
roi_rows <- list()
for (roi_name in roi_set) {
  d <- merged_long[merged_long[[roi_col]] == roi_name, , drop=FALSE]
  d <- merge(d, meta, by="subject_id")
  required_cols <- c("baseline","final","age_z","sex","intervention", if (mode_moderator) "moderator_z" else NULL)
  d <- d[complete.cases(d[, required_cols]), , drop=FALSE]
  if (nrow(d) < 15) { msg("[%s] %s: only %d complete subjects, skipping\n", opt$label, roi_name, nrow(d)); next }
  fit <- lm(as.formula(paste("final ~", roi_formula_rhs)), data=d)
  tt <- summary(fit)$coefficients
  if (!primary_term %in% rownames(tt)) { msg("[%s] %s: term '%s' not estimable, skipping\n", opt$label, roi_name, primary_term); next }
  roi_rows[[roi_name]] <- data.frame(
    roi = roi_name, n = nrow(d),
    primary_beta = tt[primary_term, "Estimate"],
    primary_se = tt[primary_term, "Std. Error"],
    primary_t = tt[primary_term, "t value"],
    primary_p = tt[primary_term, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}
roi_df <- do.call(rbind, roi_rows)
if (!is.null(roi_df) && nrow(roi_df)) {
  roi_df$primary_p_fdr <- p.adjust(roi_df$primary_p, method="fdr")
  roi_df$significant_uncorrected <- roi_df$primary_p < opt$alpha
  roi_df$significant_fdr <- roi_df$primary_p_fdr < opt$alpha
  roi_df <- roi_df[order(roi_df$primary_p), ]
  write.csv(roi_df, file.path(opt$outdir, paste0(opt$label, "_per_roi_ancova.csv")), row.names=FALSE)
}

# ------------------------------------------------------------------------
# 2. Composite: residualized-change PC1 ~ primary variable
# ------------------------------------------------------------------------
change_wide <- reshape(merged_long[, c("subject_id", roi_col, "change")], idvar="subject_id", timevar=roi_col, direction="wide")
names(change_wide) <- sub("^change\\.", "", names(change_wide))
change_cols <- setdiff(names(change_wide), "subject_id")

dat <- merge(change_wide, meta, by="subject_id")
required_cols2 <- c(change_cols, "age_z","sex","intervention", if (mode_moderator) "moderator_z" else NULL)
dat <- dat[complete.cases(dat[, required_cols2]), , drop=FALSE]

resid_mat <- sapply(change_cols, function(col) residuals(lm(as.formula(paste0("`", col, "` ~ age_z + sex")), data=dat)))
colnames(resid_mat) <- change_cols
pca_fit <- prcomp(resid_mat, center=TRUE, scale.=TRUE)
dat$PC1 <- pca_fit$x[, 1]
var1 <- pca_fit$sdev[1]^2 / sum(pca_fit$sdev^2)

fit_pc1 <- lm(as.formula(paste("PC1 ~", pc1_formula_rhs)), data=dat)
tt_pc1 <- summary(fit_pc1)$coefficients

con <- file(file.path(opt$outdir, paste0(opt$label, "_summary.txt")), open="wt")
cat(sprintf(
  "Targeted ANCOVA test (primary variable: %s): %s structural change\nBaseline -> final: %s -> %s. n(ANCOVA composite) = %d\n\n1. COMPOSITE test (primary): residualized-change PC1 (%.1f%% of variance) ~ %s\n   %s beta = %.4f, SE = %.4f, t = %.2f, p = %.4f\n\n2. Per-ROI ANCOVA (final ~ %s), FDR across %d ROIs:\n",
  primary_label, opt$label, baseline_ses, final_ses, nrow(dat),
  100*var1, pc1_formula_rhs,
  primary_label, tt_pc1[primary_term,"Estimate"], tt_pc1[primary_term,"Std. Error"], tt_pc1[primary_term,"t value"], tt_pc1[primary_term,"Pr(>|t|)"],
  roi_formula_rhs, if (!is.null(roi_df)) nrow(roi_df) else 0
), file=con)
if (!is.null(roi_df) && nrow(roi_df)) {
  for (i in seq_len(nrow(roi_df))) {
    row <- roi_df[i, ]
    cat(sprintf("   %-30s beta=%.4f p=%.4f p_fdr=%.4f%s\n", row$roi, row$primary_beta, row$primary_p, row$primary_p_fdr,
                if (row$significant_fdr) " *" else if (row$significant_uncorrected) " (nominal)" else ""), file=con)
  }
}
close(con)

msg("[%s] Composite PC1 ~ %s: beta=%.4f, p=%.4f\n", opt$label, primary_label, tt_pc1[primary_term,"Estimate"], tt_pc1[primary_term,"Pr(>|t|)"])
msg("[%s] Done. See %s/%s_summary.txt\n", opt$label, opt$outdir, opt$label)
