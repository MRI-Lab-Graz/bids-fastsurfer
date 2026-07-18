#!/usr/bin/env Rscript
#
# Multivariate pattern analysis: does dance vs control differ in the JOINT
# pattern of change across the pre-specified ROI set, even if no single ROI
# survives per-ROI FDR correction?
#
# Two complementary, non-black-box approaches (deliberately NOT machine
# learning classification -- at this sample size a MANOVA/PCA is more
# trustworthy than a classifier, see 04_exploratory_rf.R's own caveats):
#
#   1. MANOVA on the vector of per-subject, per-ROI change scores
#      (baseline -> final session). One multivariate test (Pillai's trace)
#      for "does the pattern across ROIs differ by group", instead of one
#      FDR-corrected test per ROI. More powerful when the true effect is a
#      coordinated shift distributed across several regions rather than one
#      large hit in a single region.
#
#   2. PCA on the (covariate-residualized) change-score matrix to find
#      latent components of COORDINATED change across regions, then test
#      those components by group. Components are NOT assumed to load in the
#      same direction on every ROI -- a region-specific pattern (some ROIs
#      up, others down) is exactly what PCA can capture and a naive
#      "more volume is better" composite/average score would erase.
#
# Change scores are summed across hemisphere per subfield (one value per
# subject per ROI) to keep the feature count manageable for MANOVA/PCA at
# this sample size; see --keep-hemisphere-separate to disable that.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/06_manova_pca_pattern",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP",
              help="Comma-separated pre-specified subfield names [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value to use as baseline [default: earliest session present]"),
  make_option(c("--final-session"), type="character", default=NULL,
              help="Session value to use as endpoint [default: latest session present]"),
  make_option(c("--keep-hemisphere-separate"), action="store_true", default=FALSE,
              help="Keep lh/rh as separate features instead of summing (doubles feature count)"),
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

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load data and build subject x ROI change-score matrix
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check subfield names against the tidy file")

group_cols <- if (isTRUE(opt$`keep-hemisphere-separate`)) c("subject_id","session","hemisphere","subfield") else c("subject_id","session","subfield")
agg <- aggregate(volume ~ ., data=roi_dat[, c(group_cols, "volume")], FUN=sum)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Change score: %s -> %s\n", baseline_ses, final_ses)

feature_col <- if (isTRUE(opt$`keep-hemisphere-separate`)) {
  agg$feature <- paste(agg$hemisphere, agg$subfield, sep="_")
  "feature"
} else {
  "subfield"
}

base_vals <- agg[agg$session == baseline_ses, c("subject_id", feature_col, "volume")]
final_vals <- agg[agg$session == final_ses, c("subject_id", feature_col, "volume")]
names(base_vals)[3] <- "baseline"
names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id", feature_col))
merged_long$change <- log(merged_long$final) - log(merged_long$baseline)

change_wide <- reshape(merged_long[, c("subject_id", feature_col, "change")],
                        idvar="subject_id", timevar=feature_col, direction="wide")
change_cols <- setdiff(names(change_wide), "subject_id")
names(change_wide)[names(change_wide) %in% change_cols] <- sub("^change\\.", "", change_cols)
change_cols <- setdiff(names(change_wide), "subject_id")

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

dat <- merge(change_wide, participants[, required_pcols], by="subject_id")
dat <- dat[dat$group %in% c(dance_groups, control_group), , drop=FALSE]
dat$dance <- factor(ifelse(dat$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
dat$sex <- factor(dat$sex)
dat$age_z <- as.numeric(scale(dat$age))

n_before <- nrow(dat)
dat <- dat[complete.cases(dat[, c(change_cols, "dance","age_z","sex")]), , drop=FALSE]
msg("Subjects with complete baseline+final data across all %d ROIs: %d of %d\n", length(change_cols), nrow(dat), n_before)
if (nrow(dat) < length(change_cols) + 10) {
  warning(sprintf("Only %d complete subjects for %d ROI features -- MANOVA/PCA results at this ratio should be treated cautiously", nrow(dat), length(change_cols)))
}

write.csv(dat[, c("subject_id","dance",change_cols)], file.path(opt$outdir, "change_score_matrix.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# 1. MANOVA: joint test of the ROI-change pattern by group
# ------------------------------------------------------------------------
msg("Running MANOVA on joint ROI-change pattern...\n")
Y <- as.matrix(dat[, change_cols])
manova_fit <- manova(Y ~ dance + age_z + sex, data=dat)
manova_summary <- as.data.frame(summary(manova_fit, test="Pillai")$stats)
manova_summary$term <- rownames(manova_summary)
write.csv(manova_summary, file.path(opt$outdir, "manova_results.csv"), row.names=FALSE)

dance_manova_p <- manova_summary[manova_summary$term == "dance", "Pr(>F)"]
msg("MANOVA (Pillai's trace) for 'dance' across %d ROIs jointly: p = %s\n",
    length(change_cols), if(length(dance_manova_p)) format(dance_manova_p, digits=4) else "NA")

# ------------------------------------------------------------------------
# 2. PCA on covariate-residualized change scores
# ------------------------------------------------------------------------
msg("Running PCA on covariate-residualized change-score matrix...\n")
resid_mat <- sapply(change_cols, function(col) {
  f <- as.formula(paste0("`", col, "` ~ age_z + sex"))
  residuals(lm(f, data=dat))
})
colnames(resid_mat) <- change_cols

pca_fit <- prcomp(resid_mat, center=TRUE, scale.=TRUE)
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
n_keep <- max(1, sum(pca_fit$sdev^2 > 1))  # Kaiser criterion (eigenvalue > 1)
msg("PCA: %d components retained by Kaiser criterion (eigenvalue > 1) out of %d\n", n_keep, length(change_cols))

loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop=FALSE])
loadings_df$roi <- rownames(loadings_df)
write.csv(loadings_df, file.path(opt$outdir, "pca_loadings.csv"), row.names=FALSE)
write.csv(data.frame(component=paste0("PC", seq_along(var_explained)), variance_explained=var_explained),
          file.path(opt$outdir, "pca_variance_explained.csv"), row.names=FALSE)

pc_scores <- as.data.frame(pca_fit$x[, seq_len(n_keep), drop=FALSE])
pc_scores$dance <- dat$dance

pc_test_rows <- list()
for (i in seq_len(n_keep)) {
  pc_name <- paste0("PC", i)
  fit <- lm(pc_scores[[pc_name]] ~ dance, data=pc_scores)
  tt <- summary(fit)$coefficients
  pc_test_rows[[pc_name]] <- data.frame(
    component = pc_name,
    variance_explained = var_explained[i],
    dance_estimate = if ("dancedance" %in% rownames(tt)) tt["dancedance", "Estimate"] else NA_real_,
    dance_p = if ("dancedance" %in% rownames(tt)) tt["dancedance", "Pr(>|t|)"] else NA_real_,
    stringsAsFactors = FALSE
  )
}
pc_test_df <- do.call(rbind, pc_test_rows)
pc_test_df$dance_p_fdr <- p.adjust(pc_test_df$dance_p, method="fdr")
pc_test_df$significant <- pc_test_df$dance_p_fdr < opt$alpha
write.csv(pc_test_df, file.path(opt$outdir, "pca_component_group_tests.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------
con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf(
  "Multivariate pattern analysis: %s -> %s change scores across %d ROIs\n\nSubjects: %d\n\n1. MANOVA (Pillai's trace): joint test of dance vs control across the ROI pattern\n   p = %s\n\n2. PCA (Kaiser criterion, eigenvalue > 1): %d components retained, %.1f%% variance explained\n   Per-component group tests (FDR-corrected across components):\n",
  baseline_ses, final_ses, length(change_cols), nrow(dat),
  if(length(dance_manova_p)) format(dance_manova_p, digits=4) else "NA",
  n_keep, 100*sum(var_explained[seq_len(n_keep)])
), file=con)
for (i in seq_len(nrow(pc_test_df))) {
  row <- pc_test_df[i, ]
  cat(sprintf("   - %s (%.1f%% var): dance effect = %.4f, p = %.4f, p_fdr = %.4f%s\n",
              row$component, 100*row$variance_explained, row$dance_estimate, row$dance_p, row$dance_p_fdr,
              if (row$significant) " *" else ""), file=con)
}
cat("\nSee pca_loadings.csv for which ROIs drive each component (same-sign loadings\n",
    "= coordinated same-direction change; mixed-sign loadings = a contrast pattern\n",
    "between regions, e.g. some growing while others shrink).\n", file=con)

msg("Done. See %s/summary.txt for a plain-text overview.\n", opt$outdir)
