#!/usr/bin/env Rscript
#
# Multivariate pattern analysis on hippocampal subfield THICKNESS, run
# SEPARATELY within each subfield (presubiculum, subiculum, CA1, CA2/CA3),
# using each subfield's individual GRID POINTS as the feature set rather than
# a single subfield-mean per subject (as in 29_hippo_thickness_lmm.R /
# 30_hippo_thickness_factorial.R). Mirrors
# 28_cortical_cluster_multivariate.R's per-cluster design exactly, with one
# structural difference explained below.
#
# The question is not "does subfield X's mean thickness change" (29/30
# already answer that) but "does subfield X show a COORDINATED, spatially
# patterned change across its grid points" -- a real effect could be
# confined to one end of a subfield and get averaged away in the mean.
# This is only possible here because hipsta's grid points are in point
# correspondence across subjects by construction (see
# scripts/build_hipsta_grid_subfield_map.py's docstring) -- there is no
# equivalent for Destrieux regions, which is why 28 clusters together several
# whole ROIs instead of decomposing a single ROI into points.
#
# Same three complementary methods as 06_manova_pca_pattern.R /
# 28_cortical_cluster_multivariate.R, applied per subfield:
#   1. MANOVA (Pillai's trace) on the subject x grid-point change-score matrix
#      (baseline -> final session), covarying age_z (rank-transformed) + sex.
#   2. PCA on the covariate-residualized change-score matrix -- latent
#      components of coordinated change, tested against intervention.
#   3. PLS-DA (leave-one-out cross-validated) -- out-of-sample check of
#      whether the subfield's spatial pattern discriminates intervention from
#      control.
#
# lh/rh are kept SEPARATE, not averaged (unlike 28's cortical regions): a
# grid point's (x,y) index is only comparable within the same hemisphere
# (the map from scripts/build_hipsta_grid_subfield_map.py is per-hemisphere),
# so pooling lh/rh grid points would silently mix two different anatomical
# locations that happen to share an (x,y) index.
#
# MANOVA p-values are FDR-corrected across (subfield x hemisphere) tests.

suppressPackageStartupMessages({
  library(optparse)
  library(pls)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_thickness_grid_tidy.tsv (from extract_hipsta_thickness.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, <group-column>, age, sex"),
  make_option(c("-c", "--clusters"), type="character",
              default="configs/hipsta_grid_subfields.tsv",
              help="Path to grid_id -> subfield mapping TSV (from scripts/build_hipsta_grid_subfield_map.py) [default %default]"),
  make_option(c("-o", "--outdir"), type="character", default="results/31_hippo_thickness_cluster_multivariate",
              help="Output directory [default %default]"),
  make_option(c("--group-column"), type="character", default="group_5",
              help="Column in --participants holding the group assignment [default %default]"),
  make_option(c("--intervention-groups"), type="character",
              default="alone_2w,alone_4w,smallgroup_2w,smallgroup_4w",
              help="Comma-separated group values pooled into the 'intervention' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value to use as baseline [default: earliest session present]"),
  make_option(c("--final-session"), type="character", default=NULL,
              help="Session value to use as endpoint [default: latest session present]"),
  make_option(c("--min-grid-points"), type="integer", default=15,
              help="Minimum grid points a subfield/hemisphere must have (post grid-map filtering) to be analyzed [default %default]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Significance threshold [default %default]"),
  make_option(c("--seed"), type="integer", default=129, help="RNG seed for reproducibility [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }
set.seed(opt$seed)

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))
if (!file.exists(opt$clusters)) stop(sprintf("grid-subfield map not found: %s (run scripts/build_hipsta_grid_subfield_map.py first)", opt$clusters))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)
group_col <- opt$`group-column`

# ------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------
msg("Loading grid-point tidy thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$thickness <- suppressWarnings(as.numeric(tidy$thickness))
tidy$grid_id <- paste0(tidy$hemisphere, "_x", tidy$x, "_y", tidy$y)

grid_map <- read.delim(opt$clusters, header=TRUE, sep="\t", stringsAsFactors=FALSE)
msg("Grid-subfield map: %d grid points across %d hemisphere/subfield groups\n",
    nrow(grid_map), length(unique(paste(grid_map$hemisphere, grid_map$subfield))))

# tidy's own "subfield" column is the PER-SUBJECT hsf label at that grid
# point (from extract_hipsta_thickness.py); drop it before merging in the
# grid_map's cohort-level majority-vote "subfield" -- the whole point of
# grid_map is to use a fixed, cohort-consistent assignment (see
# build_hipsta_grid_subfield_map.py's docstring) instead of each subject's
# own wobbling boundary.
tidy$subfield <- NULL
tidy <- merge(tidy, grid_map[, c("grid_id","subfield")], by="grid_id")

sessions <- sort(unique(tidy$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Change score: %s -> %s\n", baseline_ses, final_ses)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"
participants$intervention <- factor(ifelse(participants$group %in% intervention_groups, "intervention",
                                     ifelse(participants$group == control_group, "control", NA)), levels=c("control","intervention"))
participants$sex <- factor(participants$sex)
participants$age_z <- as.numeric(scale(rank(participants$age)))

# grid-point thickness per subject/session/hemisphere/subfield/grid_id
# (already one row per key in the input; no aggregation needed)
base_vals <- tidy[tidy$session == baseline_ses, c("subject_id","hemisphere","subfield","grid_id","thickness")]
names(base_vals)[5] <- "baseline"
final_vals <- tidy[tidy$session == final_ses, c("subject_id","hemisphere","subfield","grid_id","thickness")]
names(final_vals)[5] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id","hemisphere","subfield","grid_id"))
merged_long$change <- merged_long$final - merged_long$baseline

# ------------------------------------------------------------------------
# Per-(hemisphere x subfield) MANOVA + PCA + PLS-DA
# ------------------------------------------------------------------------
run_group <- function(hemi, subfield_name, grid_ids) {
  label <- sprintf("%s_%s", hemi, subfield_name)
  d_long <- merged_long[merged_long$hemisphere == hemi & merged_long$subfield == subfield_name & merged_long$grid_id %in% grid_ids, ]
  if (!length(unique(d_long$grid_id))) return(NULL)

  change_wide <- reshape(d_long[, c("subject_id","grid_id","change")],
                          idvar="subject_id", timevar="grid_id", direction="wide")
  names(change_wide) <- sub("^change\\.", "", names(change_wide))
  grid_cols <- setdiff(names(change_wide), "subject_id")

  if (length(grid_cols) < opt$`min-grid-points`) {
    msg("[%s] only %d grid points (< --min-grid-points=%d), skipping\n", label, length(grid_cols), opt$`min-grid-points`)
    return(NULL)
  }

  dat <- merge(change_wide, participants[, c("subject_id","intervention","age_z","sex")], by="subject_id")
  dat <- dat[!is.na(dat$intervention), , drop=FALSE]
  dat <- dat[complete.cases(dat[, c(grid_cols, "intervention","age_z","sex")]), , drop=FALSE]
  n <- nrow(dat)
  msg("[%s] %d grid points, n=%d complete subjects\n", label, length(grid_cols), n)
  if (n < length(grid_cols) + 10) {
    warning(sprintf("[%s] only %d complete subjects for %d grid-point features -- results should be treated cautiously", label, n, length(grid_cols)))
  }
  if (n < 8) {
    msg("[%s] fewer than 8 complete subjects, skipping\n", label)
    return(NULL)
  }

  Y <- as.matrix(dat[, grid_cols])
  write.csv(cbind(dat[, c("subject_id","intervention")], Y),
            file.path(opt$outdir, paste0(label, "_change_score_matrix.csv")), row.names=FALSE)

  # 1. MANOVA
  manova_fit <- manova(Y ~ intervention + age_z + sex, data=dat)
  manova_summary <- as.data.frame(summary(manova_fit, test="Pillai")$stats)
  manova_p <- manova_summary["intervention", "Pr(>F)"]
  manova_pillai <- manova_summary["intervention", "Pillai"]

  # 2. PCA on covariate-residualized change scores
  resid_mat <- sapply(grid_cols, function(col) residuals(lm(as.formula(paste0("`", col, "` ~ age_z + sex")), data=dat)))
  colnames(resid_mat) <- grid_cols
  pca_fit <- prcomp(resid_mat, center=TRUE, scale.=TRUE)
  var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  n_keep <- max(1, sum(pca_fit$sdev^2 > 1))

  loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop=FALSE])
  loadings_df$grid_id <- rownames(loadings_df)
  loadings_df <- merge(loadings_df, grid_map[, c("grid_id","x","y")], by="grid_id")
  write.csv(loadings_df, file.path(opt$outdir, paste0(label, "_pca_loadings.csv")), row.names=FALSE)

  pc_scores <- as.data.frame(pca_fit$x[, seq_len(n_keep), drop=FALSE])
  pc_scores$intervention <- dat$intervention
  pc_p <- sapply(seq_len(n_keep), function(i) {
    fit <- lm(pc_scores[[i]] ~ intervention, data=pc_scores)
    tt <- summary(fit)$coefficients
    if ("interventionintervention" %in% rownames(tt)) tt["interventionintervention", "Pr(>|t|)"] else NA_real_
  })
  pca_min_p <- if (length(pc_p)) min(pc_p, na.rm=TRUE) else NA_real_
  write.csv(data.frame(component=paste0("PC", seq_len(n_keep)), variance_explained=var_explained[seq_len(n_keep)], intervention_p=pc_p),
            file.path(opt$outdir, paste0(label, "_pca_component_tests.csv")), row.names=FALSE)

  # 3. PLS-DA: leave-one-out cross-validated classification
  y_bin <- as.numeric(dat$intervention) - 1
  X <- scale(Y)
  ncomp <- min(3, length(grid_cols) - 1, n - 2)
  plsda_fit <- tryCatch(plsr(y_bin ~ X, ncomp=ncomp, validation="LOO"), error=function(e) NULL)
  loo_auc <- NA_real_
  loo_acc <- NA_real_
  best_ncomp <- NA_integer_
  if (!is.null(plsda_fit)) {
    press <- as.numeric(plsda_fit$validation$PRESS)
    best_ncomp <- which.min(press)
    loo_pred <- plsda_fit$validation$pred[, 1, best_ncomp]
    loo_class <- as.numeric(loo_pred > 0.5)
    loo_acc <- mean(loo_class == y_bin)
    loo_auc <- tryCatch({
      n_pos <- sum(y_bin == 1); n_neg <- sum(y_bin == 0)
      r <- rank(loo_pred)
      (sum(r[y_bin == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    }, error=function(e) NA_real_)
  }

  list(
    group = label,
    hemisphere = hemi,
    subfield = subfield_name,
    n_grid_points = length(grid_cols),
    n_subjects = n,
    manova_pillai = manova_pillai,
    manova_p = manova_p,
    n_pca_components_kaiser = n_keep,
    pca_min_component_p = pca_min_p,
    plsda_best_ncomp = best_ncomp,
    plsda_loo_accuracy = loo_acc,
    plsda_loo_auc = loo_auc,
    plsda_chance_accuracy = max(mean(y_bin==1), mean(y_bin==0))
  )
}

groups <- unique(grid_map[, c("hemisphere","subfield")])
results <- list()
for (i in seq_len(nrow(groups))) {
  hemi <- groups$hemisphere[i]; subfield_name <- groups$subfield[i]
  grid_ids <- grid_map$grid_id[grid_map$hemisphere == hemi & grid_map$subfield == subfield_name]
  r <- run_group(hemi, subfield_name, grid_ids)
  if (!is.null(r)) results[[r$group]] <- r
}

summary_df <- do.call(rbind, lapply(results, as.data.frame))
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$manova_p_fdr <- p.adjust(summary_df$manova_p, method="fdr")
  summary_df$significant_manova_uncorrected <- summary_df$manova_p < opt$alpha
  summary_df$significant_manova_fdr <- summary_df$manova_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$manova_p), ]
  write.csv(summary_df, file.path(opt$outdir, "cluster_multivariate_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (MANOVA FDR-corrected across %d subfield/hemisphere groups) to cluster_multivariate_summary.csv\n", nrow(summary_df))
  msg("\n%-18s %10s %6s %10s %10s %10s\n", "group", "n_points", "n_sub", "manova_p", "p_fdr", "plsda_acc")
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    msg("%-18s %10d %6d %10.4f %10.4f %10.3f\n", row$group, row$n_grid_points, row$n_subjects, row$manova_p, row$manova_p_fdr, row$plsda_loo_accuracy)
  }
} else {
  warning("No subfield/hemisphere groups were successfully analyzed; no summary written")
}

msg("\nDone. Per-group change-score matrices, PCA loadings (with x/y grid coordinates for mapping back onto the surface), and component tests are in %s/\n", opt$outdir)
