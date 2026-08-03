#!/usr/bin/env Rscript
#
# Multivariate pattern analysis on cortical thickness, run SEPARATELY within
# each of 9 a priori functional clusters of Destrieux (aparc.a2009s) regions
# (configs/destrieux_functional_clusters.tsv: Sensorimotor, Visual,
# Auditory_Language, Temporal_Association, Prefrontal, Orbitofrontal,
# Cingulate_Limbic, Insula, Parietal_Association) rather than across all 74
# regions jointly -- the question is not "is there a pattern somewhere in
# cortex" (74-region multivariate model, no single interpretable pattern)
# but "does a functionally coherent set of regions show a coordinated
# intervention-related change", one cluster at a time.
#
# Same three complementary methods as 06_manova_pca_pattern.R, applied per
# cluster:
#   1. MANOVA (Pillai's trace) on the subject x region change-score matrix
#      (baseline -> final session), covarying age_z (rank-transformed) + sex.
#   2. PCA on the covariate-residualized change-score matrix -- latent
#      components of coordinated change, tested against intervention.
#   3. PLS-DA (1-component PLSR of the binary intervention label on the
#      change-score matrix, leave-one-out cross-validated) -- a model-free
#      check of whether the region-set's joint pattern discriminates
#      intervention from control out-of-sample, not just in-sample.
#
# lh/rh are averaged (not summed) per region -- thickness is not additive
# across hemisphere the way volume is (see extract_freesurfer.py's docstring
# on why eTIV/log-transform don't apply to thickness either).
#
# MANOVA p-values are FDR-corrected across the 9 clusters (one test per
# cluster, analogous to the per-ROI FDR correction elsewhere in this project).

suppressPackageStartupMessages({
  library(optparse)
  library(pls)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to aparc_a2009s_tidy_2tp_fs82.tsv (from extract_freesurfer.py --source aparc.a2009s)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-c", "--clusters"), type="character",
              default="configs/destrieux_functional_clusters.tsv",
              help="Path to region -> cluster mapping TSV [default %default]"),
  make_option(c("-o", "--outdir"), type="character", default="results/28_cortical_cluster_multivariate",
              help="Output directory [default %default]"),
  make_option(c("--intervention-groups"), type="character", default="Single_2wk,Single_4wk,Group_2wk,Group_4wk",
              help="Comma-separated group values pooled into the 'intervention' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="Control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value to use as baseline [default: earliest session present]"),
  make_option(c("--final-session"), type="character", default=NULL,
              help="Session value to use as endpoint [default: latest session present]"),
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
if (!file.exists(opt$clusters)) stop(sprintf("clusters file not found: %s", opt$clusters))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------
msg("Loading tidy cortical thickness data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$value <- suppressWarnings(as.numeric(tidy$value))

clusters <- read.delim(opt$clusters, header=TRUE, sep="\t", stringsAsFactors=FALSE)
missing_regions <- setdiff(unique(tidy$region), clusters$region)
if (length(missing_regions)) {
  warning(sprintf("Regions in tidy data but not in --clusters mapping (excluded): %s", paste(missing_regions, collapse=", ")))
}

sessions <- sort(unique(tidy$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Change score: %s -> %s\n", baseline_ses, final_ses)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))
participants$intervention <- factor(ifelse(participants$group %in% intervention_groups, "intervention",
                                     ifelse(participants$group == control_group, "control", NA)), levels=c("control","intervention"))
participants$sex <- factor(participants$sex)
participants$age_z <- as.numeric(scale(rank(participants$age)))

# region-mean thickness (across hemisphere) per subject/session/region
agg <- aggregate(value ~ subject_id + session + region, data=tidy, FUN=mean)

base_vals <- agg[agg$session == baseline_ses, c("subject_id","region","value")]; names(base_vals)[3] <- "baseline"
final_vals <- agg[agg$session == final_ses, c("subject_id","region","value")]; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id","region"))
merged_long$change <- merged_long$final - merged_long$baseline

change_wide_all <- reshape(merged_long[, c("subject_id","region","change")],
                            idvar="subject_id", timevar="region", direction="wide")
names(change_wide_all) <- sub("^change\\.", "", names(change_wide_all))

change_wide_all <- merge(change_wide_all, participants[, c("subject_id","intervention","age_z","sex")], by="subject_id")
change_wide_all <- change_wide_all[!is.na(change_wide_all$intervention), , drop=FALSE]

# ------------------------------------------------------------------------
# Per-cluster MANOVA + PCA + PLS-DA
# ------------------------------------------------------------------------
run_cluster <- function(cluster_name, region_cols) {
  region_cols <- intersect(region_cols, names(change_wide_all))
  if (length(region_cols) < 2) {
    msg("[%s] fewer than 2 regions available, skipping\n", cluster_name)
    return(NULL)
  }
  dat <- change_wide_all[complete.cases(change_wide_all[, c(region_cols, "intervention","age_z","sex")]), , drop=FALSE]
  n <- nrow(dat)
  msg("[%s] %d regions, n=%d complete subjects\n", cluster_name, length(region_cols), n)
  if (n < length(region_cols) + 10) {
    warning(sprintf("[%s] only %d complete subjects for %d region features -- results should be treated cautiously", cluster_name, n, length(region_cols)))
  }

  Y <- as.matrix(dat[, region_cols])
  write.csv(cbind(dat[, c("subject_id","intervention")], Y),
            file.path(opt$outdir, paste0(cluster_name, "_change_score_matrix.csv")), row.names=FALSE)

  # 1. MANOVA
  manova_fit <- manova(Y ~ intervention + age_z + sex, data=dat)
  manova_summary <- as.data.frame(summary(manova_fit, test="Pillai")$stats)
  manova_p <- manova_summary["intervention", "Pr(>F)"]
  manova_pillai <- manova_summary["intervention", "Pillai"]

  # 2. PCA on covariate-residualized change scores
  resid_mat <- sapply(region_cols, function(col) residuals(lm(as.formula(paste0("`", col, "` ~ age_z + sex")), data=dat)))
  colnames(resid_mat) <- region_cols
  pca_fit <- prcomp(resid_mat, center=TRUE, scale.=TRUE)
  var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  n_keep <- max(1, sum(pca_fit$sdev^2 > 1))

  loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop=FALSE])
  loadings_df$region <- rownames(loadings_df)
  write.csv(loadings_df, file.path(opt$outdir, paste0(cluster_name, "_pca_loadings.csv")), row.names=FALSE)

  pc_scores <- as.data.frame(pca_fit$x[, seq_len(n_keep), drop=FALSE])
  pc_scores$intervention <- dat$intervention
  pc_p <- sapply(seq_len(n_keep), function(i) {
    fit <- lm(pc_scores[[i]] ~ intervention, data=pc_scores)
    tt <- summary(fit)$coefficients
    if ("interventionintervention" %in% rownames(tt)) tt["interventionintervention", "Pr(>|t|)"] else NA_real_
  })
  pca_min_p <- if (length(pc_p)) min(pc_p, na.rm=TRUE) else NA_real_
  write.csv(data.frame(component=paste0("PC", seq_len(n_keep)), variance_explained=var_explained[seq_len(n_keep)], intervention_p=pc_p),
            file.path(opt$outdir, paste0(cluster_name, "_pca_component_tests.csv")), row.names=FALSE)

  # 3. PLS-DA: leave-one-out cross-validated classification of intervention
  # vs control from the region-change pattern (out-of-sample check, not just
  # in-sample fit -- a coherent multivariate pattern should generalize).
  y_bin <- as.numeric(dat$intervention) - 1
  X <- scale(Y)
  ncomp <- min(3, length(region_cols) - 1, n - 2)
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
    # Rank-sum (Mann-Whitney U) AUC: robust closed-form estimator, no manual
    # ROC integration needed. AUC = P(score of a random positive > score of
    # a random negative), estimated from the LOO-predicted scores' ranks.
    loo_auc <- tryCatch({
      n_pos <- sum(y_bin == 1); n_neg <- sum(y_bin == 0)
      r <- rank(loo_pred)
      (sum(r[y_bin == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    }, error=function(e) NA_real_)
  }

  list(
    cluster = cluster_name,
    n_regions = length(region_cols),
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

cluster_names <- sort(unique(clusters$cluster))
results <- list()
for (cl in cluster_names) {
  region_cols <- clusters$region[clusters$cluster == cl]
  r <- run_cluster(cl, region_cols)
  if (!is.null(r)) results[[cl]] <- r
}

summary_df <- do.call(rbind, lapply(results, as.data.frame))
if (!is.null(summary_df) && nrow(summary_df)) {
  summary_df$manova_p_fdr <- p.adjust(summary_df$manova_p, method="fdr")
  summary_df$significant_manova_uncorrected <- summary_df$manova_p < opt$alpha
  summary_df$significant_manova_fdr <- summary_df$manova_p_fdr < opt$alpha
  summary_df <- summary_df[order(summary_df$manova_p), ]
  write.csv(summary_df, file.path(opt$outdir, "cluster_multivariate_summary.csv"), row.names=FALSE)
  msg("\nWrote combined summary (MANOVA FDR-corrected across %d clusters) to cluster_multivariate_summary.csv\n", nrow(summary_df))
  msg("\n%-22s %8s %6s %10s %10s %10s\n", "cluster", "n_regio", "n_sub", "manova_p", "p_fdr", "plsda_acc")
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    msg("%-22s %8d %6d %10.4f %10.4f %10.3f\n", row$cluster, row$n_regions, row$n_subjects, row$manova_p, row$manova_p_fdr, row$plsda_loo_accuracy)
  }
  con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
  cat(sprintf(
    "Cortical thickness cluster multivariate pattern analysis: %s -> %s change scores\n9 a priori functional clusters (configs/destrieux_functional_clusters.tsv), tested separately.\n\nFor each cluster: MANOVA (Pillai's trace, joint test of intervention vs control across the\ncluster's regions), PCA (latent components of coordinated change, Kaiser criterion),\nand PLS-DA (leave-one-out cross-validated classification of intervention vs control).\n\n", baseline_ses, final_ses), file=con)
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    cat(sprintf(
      "%s (%d regions, n=%d)\n  MANOVA: Pillai=%.4f, p=%.4f, p_fdr=%.4f%s\n  PCA: %d components (Kaiser), smallest per-component intervention p=%.4f\n  PLS-DA: LOO accuracy=%.3f (chance=%.3f), LOO AUC=%.3f\n\n",
      row$cluster, row$n_regions, row$n_subjects, row$manova_pillai, row$manova_p, row$manova_p_fdr,
      if (row$significant_manova_fdr) " *" else "",
      row$n_pca_components_kaiser, row$pca_min_component_p,
      row$plsda_loo_accuracy, row$plsda_chance_accuracy, row$plsda_loo_auc
    ), file=con)
  }
  close(con)
} else {
  warning("No clusters were successfully analyzed; no summary written")
}

msg("\nDone. Per-cluster change-score matrices, PCA loadings, and component tests are in %s/\n", opt$outdir)
