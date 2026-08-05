#!/usr/bin/env Rscript
#
# UNSUPERVISED structure discovery in longitudinal (baseline -> final session)
# change, deliberately WITHOUT using the group label to define the analysis --
# the confirmatory battery (01/05/08/10/24/25/27/28) all condition on
# intervention-vs-control from the start; this script asks a different
# question: are there natural, data-driven modes/subgroups of brain change in
# this cohort at all, and if so, do they happen to align with anything we
# already measured (group, age, sex, baseline moderators)? Generic over any
# domain's tidy table (hippo subfields, amygdala/thalamic nuclei, brainstem
# structures, or Destrieux cortical thickness) via --roi-col/--value-col.
#
# Method:
#   1. Build the subject x ROI change-score matrix (final - baseline, or
#      log-ratio for volume domains), z-scored per ROI.
#   2. PCA for visualization/variance structure (fully unsupervised).
#   3. Hierarchical (Ward.D2) + k-means clustering on the standardized
#      matrix; k chosen by average silhouette width over k=2..6 (data-driven,
#      not assumed).
#   4. Post-hoc profiling of the resulting clusters against group, age, sex,
#      and (optionally) baseline psychological moderators.
#
# IMPORTANT: step 4 is EXPLORATORY / HYPOTHESIS-GENERATING ONLY. Clusters are
# discovered blind to these variables, but testing several of them against
# several emergent clusters after the fact is still a multiple-comparisons
# fishing expedition -- p-values here are reported uncorrected and are NOT a
# confirmatory claim. Treat any hit as something to pre-register and test in
# an independent sample, not as a publishable effect on its own.

suppressPackageStartupMessages({
  library(optparse)
  library(cluster)   # silhouette
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to a domain tidy TSV"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-m", "--moderators"), type="character", default=NULL,
              help="Optional path to a subject-level baseline moderators TSV"),
  make_option(c("--roi-col"), type="character", default="subfield",
              help="Column name identifying the ROI (subfield/nucleus/structure/region) [default %default]"),
  make_option(c("--value-col"), type="character", default="volume",
              help="Column name holding the measurement (volume/value) [default %default]"),
  make_option(c("--roi-set"), type="character", default="all",
              help="Comma-separated ROI names, or 'all' [default %default]"),
  make_option(c("--hemisphere-agg"), type="character", default="sum",
              help="'sum' (volume domains) or 'mean' (thickness) when pooling lh/rh [default %default]"),
  make_option(c("--log-change"), action="store_true", default=FALSE,
              help="Use log(final)-log(baseline) instead of a raw difference (volume domains) [default FALSE]"),
  make_option(c("-o", "--outdir"), type="character", default="results/29_unsupervised_change_clustering",
              help="Output directory [default %default]"),
  make_option(c("--label"), type="character", default="domain", help="Label used in output filenames/summary [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--final-session"), type="character", default=NULL),
  make_option(c("--seed"), type="integer", default=129, help="RNG seed [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }
set.seed(opt$seed)

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-col`
value_col <- opt$`value-col`

# ------------------------------------------------------------------------
# Build subject x ROI change-score matrix
# ------------------------------------------------------------------------
msg("[%s] Loading tidy data from %s...\n", opt$label, opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-col '%s' not found in tidy file (columns: %s)", roi_col, paste(names(tidy), collapse=", ")))
if (!value_col %in% names(tidy)) stop(sprintf("--value-col '%s' not found in tidy file (columns: %s)", value_col, paste(names(tidy), collapse=", ")))
tidy[[value_col]] <- suppressWarnings(as.numeric(tidy[[value_col]]))

roi_set <- if (identical(opt$`roi-set`, "all")) sort(unique(tidy[[roi_col]])) else trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set")

agg_fun <- if (identical(opt$`hemisphere-agg`, "mean")) mean else sum
agg <- aggregate(as.formula(paste0(value_col, " ~ subject_id + session + `", roi_col, "`")), data=roi_dat, FUN=agg_fun)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("[%s] Change score: %s -> %s\n", opt$label, baseline_ses, final_ses)

base_vals <- agg[agg$session == baseline_ses, c("subject_id", roi_col, value_col)]; names(base_vals)[3] <- "baseline"
final_vals <- agg[agg$session == final_ses, c("subject_id", roi_col, value_col)]; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id", roi_col))
merged_long$change <- if (isTRUE(opt$`log-change`)) log(merged_long$final) - log(merged_long$baseline) else merged_long$final - merged_long$baseline

change_wide <- reshape(merged_long[, c("subject_id", roi_col, "change")], idvar="subject_id", timevar=roi_col, direction="wide")
names(change_wide) <- sub("^change\\.", "", names(change_wide))
change_cols <- setdiff(names(change_wide), "subject_id")

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat <- merge(change_wide, participants[, c("subject_id","group","age","sex")], by="subject_id")
if (!is.null(opt$moderators)) {
  moderators <- read.delim(opt$moderators, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  dat <- merge(dat, moderators, by="subject_id", all.x=TRUE)
}

dat <- dat[complete.cases(dat[, change_cols]), , drop=FALSE]
n <- nrow(dat)
msg("[%s] %d ROIs, n=%d complete subjects\n", opt$label, length(change_cols), n)
if (n < length(change_cols) + 10) {
  warning(sprintf("[%s] only %d complete subjects for %d ROI features -- clustering at this ratio should be treated cautiously", opt$label, n, length(change_cols)))
}

X <- scale(as.matrix(dat[, change_cols]))
rownames(X) <- dat$subject_id

# ------------------------------------------------------------------------
# 1. PCA (fully unsupervised, no group label used)
# ------------------------------------------------------------------------
pca_fit <- prcomp(X, center=FALSE, scale.=FALSE)  # X already standardized
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
n_keep <- max(2, min(sum(pca_fit$sdev^2 > 1), length(change_cols)))

loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop=FALSE])
loadings_df$roi <- rownames(loadings_df)
write.csv(loadings_df, file.path(opt$outdir, paste0(opt$label, "_pca_loadings.csv")), row.names=FALSE)
write.csv(data.frame(component=paste0("PC", seq_along(var_explained)), variance_explained=var_explained),
          file.path(opt$outdir, paste0(opt$label, "_pca_variance_explained.csv")), row.names=FALSE)

pc_scores <- as.data.frame(pca_fit$x)
pc_scores$subject_id <- dat$subject_id

# ------------------------------------------------------------------------
# 2. Unsupervised clustering: k chosen by average silhouette width, k=2..6
# ------------------------------------------------------------------------
dist_mat <- dist(X)
hc <- hclust(dist_mat, method="ward.D2")

k_range <- 2:min(6, n - 1)
sil_widths <- sapply(k_range, function(k) {
  cl <- cutree(hc, k=k)
  if (length(unique(cl)) < 2) return(NA_real_)
  mean(silhouette(cl, dist_mat)[, "sil_width"])
})
best_k <- k_range[which.max(sil_widths)]
msg("[%s] Silhouette widths by k: %s\n", opt$label, paste(sprintf("k=%d:%.3f", k_range, sil_widths), collapse="  "))
msg("[%s] Best k (hierarchical, Ward.D2): %d (avg silhouette = %.3f)\n", opt$label, best_k, max(sil_widths, na.rm=TRUE))

hc_clusters <- cutree(hc, k=best_k)
km_fit <- kmeans(X, centers=best_k, nstart=50)
km_clusters <- km_fit$cluster

# hclust/kmeans cluster labels are each arbitrary integers (1..k), not
# aligned to each other -- "agreement" must be evaluated over the best
# label permutation, not the raw diagonal (which would understate agreement
# by an arbitrary amount whenever the two methods happen to number an
# identical split in opposite order).
combinat_permutations <- function(x) {
  if (length(x) == 1) return(list(x))
  do.call(c, lapply(seq_along(x), function(i) {
    lapply(combinat_permutations(x[-i]), function(rest) c(x[i], rest))
  }))
}
# hclust/kmeans cluster labels are each arbitrary integers (1..k), not
# aligned to each other -- "agreement" must be evaluated over the best
# label permutation, not the raw diagonal (which would understate agreement
# by an arbitrary amount whenever the two methods happen to number an
# identical split in opposite order).
agreement <- table(hclust=hc_clusters, kmeans=km_clusters)
best_label_match_agreement <- if (best_k == 2) {
  max(sum(diag(agreement)), sum(agreement) - sum(diag(agreement)))
} else {
  perms <- combinat_permutations(seq_len(best_k))
  max(sapply(perms, function(p) sum(agreement[cbind(seq_len(best_k), p)])))
}

cluster_df <- data.frame(subject_id=dat$subject_id, hclust_cluster=hc_clusters, kmeans_cluster=km_clusters)
cluster_df <- merge(cluster_df, pc_scores[, c("subject_id","PC1","PC2")], by="subject_id")
write.csv(cluster_df, file.path(opt$outdir, paste0(opt$label, "_cluster_assignments.csv")), row.names=FALSE)
saveRDS(hc, file.path(opt$outdir, paste0(opt$label, "_hclust.rds")))

write.csv(data.frame(k=k_range, avg_silhouette=sil_widths),
          file.path(opt$outdir, paste0(opt$label, "_silhouette_by_k.csv")), row.names=FALSE)

# ------------------------------------------------------------------------
# 3. EXPLORATORY post-hoc profiling of the emergent (hclust) clusters
# ------------------------------------------------------------------------
dat$hclust_cluster <- factor(hc_clusters)
profile_rows <- list()

# vs. group (categorical)
tab_group <- table(dat$hclust_cluster, dat$group)
p_group <- tryCatch(fisher.test(tab_group, simulate.p.value=TRUE, B=10000)$p.value, error=function(e) NA_real_)
profile_rows[["group"]] <- data.frame(variable="group", test="fisher.test (simulated p)", p_value=p_group, stringsAsFactors=FALSE)

# vs. sex (categorical)
tab_sex <- table(dat$hclust_cluster, dat$sex)
p_sex <- tryCatch(fisher.test(tab_sex, simulate.p.value=TRUE, B=10000)$p.value, error=function(e) NA_real_)
profile_rows[["sex"]] <- data.frame(variable="sex", test="fisher.test (simulated p)", p_value=p_sex, stringsAsFactors=FALSE)

# vs. age (continuous)
p_age <- tryCatch(kruskal.test(age ~ hclust_cluster, data=dat)$p.value, error=function(e) NA_real_)
profile_rows[["age"]] <- data.frame(variable="age", test="kruskal.test", p_value=p_age, stringsAsFactors=FALSE)

# vs. baseline moderators, if supplied
if (!is.null(opt$moderators)) {
  mod_cols <- setdiff(names(moderators), "subject_id")
  for (col in mod_cols) {
    if (!col %in% names(dat)) next
    p_mod <- tryCatch(kruskal.test(as.formula(paste0("`", col, "` ~ hclust_cluster")), data=dat)$p.value, error=function(e) NA_real_)
    profile_rows[[col]] <- data.frame(variable=col, test="kruskal.test", p_value=p_mod, stringsAsFactors=FALSE)
  }
}

profile_df <- do.call(rbind, profile_rows)
profile_df$note <- "EXPLORATORY -- uncorrected, hypothesis-generating only (clusters discovered blind to these variables, but many variables tested post-hoc)"
write.csv(profile_df, file.path(opt$outdir, paste0(opt$label, "_cluster_profile_exploratory.csv")), row.names=FALSE)

cluster_sizes <- table(dat$hclust_cluster)

# ------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------
con <- file(file.path(opt$outdir, paste0(opt$label, "_summary.txt")), open="wt")
cat(sprintf(
  "UNSUPERVISED change-pattern discovery: %s (%s -> %s), %d ROIs, n=%d\n\n1. PCA: %d components retained (Kaiser), top-2 explain %.1f%% of variance\n\n2. Clustering (Ward.D2 hierarchical, k chosen by avg silhouette over k=2..6):\n   Best k = %d (avg silhouette = %.3f)\n   Cluster sizes: %s\n   Agreement with k-means (same k): %d/%d subjects assigned to the same-labeled cluster by both methods\n\n3. EXPLORATORY post-hoc profiling of the %d emergent clusters (uncorrected,\n   hypothesis-generating only -- see the note in the CSV):\n",
  opt$label, baseline_ses, final_ses, length(change_cols), n,
  n_keep, 100*sum(var_explained[1:min(2,length(var_explained))]),
  best_k, max(sil_widths, na.rm=TRUE), paste(sprintf("cluster %s: n=%d", names(cluster_sizes), cluster_sizes), collapse=", "),
  best_label_match_agreement, n, best_k
), file=con)
for (i in seq_len(nrow(profile_df))) {
  row <- profile_df[i, ]
  cat(sprintf("   - %s (%s): p = %.4f%s\n", row$variable, row$test, row$p_value, if (!is.na(row$p_value) && row$p_value < 0.05) " (nominal)" else ""), file=con)
}
cat("\nSee _pca_loadings.csv for which ROIs drive each component, and\n_cluster_assignments.csv for per-subject cluster membership + PC1/PC2 scores\n(useful for a scatter plot colored by cluster and, separately, by group).\n", file=con)
close(con)

msg("[%s] Done. See %s_summary.txt\n", opt$label, opt$label)
