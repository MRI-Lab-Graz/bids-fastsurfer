#!/usr/bin/env Rscript
#
# Unsupervised, data-driven trajectory clustering: rather than testing a
# pre-specified group (dance/running vs control) or a post-hoc "responder"
# split defined from the same outcome being tested (both circular traps
# this project has deliberately avoided), this asks a genuinely different
# question -- does the data itself contain distinct SHAPES of change over
# time, discovered without using the group label at all?
#
# Uses mclust (Gaussian finite mixture modelling with BIC-based model
# selection) specifically because it is allowed to conclude "1 component"
# -- i.e. NO real subgroups, just continuous individual variation around a
# single mean -- rather than a method like k-means that always forces a
# split regardless of whether real cluster structure exists. Reporting "no
# clusters found" honestly is as valid an outcome as finding some.
#
# Clustering features capture TRAJECTORY SHAPE, not just endpoint change:
# for each subject, early change (baseline -> mid) and late change (mid ->
# final) on the two headline summary measures (Whole_hippocampus,
# Whole_amygdala) -- a 4-dimensional space distinguishing e.g. "grew early
# then plateaued" from "grew steadily" from "declined throughout", across
# both regions at once.
#
# The group label is used ONLY after clustering, to test whether cluster
# membership associates with dance/running vs control -- discovered
# structure tested against group, never the reverse.

suppressPackageStartupMessages({
  library(optparse)
  library(mclust)
})

option_list <- list(
  make_option(c("--hippo-tidy"), type="character", help="Path to combined hippocampus tidy TSV"),
  make_option(c("--amygdala-tidy"), type="character", help="Path to combined amygdala tidy TSV"),
  make_option(c("-p", "--participants"), type="character", help="Path to combined participants TSV"),
  make_option(c("-o", "--outdir"), type="character", default="results/17_unsupervised_trajectory_clusters"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--mid-session"), type="character", default=NULL),
  make_option(c("--final-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$`hippo-tidy`) || is.null(opt$`amygdala-tidy`) || is.null(opt$participants)) {
  stop("--hippo-tidy, --amygdala-tidy, and --participants are all required")
}
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ------------------------------------------------------------------------
# Build the 4D feature matrix: subject x {hippo_early, hippo_late,
# amyg_early, amyg_late}
# ------------------------------------------------------------------------
summary_change <- function(tidy_path, roi_name, roi_col, roi_val) {
  tidy <- read.delim(tidy_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
  d <- tidy[tidy[[roi_col]] == roi_val, , drop=FALSE]
  agg <- aggregate(volume ~ subject_id + session, data=d, FUN=sum)
  sessions <- sort(unique(agg$session))
  baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
  final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
  mid_ses <- if (!is.null(opt$`mid-session`)) opt$`mid-session` else setdiff(sessions, c(baseline_ses, final_ses))[1]

  w <- reshape(agg[agg$session %in% c(baseline_ses, mid_ses, final_ses), ], idvar="subject_id", timevar="session", direction="wide")
  names(w) <- sub("^volume\\.", "", names(w))
  needed <- c("subject_id", baseline_ses, mid_ses, final_ses)
  if (!all(needed %in% names(w))) stop(sprintf("Missing expected sessions for %s: have %s", roi_name, paste(names(w), collapse=", ")))
  w <- w[complete.cases(w[, needed]), needed]
  names(w) <- c("subject_id", "baseline", "mid", "final")
  w[[paste0(roi_name, "_early")]] <- log(w$mid) - log(w$baseline)
  w[[paste0(roi_name, "_late")]] <- log(w$final) - log(w$mid)
  w[, c("subject_id", paste0(roi_name, "_early"), paste0(roi_name, "_late"))]
}

msg("Building trajectory features (early/late change) for Whole_hippocampus and Whole_amygdala...\n")
hippo_feat <- summary_change(opt$`hippo-tidy`, "hippo", "subfield", "Whole_hippocampus")
amyg_feat <- summary_change(opt$`amygdala-tidy`, "amyg", "nucleus", "Whole_amygdala")
feat <- merge(hippo_feat, amyg_feat, by="subject_id")

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
feat <- merge(feat, participants[, c("subject_id","study","group")], by="subject_id")
feat$dance <- factor(ifelse(feat$group == "intervention", "intervention", "control"), levels=c("control","intervention"))

feature_cols <- c("hippo_early","hippo_late","amyg_early","amyg_late")
X <- scale(as.matrix(feat[, feature_cols]))
msg("Subjects with complete 4D trajectory features: %d\n", nrow(X))

write.csv(cbind(subject_id=feat$subject_id, feat[, feature_cols]), file.path(opt$outdir, "trajectory_features.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# Gaussian mixture model with BIC-based selection of the number of
# components (mclust tries G=1..9 and multiple covariance shapes; G=1
# means "no real clusters" and is a legitimate, honest outcome)
# ------------------------------------------------------------------------
msg("Fitting Gaussian mixture models (G=1..9), selecting by BIC...\n")
set.seed(134)
mc <- Mclust(X, G=1:9)
msg("Selected model: G=%d components, covariance model '%s', BIC=%.1f\n", mc$G, mc$modelName, mc$bic)

tryCatch(
  write.csv(as.data.frame(unclass(mc$BIC)), file.path(opt$outdir, "bic_by_ncomponents.csv")),
  error = function(e) msg("(Skipping BIC-by-G export: %s)\n", conditionMessage(e))
)

con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf("Unsupervised trajectory clustering (n=%d, 4D: hippo/amygdala x early/late change)\n\n", nrow(X)), file=con)
cat(sprintf("BIC-selected model: G=%d component(s), covariance shape '%s'\n\n", mc$G, mc$modelName), file=con)

if (mc$G == 1) {
  cat("RESULT: BIC selects a single component -- the data-driven search finds NO evidence\n", file=con)
  cat("of distinct trajectory subgroups. Individual variation in the timing/shape of\n", file=con)
  cat("hippocampal and amygdala change looks continuous, not clustered.\n", file=con)
  msg("\nNo real clusters found (G=1 selected by BIC). See %s/summary.txt\n", opt$outdir)
} else {
  feat$cluster <- factor(mc$classification)
  write.csv(feat, file.path(opt$outdir, "cluster_assignments.csv"), row.names=FALSE)

  # Describe clusters: size + mean trajectory shape per cluster
  cluster_means <- aggregate(feat[, feature_cols], by=list(cluster=feat$cluster), FUN=mean)
  cluster_sizes <- table(feat$cluster)
  cat("Cluster sizes:", paste(names(cluster_sizes), cluster_sizes, sep="=", collapse=", "), "\n\n", file=con)
  cat("Cluster mean trajectories (log-volume change per interval):\n", file=con)
  write.csv(cluster_means, file.path(opt$outdir, "cluster_mean_trajectories.csv"), row.names=FALSE)
  for (i in seq_len(nrow(cluster_means))) {
    row <- cluster_means[i, ]
    cat(sprintf("  Cluster %s: hippo_early=%.4f, hippo_late=%.4f, amyg_early=%.4f, amyg_late=%.4f\n",
                row$cluster, row$hippo_early, row$hippo_late, row$amyg_early, row$amyg_late), file=con)
  }

  # Test cluster x group association -- discovered structure tested
  # against group, never the reverse
  ctab <- table(feat$cluster, feat$dance)
  cat("\nCluster x intervention-group contingency table:\n", file=con)
  capture.output(print(ctab), file=con, append=TRUE)
  fisher_test <- fisher.test(ctab, simulate.p.value = (nrow(ctab) > 2 || ncol(ctab) > 2))
  cat(sprintf("\nFisher's exact test (cluster x group independence): p = %.4f\n", fisher_test$p.value), file=con)

  # LOSO for the cluster x group association specifically
  msg("\nRunning LOSO for the cluster x group association test...\n")
  loso_p <- sapply(seq_len(nrow(feat)), function(i) {
    ctab_i <- table(feat$cluster[-i], feat$dance[-i])
    tryCatch(fisher.test(ctab_i, simulate.p.value = (nrow(ctab_i) > 2 || ncol(ctab_i) > 2))$p.value, error=function(e) NA_real_)
  })
  loso_p <- loso_p[!is.na(loso_p)]
  n_flips <- sum((loso_p < opt$alpha) != (fisher_test$p.value < opt$alpha))
  cat(sprintf("LOSO: range [%.4f, %.4f], %d/%d exclusions flip significance\n",
              min(loso_p), max(loso_p), n_flips, length(loso_p)), file=con)
  write.csv(data.frame(loso_p=loso_p), file.path(opt$outdir, "loso_cluster_group_association.csv"), row.names=FALSE)

  msg("Clusters found: G=%d. See %s/summary.txt for full description and LOSO validation.\n", mc$G, opt$outdir)
}
