#!/usr/bin/env Rscript
#
# Vertex-wise (grid-point-wise) mass-univariate analysis of hippocampal
# subfield THICKNESS change, with cluster-based permutation inference. This
# is the point where hipsta's output is used for something genuinely
# different from the volume/cortical-thickness pipelines elsewhere in this
# repo: those test pre-defined ROIs; this tests every one of the ~860 grid
# points independently and then asks whether any spatially CONTIGUOUS set of
# points shows a coordinated effect -- which can catch a focal effect (e.g.
# thinning confined to the head of CA1) that a whole-subfield mean
# (29/30) or even a whole-subfield multivariate pattern (31) would still
# partly average away.
#
# This only works because hipsta's cube parametrization puts grid index
# (x, y) at the same relative anatomical position for every subject with NO
# separate registration step (see build_hipsta_grid_subfield_map.py's
# docstring) -- the grid's own (x, y) index IS the correspondence, and also
# gives a trivial spatial adjacency graph (up/down/left/right neighbours) to
# form clusters on, with no cortical-surface mesh or atlas required.
#
# Method (analogous to FSL PALM / SPM's cluster-based permutation inference,
# implemented here directly in base R rather than via an external
# neuroimaging stats package, since the "image" is just a small 41x21 table):
#
#   1. Reduce each grid point to one baseline->final CHANGE score per subject
#      (as in 28/31), and fit, independently at every grid point, the SAME
#      linear model:
#         change ~ intervention + age_z + sex
#      Fit as one vectorized OLS across all grid points at once (shared
#      design matrix X, one column of Y per grid point) rather than 861
#      separate lm() calls -- algebraically identical, just fast.
#   2. Threshold the resulting t-map at --cluster-forming-p (two-sided) and
#      find spatially connected clusters (separately for the positive and
#      negative tail) via breadth-first search over 4-connected grid
#      neighbours. Cluster MASS = sum(|t| - threshold) over the cluster's
#      points (more sensitive than cluster size to a few strong points).
#   3. Permute the intervention label across subjects --n-permutations times
#      (age_z, sex, and each grid point's Y are held fixed), repeat steps 1-2
#      for each permutation, and record the single largest cluster mass
#      (either polarity) per permutation as the null distribution.
#   4. Each OBSERVED cluster's corrected p-value is the fraction of
#      permutations whose max cluster mass equalled or exceeded it -- this
#      controls the family-wise error rate across the whole grid without
#      needing FDR across 861 near-duplicate (spatially correlated) tests.
#
# Caveat (documented, not hidden): step 3 permutes the label of the tested
# covariate directly rather than using Freedman-Lane residualization of
# nuisance covariates before permuting -- simpler, and standard practice for
# a single, cleanly-defined between-subject group contrast, but technically
# an approximation when nuisance covariates (age, sex) are present. For this
# design (group assigned independent of age/sex at randomization) the
# approximation is expected to be good; if age/sex confounding with group
# turns out to be non-trivial in this cohort, upgrading to Freedman-Lane
# permutation would be the next step.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_thickness_grid_tidy.tsv (from extract_hipsta_thickness.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, <group-column>, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/32_hippo_thickness_vertexwise",
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
  make_option(c("--cluster-forming-p"), type="double", default=0.01,
              help="Two-sided uncorrected p-value threshold for cluster formation [default %default]"),
  make_option(c("--n-permutations"), type="integer", default=1000,
              help="Number of label permutations for the cluster-mass null distribution [default %default]"),
  make_option(c("--min-subjects"), type="integer", default=15,
              help="Minimum complete subjects required to run a hemisphere [default %default]"),
  make_option(c("--alpha"), type="double", default=0.05, help="Cluster-level significance threshold [default %default]"),
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

sessions <- sort(unique(tidy$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Change score: %s -> %s\n", baseline_ses, final_ses)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id", group_col, "age", "sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s (--group-column='%s')", paste(missing_pcols, collapse=", "), group_col))
names(participants)[names(participants) == group_col] <- "group"
participants$intervention <- ifelse(participants$group %in% intervention_groups, 1L,
                              ifelse(participants$group == control_group, 0L, NA_integer_))
participants$sex <- factor(participants$sex)
participants$age_z <- as.numeric(scale(rank(participants$age)))

base_vals <- tidy[tidy$session == baseline_ses, c("subject_id","hemisphere","x","y","thickness")]
names(base_vals)[5] <- "baseline"
final_vals <- tidy[tidy$session == final_ses, c("subject_id","hemisphere","x","y","thickness")]
names(final_vals)[5] <- "final"
merged <- merge(base_vals, final_vals, by=c("subject_id","hemisphere","x","y"))
merged$change <- merged$final - merged$baseline

# ------------------------------------------------------------------------
# Vectorized OLS across all grid points at once: Y is n_subjects x n_points,
# X is the (shared) design matrix. beta = (X'X)^-1 X'Y for every point in one
# matrix multiply, instead of one lm() call per point.
# ------------------------------------------------------------------------
fit_tmap <- function(Y, X, intervention_col) {
  XtX_inv <- solve(crossprod(X))
  beta <- XtX_inv %*% t(X) %*% Y
  resid <- Y - X %*% beta
  df <- nrow(X) - ncol(X)
  sigma2 <- colSums(resid^2) / df
  se <- sqrt(sigma2 * XtX_inv[intervention_col, intervention_col])
  beta[intervention_col, ] / se
}

# ------------------------------------------------------------------------
# Grid clustering: BFS over 4-connected (x,y) neighbours among points whose
# |t| exceeds the cluster-forming threshold, split by sign.
# ------------------------------------------------------------------------
find_clusters <- function(x, y, t, threshold) {
  idx_pos <- which(t > threshold)
  idx_neg <- which(t < -threshold)
  coord_key <- function(i) paste(x[i], y[i])
  all_keys <- coord_key(seq_along(x))

  bfs_clusters <- function(idx, polarity) {
    if (!length(idx)) return(list())
    keyset <- coord_key(idx)
    key_to_idx <- setNames(idx, keyset)
    visited <- setNames(rep(FALSE, length(idx)), keyset)
    neighbours_of <- function(k) {
      parts <- as.integer(strsplit(k, " ")[[1]])
      xi <- parts[1]; yi <- parts[2]
      c(paste(xi+1, yi), paste(xi-1, yi), paste(xi, yi+1), paste(xi, yi-1))
    }
    clusters <- list()
    for (k0 in keyset) {
      if (visited[k0]) next
      queue <- k0; visited[k0] <- TRUE
      members <- c()
      while (length(queue)) {
        k <- queue[1]; queue <- queue[-1]
        members <- c(members, k)
        for (nb in neighbours_of(k)) {
          if (!is.na(visited[nb]) && !visited[nb]) {
            visited[nb] <- TRUE
            queue <- c(queue, nb)
          }
        }
      }
      member_idx <- key_to_idx[members]
      clusters[[length(clusters) + 1]] <- list(
        polarity = polarity,
        idx = member_idx,
        size = length(member_idx),
        mass = sum(abs(t[member_idx]) - threshold)
      )
    }
    clusters
  }

  c(bfs_clusters(idx_pos, "positive"), bfs_clusters(idx_neg, "negative"))
}

# ------------------------------------------------------------------------
# Per-hemisphere analysis
# ------------------------------------------------------------------------
overall_cluster_rows <- list()

for (hemi in sort(unique(merged$hemisphere))) {
  msg("\n=== Hemisphere: %s ===\n", hemi)
  d_hemi <- merged[merged$hemisphere == hemi, , drop=FALSE]

  # reshape()'s multi-key timevar naming isn't reliable across R versions;
  # build the subject x grid-point wide matrix explicitly from x/y instead.
  point_coords <- unique(d_hemi[, c("x","y")])
  point_coords <- point_coords[order(point_coords$x, point_coords$y), ]
  point_cols <- paste0("x", point_coords$x, "_y", point_coords$y)

  subj_ids <- sort(unique(d_hemi$subject_id))
  Y_full <- matrix(NA_real_, nrow=length(subj_ids), ncol=nrow(point_coords),
                    dimnames=list(subj_ids, point_cols))
  for (i in seq_len(nrow(d_hemi))) {
    r <- d_hemi[i, ]
    Y_full[as.character(r$subject_id), paste0("x", r$x, "_y", r$y)] <- r$change
  }

  dat <- merge(data.frame(subject_id=subj_ids, stringsAsFactors=FALSE),
               participants[, c("subject_id","intervention","age_z","sex")], by="subject_id")
  dat <- dat[!is.na(dat$intervention), , drop=FALSE]
  # participants$sex was factored over the full participants file (which has
  # a third "n/a" level for a few subjects with no recorded sex) -- if this
  # hemisphere's actual complete-case subset happens to contain none of
  # those subjects, the stale "n/a" level survives as an all-zero column in
  # model.matrix() below, making X'X exactly singular. Drop unused levels
  # against the subset actually being modelled, as 29/30 already do.
  dat$sex <- droplevels(dat$sex)
  Y_full <- Y_full[dat$subject_id, , drop=FALSE]
  complete_points <- colSums(is.na(Y_full)) == 0
  Y <- Y_full[, complete_points, drop=FALSE]
  pts <- point_coords[complete_points, , drop=FALSE]

  n <- nrow(dat)
  msg("%d complete subjects, %d grid points with full coverage\n", n, ncol(Y))
  if (n < opt$`min-subjects`) {
    msg("(%s) fewer than --min-subjects=%d complete subjects, skipping\n", hemi, opt$`min-subjects`)
    next
  }
  if (!ncol(Y)) {
    msg("(%s) no grid points with complete coverage across all subjects, skipping\n", hemi)
    next
  }

  X <- model.matrix(~ intervention + age_z + sex, data=dat)
  intervention_col <- which(colnames(X) == "intervention")
  df_resid <- n - ncol(X)

  t_obs <- fit_tmap(Y, X, intervention_col)
  p_obs <- 2 * pt(-abs(t_obs), df=df_resid)
  t_crit <- qt(1 - opt$`cluster-forming-p` / 2, df=df_resid)

  tmap_df <- data.frame(x=pts$x, y=pts$y, t=as.numeric(t_obs), p_uncorrected=as.numeric(p_obs))
  write.csv(tmap_df, file.path(opt$outdir, paste0(hemi, "_tmap.csv")), row.names=FALSE)

  clusters_obs <- find_clusters(pts$x, pts$y, t_obs, t_crit)
  msg("Cluster-forming threshold |t| > %.3f (uncorrected p < %.3f): %d clusters found\n",
      t_crit, opt$`cluster-forming-p`, length(clusters_obs))

  if (!length(clusters_obs)) {
    msg("(%s) no suprathreshold clusters at the cluster-forming threshold; nothing to test\n", hemi)
    next
  }

  # ---- permutation null: max cluster mass per permutation ----
  msg("Running %d label permutations...\n", opt$`n-permutations`)
  null_max_mass <- numeric(opt$`n-permutations`)
  for (perm in seq_len(opt$`n-permutations`)) {
    X_perm <- X
    X_perm[, intervention_col] <- sample(X[, intervention_col])
    t_perm <- fit_tmap(Y, X_perm, intervention_col)
    clusters_perm <- find_clusters(pts$x, pts$y, t_perm, t_crit)
    null_max_mass[perm] <- if (length(clusters_perm)) max(vapply(clusters_perm, function(cl) cl$mass, numeric(1))) else 0
  }

  cluster_rows <- lapply(seq_along(clusters_obs), function(i) {
    cl <- clusters_obs[[i]]
    p_corrected <- (1 + sum(null_max_mass >= cl$mass)) / (1 + opt$`n-permutations`)
    data.frame(
      hemisphere = hemi,
      cluster_id = i,
      polarity = cl$polarity,
      size = cl$size,
      mass = cl$mass,
      mean_t = mean(t_obs[cl$idx]),
      centroid_x = mean(pts$x[cl$idx]),
      centroid_y = mean(pts$y[cl$idx]),
      p_corrected = p_corrected,
      significant = p_corrected < opt$alpha,
      stringsAsFactors = FALSE
    )
  })
  cluster_df <- do.call(rbind, cluster_rows)
  cluster_df <- cluster_df[order(cluster_df$p_corrected), ]
  write.csv(cluster_df, file.path(opt$outdir, paste0(hemi, "_clusters.csv")), row.names=FALSE)
  overall_cluster_rows[[hemi]] <- cluster_df

  # per-cluster member point lists (for mapping back onto the surface)
  members_df <- do.call(rbind, lapply(seq_along(clusters_obs), function(i) {
    cl <- clusters_obs[[i]]
    data.frame(cluster_id=i, polarity=cl$polarity, x=pts$x[cl$idx], y=pts$y[cl$idx], t=t_obs[cl$idx])
  }))
  write.csv(members_df, file.path(opt$outdir, paste0(hemi, "_cluster_members.csv")), row.names=FALSE)

  msg("%s: %d/%d clusters significant at p_corrected < %.2f\n", hemi, sum(cluster_df$significant), nrow(cluster_df), opt$alpha)
}

# ------------------------------------------------------------------------
# Combined summary
# ------------------------------------------------------------------------
all_clusters <- do.call(rbind, overall_cluster_rows)
con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf(
  "Hippocampal thickness vertex-wise cluster-based permutation analysis\nChange score: %s -> %s\nCluster-forming p < %.3f, %d permutations, cluster-level alpha = %.2f\n\n",
  baseline_ses, final_ses, opt$`cluster-forming-p`, opt$`n-permutations`, opt$alpha
), file=con)
if (!is.null(all_clusters) && nrow(all_clusters)) {
  cat(sprintf("%d total clusters across hemispheres, %d significant.\n\n", nrow(all_clusters), sum(all_clusters$significant)), file=con)
  for (i in seq_len(nrow(all_clusters))) {
    row <- all_clusters[i, ]
    cat(sprintf("%s cluster %d [%s]: size=%d points, mass=%.2f, mean_t=%.2f, centroid=(x=%.1f,y=%.1f), p_corrected=%.4f%s\n",
                row$hemisphere, row$cluster_id, row$polarity, row$size, row$mass, row$mean_t,
                row$centroid_x, row$centroid_y, row$p_corrected, if (row$significant) " *" else ""), file=con)
  }
} else {
  cat("No suprathreshold clusters found in either hemisphere.\n", file=con)
}
msg("\nDone. Per-hemisphere t-maps, cluster tables, and cluster member point lists are in %s/\n", opt$outdir)
msg("See %s/summary.txt for the combined report.\n", opt$outdir)
