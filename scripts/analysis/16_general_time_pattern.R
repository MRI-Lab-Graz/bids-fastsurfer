#!/usr/bin/env Rscript
#
# General time-pattern analysis: now that the combined dance+running vs.
# control test is a clean, well-powered null (n~205, LOSO-robust,
# consistent across studies), a different question becomes worth asking
# with this much power -- is there a GENERAL temporal pattern in subfield
# volume, present regardless of group (i.e. a population-level trajectory
# over the study window, not an intervention effect)? A true group effect
# and a general time effect are different things a null on 15 does not
# rule out.
#
# Three tests per ROI, ignoring group (dance is only ever a covariate,
# never the term of interest here):
#   1. time_effect: does log(volume) change from baseline to final
#      session, pooling everyone? (one-sample-style test via the
#      followup_f fixed effect in a no-dance-interaction model)
#   2. Joint/multivariate: does the VECTOR of per-subject change scores
#      across all pre-specified ROIs jointly differ from zero (Hotelling's
#      one-sample T^2 via MANOVA), i.e. is there a coordinated cross-region
#      pattern rather than independent per-ROI noise?
#   3. PCA on the change-score matrix: purely descriptive characterisation
#      of the general trajectory SHAPE (which regions move together, which
#      move in opposite directions) -- not tested against group, just
#      described.
#
# LOSO robustness is run on any ROI whose time_effect is nominally
# interesting, same practice as everywhere else in this project.

suppressPackageStartupMessages({
  library(optparse)
  library(lme4)
  library(lmerTest)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to combined tidy volumes TSV"),
  make_option(c("-p", "--participants"), type="character", help="Path to combined participants TSV"),
  make_option(c("-o", "--outdir"), type="character", default="results/16_general_time_pattern"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,GC-ML-DG,CA1,CA3,CA4,subiculum,molecular_layer_HP"),
  make_option(c("--roi-column"), type="character", default="subfield"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--final-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy) || is.null(opt$participants)) stop("--tidy and --participants are required")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_col <- opt$`roi-column`

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
if (!"study" %in% names(tidy)) stop("tidy file must have a 'study' column")
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set")

agg <- aggregate(as.formula(paste("volume ~ subject_id + session + hemisphere + study +", roi_col)), data=roi_dat, FUN=sum)
names(agg)[names(agg) == roi_col] <- "subfield"
etiv_lookup <- unique(tidy[, c("subject_id","etiv")]); etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
agg <- merge(agg, participants[, c("subject_id","group","age","sex")], by="subject_id")
agg$dance <- factor(ifelse(agg$group == "intervention", "intervention", "control"), levels=c("control","intervention"))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$study <- factor(agg$study)
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
followup_sessions <- setdiff(sessions, baseline_ses)
msg("Subjects: %d. Baseline: %s. Follow-ups: %s\n", length(unique(agg$subject_id)), baseline_ses, paste(followup_sessions, collapse=", "))

build_change_data <- function(d_roi) {
  base <- d_roi[d_roi$session == baseline_ses, c("subject_id","hemisphere","log_volume")]
  names(base)[3] <- "log_baseline"
  fu <- d_roi[d_roi$session %in% followup_sessions, , drop=FALSE]
  merged <- merge(fu, base, by=c("subject_id","hemisphere"))
  merged$log_change <- merged$log_volume - merged$log_baseline
  merged$log_baseline_z <- as.numeric(scale(merged$log_baseline))
  merged$followup_f <- factor(merged$session, levels=followup_sessions)
  merged
}

# ------------------------------------------------------------------------
# 1. Per-ROI general time effect (dance/study as covariates, not tested)
# ------------------------------------------------------------------------
msg("\n--- Test 1: general time effect per ROI ---\n")
time_formula <- log_change ~ followup_f + dance + study + log_baseline_z + hemisphere + age_z + sex + etiv_z
time_rows <- list()
change_by_roi <- list()

for (roi_name in unique(agg$subfield)) {
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  change_dat <- build_change_data(d_roi)
  change_by_roi[[roi_name]] <- change_dat
  model <- tryCatch(lmerTest::lmer(update(time_formula, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE), error=function(e) NULL)
  if (is.null(model)) next
  at <- tryCatch(anova(model), error=function(e) NULL)
  p_val <- if (!is.null(at) && "followup_f" %in% rownames(at)) at["followup_f", "Pr(>F)"] else NA_real_
  fe <- lme4::fixef(model)
  est <- if ("followup_fses-3" %in% names(fe)) fe[["followup_fses-3"]] else if (length(followup_sessions)==1) fe[[paste0("followup_f", followup_sessions[1])]] else NA_real_
  time_rows[[roi_name]] <- data.frame(roi=roi_name, n_obs=nrow(change_dat), time_estimate=est, time_p=p_val, stringsAsFactors=FALSE)
  msg("  %s: estimate=%.5f, p=%.4f\n", roi_name, est, p_val)
}
time_df <- do.call(rbind, time_rows)
time_df$time_p_fdr <- p.adjust(time_df$time_p, method="fdr")
time_df$significant <- time_df$time_p_fdr < opt$alpha
time_df <- time_df[order(time_df$time_p_fdr), ]
write.csv(time_df, file.path(opt$outdir, "general_time_effect.csv"), row.names=FALSE)

# LOSO for nominally interesting time effects
screen <- time_df[!is.na(time_df$time_p) & time_df$time_p < 0.10, ]
loso_time_rows <- list()
if (nrow(screen)) {
  msg("\nRunning LOSO for %d general-time hits...\n", nrow(screen))
  for (roi_name in screen$roi) {
    change_dat <- change_by_roi[[roi_name]]
    full_model <- lmerTest::lmer(update(time_formula, ". ~ . + (1 | subject_id)"), data=change_dat, REML=TRUE)
    full_p <- anova(full_model)["followup_f", "Pr(>F)"]
    for (excl_subj in unique(change_dat$subject_id)) {
      d_sub <- change_dat[change_dat$subject_id != excl_subj, ]
      fit <- tryCatch(lmerTest::lmer(update(time_formula, ". ~ . + (1 | subject_id)"), data=d_sub, REML=TRUE), error=function(e) NULL)
      p_val <- if (!is.null(fit)) { at <- tryCatch(anova(fit), error=function(e) NULL); if (!is.null(at) && "followup_f" %in% rownames(at)) at["followup_f","Pr(>F)"] else NA_real_ } else NA_real_
      loso_time_rows[[length(loso_time_rows)+1]] <- data.frame(roi=roi_name, excluded_subject=excl_subj, loso_p=p_val, full_p=full_p, stringsAsFactors=FALSE)
    }
  }
  loso_time_df <- do.call(rbind, loso_time_rows)
  write.csv(loso_time_df, file.path(opt$outdir, "loso_general_time.csv"), row.names=FALSE)
}

# ------------------------------------------------------------------------
# 2 & 3. Joint pattern across ROIs: one-sample MANOVA + PCA on the
# subject x ROI change-score matrix (bilateral-summed, ses baseline->final)
# ------------------------------------------------------------------------
msg("\n--- Test 2/3: joint cross-region pattern ---\n")
group_cols <- c("subject_id","session","subfield")
agg_sum <- aggregate(as.formula(paste("volume ~ subject_id + session +", "subfield")),
                      data=agg[, c("subject_id","session","subfield","volume")], FUN=sum)
base_vals <- agg_sum[agg_sum$session == baseline_ses, c("subject_id","subfield","volume")]
final_vals <- agg_sum[agg_sum$session == final_ses, c("subject_id","subfield","volume")]
names(base_vals)[3] <- "baseline"; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id","subfield"))
merged_long$change <- log(merged_long$final) - log(merged_long$baseline)
change_wide <- reshape(merged_long[, c("subject_id","subfield","change")], idvar="subject_id", timevar="subfield", direction="wide")
change_cols <- setdiff(names(change_wide), "subject_id")
names(change_wide)[names(change_wide) %in% change_cols] <- sub("^change\\.", "", change_cols)
change_cols <- setdiff(names(change_wide), "subject_id")
change_wide <- change_wide[complete.cases(change_wide[, change_cols]), ]
msg("Subjects with complete change data across all %d ROIs: %d\n", length(change_cols), nrow(change_wide))

Y <- as.matrix(change_wide[, change_cols])
# One-sample Hotelling's T^2: tests whether the mean change VECTOR
# (jointly across ROIs) differs from zero -- i.e. is there a real,
# coordinated population-level trajectory, not just per-ROI noise.
n <- nrow(Y); p <- ncol(Y)
Ybar <- colMeans(Y)
S <- cov(Y)
T2 <- n * t(Ybar) %*% solve(S) %*% Ybar
Fstat <- (n - p) / (p * (n - 1)) * T2
pval_joint <- pf(Fstat, p, n - p, lower.tail = FALSE)
msg("Hotelling's one-sample T^2 (mean change vector != 0, jointly across %d ROIs): F=%.3f, p=%s\n",
    p, Fstat, format(pval_joint, digits=4))

writeLines(sprintf("Hotelling's one-sample T^2: n=%d, p_roi=%d, F=%.4f, df1=%d, df2=%d, p=%s",
                    n, p, Fstat, p, n-p, format(pval_joint, digits=6)),
           file.path(opt$outdir, "joint_pattern_test.txt"))

# Per-ROI mean change (for interpreting the joint test's direction)
mean_change_df <- data.frame(roi=change_cols, mean_change=colMeans(Y), sd_change=apply(Y,2,sd))
mean_change_df$t <- mean_change_df$mean_change / (mean_change_df$sd_change/sqrt(n))
mean_change_df$p_uncorrected <- 2*pt(-abs(mean_change_df$t), df=n-1)
mean_change_df$p_fdr <- p.adjust(mean_change_df$p_uncorrected, method="fdr")
mean_change_df <- mean_change_df[order(mean_change_df$p_fdr), ]
write.csv(mean_change_df, file.path(opt$outdir, "per_roi_mean_change.csv"), row.names=FALSE)

# PCA: descriptive shape of coordinated change
pca_fit <- prcomp(Y, center=TRUE, scale.=TRUE)
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
n_keep <- max(1, sum(pca_fit$sdev^2 > 1))
loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop=FALSE])
loadings_df$roi <- rownames(loadings_df)
write.csv(loadings_df, file.path(opt$outdir, "pca_loadings.csv"), row.names=FALSE)
write.csv(data.frame(component=paste0("PC",seq_along(var_explained)), variance_explained=var_explained),
          file.path(opt$outdir, "pca_variance_explained.csv"), row.names=FALSE)

con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf("General time-pattern analysis (%s -> %s, n=%d, pooled across groups/studies)\n\n", baseline_ses, final_ses, n), file=con)
cat("1. Per-ROI general time effect (mixed model, dance/study as covariates):\n", file=con)
for (i in seq_len(nrow(time_df))) {
  row <- time_df[i,]
  cat(sprintf("   %s: estimate=%.5f, p=%.4f, p_fdr=%.4f%s\n", row$roi, row$time_estimate, row$time_p, row$time_p_fdr, if(row$significant)" *" else ""), file=con)
}
cat(sprintf("\n2. Joint cross-region pattern (Hotelling's T^2): F=%.3f, p=%s%s\n", Fstat, format(pval_joint,digits=4), if(pval_joint<opt$alpha) " *" else ""), file=con)
cat("\n3. Per-ROI mean change (simple one-sample t-test, descriptive):\n", file=con)
for (i in seq_len(nrow(mean_change_df))) {
  row <- mean_change_df[i,]
  cat(sprintf("   %s: mean=%.5f, p=%.4f, p_fdr=%.4f%s\n", row$roi, row$mean_change, row$p_uncorrected, row$p_fdr, if(row$p_fdr<opt$alpha)" *" else ""), file=con)
}
cat(sprintf("\n4. PCA: %d components retained (Kaiser), %.1f%% variance explained by PC1-%d\n", n_keep, 100*sum(var_explained[seq_len(n_keep)]), n_keep), file=con)
cat("   See pca_loadings.csv -- same-sign loadings = regions moving together;\n   mixed-sign = a contrast pattern (some growing while others shrink).\n", file=con)

msg("\nDone. See %s/summary.txt\n", opt$outdir)
