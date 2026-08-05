#!/usr/bin/env Rscript
#
# Multivariate pattern analysis for the amygdala, closing a coverage gap:
# hippocampal volume (06), Destrieux cortical thickness (28), and
# hippocampal thickness (31) each have a committed MANOVA+PCA(+PLS-DA)
# script; amygdala never did (it only ever got a one-off scratch analysis).
# This mirrors 06_manova_pca_pattern.R's MANOVA/PCA design exactly (same
# change-score construction, same covariate-residualized PCA, same FDR
# handling), plus the PLS-DA leave-one-out classification check that 28/31
# added on top of 06's original two methods, for full parity across domains.
#
#   1. MANOVA (Pillai's trace) on the subject x nucleus change-score matrix
#      (baseline -> final session, hemispheres summed), covarying
#      age_z (rank-transformed) + sex.
#   2. PCA on the covariate-residualized change-score matrix -- latent
#      components of coordinated change across nuclei, tested against group.
#   3. PLS-DA (leave-one-out cross-validated) -- out-of-sample check of
#      whether the whole-amygdala nucleus pattern discriminates intervention
#      from control.

suppressPackageStartupMessages({
  library(optparse)
  library(pls)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to amygdala_tidy_2tp_fs82.tsv"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/40_amygdala_manova_pca_pattern",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character", default="all",
              help="Comma-separated nucleus names, or 'all' [default %default]"),
  make_option(c("--intervention-groups"), type="character", default="Single_2wk,Single_4wk,Group_2wk,Group_4wk"),
  make_option(c("--control-group"), type="character", default="Control"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--final-session"), type="character", default=NULL),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--seed"), type="integer", default=129),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }
set.seed(opt$seed)

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

msg("Loading tidy amygdala data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_set <- if (identical(opt$`roi-set`, "all")) sort(unique(tidy$nucleus)) else trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_dat <- tidy[tidy$nucleus %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check nucleus names against the tidy file")

agg <- aggregate(volume ~ subject_id + session + nucleus, data=roi_dat, FUN=sum)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Change score: %s -> %s\n", baseline_ses, final_ses)

base_vals <- agg[agg$session == baseline_ses, c("subject_id","nucleus","volume")]; names(base_vals)[3] <- "baseline"
final_vals <- agg[agg$session == final_ses, c("subject_id","nucleus","volume")]; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id","nucleus"))
merged_long$change <- log(merged_long$final) - log(merged_long$baseline)

change_wide <- reshape(merged_long[, c("subject_id","nucleus","change")], idvar="subject_id", timevar="nucleus", direction="wide")
change_cols <- setdiff(names(change_wide), "subject_id")
names(change_wide)[names(change_wide) %in% change_cols] <- sub("^change\\.", "", change_cols)
change_cols <- setdiff(names(change_wide), "subject_id")

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

dat <- merge(change_wide, participants[, required_pcols], by="subject_id")
dat <- dat[dat$group %in% c(intervention_groups, control_group), , drop=FALSE]
dat$intervention <- factor(ifelse(dat$group %in% intervention_groups, "intervention", "control"), levels=c("control","intervention"))
dat$sex <- factor(dat$sex)
dat$age_z <- as.numeric(scale(rank(dat$age)))

n_before <- nrow(dat)
dat <- dat[complete.cases(dat[, c(change_cols, "intervention","age_z","sex")]), , drop=FALSE]
n <- nrow(dat)
msg("Subjects with complete baseline+final data across all %d nuclei: %d of %d\n", length(change_cols), n, n_before)
if (n < length(change_cols) + 10) {
  warning(sprintf("Only %d complete subjects for %d nucleus features -- MANOVA/PCA results at this ratio should be treated cautiously", n, length(change_cols)))
}

write.csv(dat[, c("subject_id","intervention",change_cols)], file.path(opt$outdir, "change_score_matrix.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# 1. MANOVA
# ------------------------------------------------------------------------
msg("Running MANOVA on joint nucleus-change pattern...\n")
Y <- as.matrix(dat[, change_cols])
manova_fit <- manova(Y ~ intervention + age_z + sex, data=dat)
manova_summary <- as.data.frame(summary(manova_fit, test="Pillai")$stats)
manova_summary$term <- rownames(manova_summary)
write.csv(manova_summary, file.path(opt$outdir, "manova_results.csv"), row.names=FALSE)
intervention_manova_p <- manova_summary[manova_summary$term == "intervention", "Pr(>F)"]
msg("MANOVA (Pillai's trace) for 'intervention' across %d nuclei jointly: p = %s\n",
    length(change_cols), if(length(intervention_manova_p)) format(intervention_manova_p, digits=4) else "NA")

# ------------------------------------------------------------------------
# 2. PCA on covariate-residualized change scores
# ------------------------------------------------------------------------
msg("Running PCA on covariate-residualized change-score matrix...\n")
resid_mat <- sapply(change_cols, function(col) residuals(lm(as.formula(paste0("`", col, "` ~ age_z + sex")), data=dat)))
colnames(resid_mat) <- change_cols
pca_fit <- prcomp(resid_mat, center=TRUE, scale.=TRUE)
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
n_keep <- max(1, sum(pca_fit$sdev^2 > 1))
msg("PCA: %d components retained by Kaiser criterion (eigenvalue > 1) out of %d\n", n_keep, length(change_cols))

loadings_df <- as.data.frame(pca_fit$rotation[, seq_len(n_keep), drop=FALSE])
loadings_df$nucleus <- rownames(loadings_df)
write.csv(loadings_df, file.path(opt$outdir, "pca_loadings.csv"), row.names=FALSE)
write.csv(data.frame(component=paste0("PC", seq_along(var_explained)), variance_explained=var_explained),
          file.path(opt$outdir, "pca_variance_explained.csv"), row.names=FALSE)

pc_scores <- as.data.frame(pca_fit$x[, seq_len(n_keep), drop=FALSE])
pc_scores$intervention <- dat$intervention
pc_test_rows <- list()
for (i in seq_len(n_keep)) {
  pc_name <- paste0("PC", i)
  fit <- lm(pc_scores[[pc_name]] ~ intervention, data=pc_scores)
  tt <- summary(fit)$coefficients
  pc_test_rows[[pc_name]] <- data.frame(
    component = pc_name, variance_explained = var_explained[i],
    intervention_estimate = if ("interventionintervention" %in% rownames(tt)) tt["interventionintervention", "Estimate"] else NA_real_,
    intervention_p = if ("interventionintervention" %in% rownames(tt)) tt["interventionintervention", "Pr(>|t|)"] else NA_real_,
    stringsAsFactors = FALSE
  )
}
pc_test_df <- do.call(rbind, pc_test_rows)
pc_test_df$intervention_p_fdr <- p.adjust(pc_test_df$intervention_p, method="fdr")
pc_test_df$significant <- pc_test_df$intervention_p_fdr < opt$alpha
write.csv(pc_test_df, file.path(opt$outdir, "pca_component_group_tests.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# 3. PLS-DA: leave-one-out cross-validated classification (parity with 28/31)
# ------------------------------------------------------------------------
msg("Running PLS-DA (leave-one-out) classification check...\n")
y_bin <- as.numeric(dat$intervention) - 1
X <- scale(Y)
ncomp <- min(5, length(change_cols) - 1, n - 2)
plsda_fit <- tryCatch(plsr(y_bin ~ X, ncomp=ncomp, validation="LOO"), error=function(e) NULL)
loo_acc <- NA_real_; loo_auc <- NA_real_; best_ncomp <- NA_integer_
if (!is.null(plsda_fit)) {
  press <- as.numeric(plsda_fit$validation$PRESS)
  best_ncomp <- which.min(press)
  loo_pred <- plsda_fit$validation$pred[, 1, best_ncomp]
  loo_class <- as.numeric(loo_pred > 0.5)
  loo_acc <- mean(loo_class == y_bin)
  n_pos <- sum(y_bin == 1); n_neg <- sum(y_bin == 0)
  r <- rank(loo_pred)
  loo_auc <- (sum(r[y_bin == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
chance_acc <- max(mean(y_bin==1), mean(y_bin==0))
plsda_df <- data.frame(best_ncomp=best_ncomp, loo_accuracy=loo_acc, loo_auc=loo_auc, chance_accuracy=chance_acc)
write.csv(plsda_df, file.path(opt$outdir, "plsda_loo_results.csv"), row.names=FALSE)

# ------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------
con <- file(file.path(opt$outdir, "summary.txt"), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf(
  "Amygdala multivariate pattern analysis: %s -> %s change scores across %d nuclei\n\nSubjects: %d\n\n1. MANOVA (Pillai's trace): joint test of intervention vs control across the nucleus pattern\n   p = %s\n\n2. PCA (Kaiser criterion, eigenvalue > 1): %d components retained, %.1f%% variance explained\n   Per-component group tests (FDR-corrected across components):\n",
  baseline_ses, final_ses, length(change_cols), n,
  if(length(intervention_manova_p)) format(intervention_manova_p, digits=4) else "NA",
  n_keep, 100*sum(var_explained[seq_len(n_keep)])
), file=con)
for (i in seq_len(nrow(pc_test_df))) {
  row <- pc_test_df[i, ]
  cat(sprintf("   - %s (%.1f%% var): intervention effect = %.4f, p = %.4f, p_fdr = %.4f%s\n",
              row$component, 100*row$variance_explained, row$intervention_estimate, row$intervention_p, row$intervention_p_fdr,
              if (row$significant) " *" else ""), file=con)
}
cat(sprintf("\n3. PLS-DA (leave-one-out, %d components): LOO accuracy = %.3f (chance = %.3f), LOO AUC = %.3f\n",
            best_ncomp, loo_acc, chance_acc, loo_auc), file=con)
cat("\nSee pca_loadings.csv for which nuclei drive each component.\n", file=con)

msg("Done. See %s/summary.txt for a plain-text overview.\n", opt$outdir)
