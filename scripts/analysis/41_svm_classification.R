#!/usr/bin/env Rscript
#
# Generic SVM classification check (intervention vs control), age/sex-
# corrected, for any domain's tidy table. Closes two gaps at once:
#   1. The original SVM check (results/plots/svm_check.R, hippocampal
#      volume only) never controlled for age/sex -- it classified on raw
#      change scores. This version residualizes every feature against
#      age_z (rank-transformed) + sex BEFORE classification, exactly as
#      06/28/31/40's PCA steps already do, so a "hit" here can't just be
#      the classifier rediscovering an age or sex imbalance between groups.
#   2. SVM was never run at all for amygdala, Destrieux/cortical thickness,
#      or hippocampal thickness (hipsta) -- this is domain-agnostic via
#      --roi-col/--value-col/--hemisphere-agg, so one script covers all of
#      them (see run_pipeline.sh for the per-domain invocations).
#
# Method (same nested-CV + permutation-test design as the original
# svm_check.R, generalized): 5-fold outer CV for the generalization
# estimate, e1071::tune.svm (5-fold inner CV) for cost tuning, linear AND
# radial kernels, then a label-permutation null (linear kernel) for a
# significance test on the observed accuracy -- avoids caret's
# svmLinear/svmRadial wrappers, which pull in the kernlab package.

suppressPackageStartupMessages({
  library(optparse)
  library(caret)
  library(e1071)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character", help="Path to a domain tidy TSV"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("--roi-col"), type="character", default="subfield"),
  make_option(c("--value-col"), type="character", default="volume"),
  make_option(c("--roi-set"), type="character", default="all"),
  make_option(c("--hemisphere-agg"), type="character", default="sum", help="'sum' or 'mean' [default %default]"),
  make_option(c("--keep-hemisphere-separate"), action="store_true", default=FALSE,
              help="Keep lh/rh as separate features instead of aggregating [default FALSE]"),
  make_option(c("--log-change"), action="store_true", default=FALSE,
              help="Use log(final)-log(baseline) instead of a raw difference (volume domains)"),
  make_option(c("--intervention-groups"), type="character", default="Single_2wk,Single_4wk,Group_2wk,Group_4wk"),
  make_option(c("--control-group"), type="character", default="Control"),
  make_option(c("--baseline-session"), type="character", default=NULL),
  make_option(c("--final-session"), type="character", default=NULL),
  make_option(c("--n-permutations"), type="integer", default=200),
  make_option(c("--n-folds"), type="integer", default=5),
  make_option(c("-o", "--outdir"), type="character", default="results/41_svm_classification"),
  make_option(c("--label"), type="character", default="domain"),
  make_option(c("--seed"), type="integer", default=134),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }
set.seed(opt$seed)

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-col`; value_col <- opt$`value-col`
intervention_groups <- trimws(strsplit(opt$`intervention-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy[[value_col]] <- suppressWarnings(as.numeric(tidy[[value_col]]))
roi_set <- if (identical(opt$`roi-set`, "all")) sort(unique(tidy[[roi_col]])) else trimws(strsplit(opt$`roi-set`, ",")[[1]])
roi_dat <- tidy[tidy[[roi_col]] %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set")

group_cols <- if (isTRUE(opt$`keep-hemisphere-separate`)) c("subject_id","session","hemisphere",roi_col) else c("subject_id","session",roi_col)
agg_fun <- if (identical(opt$`hemisphere-agg`, "mean")) mean else sum
agg <- aggregate(as.formula(paste(value_col, "~ .")), data=roi_dat[, c(group_cols, value_col)], FUN=agg_fun)

sessions <- sort(unique(agg$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("[%s] Change score: %s -> %s\n", opt$label, baseline_ses, final_ses)

feature_col <- if (isTRUE(opt$`keep-hemisphere-separate`)) {
  agg$feature <- paste(agg$hemisphere, agg[[roi_col]], sep="_")
  "feature"
} else roi_col

base_vals <- agg[agg$session == baseline_ses, c("subject_id", feature_col, value_col)]; names(base_vals)[3] <- "baseline"
final_vals <- agg[agg$session == final_ses, c("subject_id", feature_col, value_col)]; names(final_vals)[3] <- "final"
merged_long <- merge(base_vals, final_vals, by=c("subject_id", feature_col))
merged_long$change <- if (isTRUE(opt$`log-change`)) log(merged_long$final) - log(merged_long$baseline) else merged_long$final - merged_long$baseline

change_wide <- reshape(merged_long[, c("subject_id", feature_col, "change")], idvar="subject_id", timevar=feature_col, direction="wide")
names(change_wide) <- sub("^change\\.", "", names(change_wide))
change_cols <- setdiff(names(change_wide), "subject_id")

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
participants$intervention <- factor(ifelse(participants$group %in% intervention_groups, "intervention",
                                     ifelse(participants$group == control_group, "control", NA)), levels=c("control","intervention"))
participants$sex <- factor(participants$sex)
participants$age_z <- as.numeric(scale(rank(participants$age)))

dat <- merge(change_wide, participants[, c("subject_id","intervention","age_z","sex")], by="subject_id")
dat <- dat[complete.cases(dat[, c(change_cols, "intervention","age_z","sex")]), , drop=FALSE]
n <- nrow(dat)
msg("[%s] %d features, n=%d complete subjects (%d intervention, %d control)\n",
    opt$label, length(change_cols), n, sum(dat$intervention=="intervention"), sum(dat$intervention=="control"))
if (n < length(change_cols) + 10) {
  warning(sprintf("[%s] only %d complete subjects for %d features -- classification at this ratio should be treated cautiously", opt$label, n, length(change_cols)))
}

# Age/sex correction: residualize each feature against age_z + sex BEFORE
# classification (matches the residualization already used for PCA
# elsewhere in this project), so the classifier can only use age/sex-
# independent variance in brain change.
X <- sapply(change_cols, function(col) residuals(lm(as.formula(paste0("`", col, "` ~ age_z + sex")), data=dat)))
colnames(X) <- change_cols
y <- dat$intervention

run_nested_cv_svm <- function(X, y, kernel, k) {
  outer_folds <- caret::createFolds(y, k=k, list=TRUE)
  fold_acc <- numeric(length(outer_folds))
  for (i in seq_along(outer_folds)) {
    test_idx <- outer_folds[[i]]
    train_X <- as.matrix(X[-test_idx, , drop=FALSE]); train_y <- y[-test_idx]
    test_X <- as.matrix(X[test_idx, , drop=FALSE]); test_y <- y[test_idx]
    tuned <- e1071::tune.svm(x=train_X, y=train_y, kernel=kernel,
                              cost=c(0.01, 0.1, 1, 10), scale=TRUE,
                              tunecontrol=e1071::tune.control(cross=min(5, length(train_y))))
    pred <- predict(tuned$best.model, test_X)
    fold_acc[i] <- mean(pred == test_y)
  }
  mean(fold_acc)
}

msg("[%s] Running nested CV (linear + radial kernels)...\n", opt$label)
observed_linear <- run_nested_cv_svm(X, y, "linear", opt$`n-folds`)
observed_radial <- run_nested_cv_svm(X, y, "radial", opt$`n-folds`)
chance <- max(prop.table(table(y)))
msg("[%s] Observed accuracy -- linear: %.3f, radial: %.3f (chance/majority class: %.3f)\n",
    opt$label, observed_linear, observed_radial, chance)

msg("[%s] Running %d-permutation null (linear kernel)...\n", opt$label, opt$`n-permutations`)
perm_acc <- numeric(opt$`n-permutations`)
for (p in seq_len(opt$`n-permutations`)) {
  y_perm <- sample(y)
  folds_p <- caret::createFolds(y_perm, k=opt$`n-folds`, list=TRUE)
  acc_p <- numeric(length(folds_p))
  for (i in seq_along(folds_p)) {
    test_idx <- folds_p[[i]]
    fit <- e1071::svm(x=as.matrix(X[-test_idx, , drop=FALSE]), y=y_perm[-test_idx], kernel="linear", scale=TRUE)
    acc_p[i] <- mean(predict(fit, as.matrix(X[test_idx, , drop=FALSE])) == y_perm[test_idx])
  }
  perm_acc[p] <- mean(acc_p)
}
perm_p <- (1 + sum(perm_acc >= observed_linear)) / (1 + opt$`n-permutations`)

result <- data.frame(
  label = opt$label, n = n, n_features = length(change_cols),
  chance_accuracy = chance,
  linear_accuracy = observed_linear, radial_accuracy = observed_radial,
  permutation_p = perm_p, n_permutations = opt$`n-permutations`
)
write.csv(result, file.path(opt$outdir, paste0(opt$label, "_svm_result.csv")), row.names=FALSE)

con <- file(file.path(opt$outdir, paste0(opt$label, "_summary.txt")), open="wt")
cat(sprintf(
  "SVM classification (intervention vs control), age/sex-corrected: %s\n%s -> %s change scores, %d features, n=%d (%d intervention, %d control)\n\nObserved 5-fold nested CV accuracy:\n  linear SVM: %.3f\n  radial SVM: %.3f\n  chance/majority class: %.3f\n\nPermutation test (linear SVM, %d permutations): p = %.4f\n",
  opt$label, baseline_ses, final_ses, length(change_cols), n, sum(y=="intervention"), sum(y=="control"),
  observed_linear, observed_radial, chance, opt$`n-permutations`, perm_p
), file=con)
close(con)

msg("[%s] Permutation p = %.4f. Done.\n", opt$label, perm_p)
