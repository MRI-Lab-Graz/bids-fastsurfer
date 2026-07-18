#!/usr/bin/env Rscript
#
# CONDITIONAL / OPTIONAL exploratory analysis: random forest classification
# of dance vs control from subfield change scores.
#
# Run this ONLY if 01_primary_lmm.R shows a signal. It is explicitly
# hypothesis-generating, not confirmatory: at the sample sizes typical of
# this kind of study (n ~ 60-120), cross-validated classification accuracy
# is unstable and prone to overfitting with many correlated features. The
# sibling 129run project ran 11 ML algorithms at comparable n and got
# 44-60% accuracy (barely above chance) -- the realistic prior here is that
# this analysis will be uninformative, and that should be stated plainly in
# any report rather than mined for a flattering number.
#
# What this script deliberately does NOT do: hyperparameter search, testing
# many algorithms, or reporting the best of several accuracy numbers. One
# algorithm (random forest), one nested CV estimate, permutation test for
# whether accuracy exceeds chance.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/04_exploratory_rf",
              help="Output directory [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--baseline-session"), type="character", default=NULL,
              help="Session value to use as baseline for change scores [default: earliest session present]"),
  make_option(c("--final-session"), type="character", default=NULL,
              help="Session value to use as endpoint for change scores [default: latest session present]"),
  make_option(c("--n-permutations"), type="integer", default=1000,
              help="Permutations for the chance-accuracy null distribution [default %default]"),
  make_option(c("--seed"), type="integer", default=134, help="Random seed [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (!requireNamespace("randomForest", quietly=TRUE)) stop("Package 'randomForest' is required")
if (!requireNamespace("caret", quietly=TRUE)) stop("Package 'caret' is required")
suppressPackageStartupMessages({
  library(randomForest)
  library(caret)
})

set.seed(opt$seed)
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$is_composite <- as.logical(tidy$is_composite)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
tidy <- tidy[!tidy$is_composite & tidy$subfield_raw != "hippocampal-fissure", , drop=FALSE]

sessions <- sort(unique(tidy$session))
baseline_ses <- if (!is.null(opt$`baseline-session`)) opt$`baseline-session` else sessions[1]
final_ses <- if (!is.null(opt$`final-session`)) opt$`final-session` else sessions[length(sessions)]
msg("Change score: %s -> %s\n", baseline_ses, final_ses)

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=tidy, FUN=sum)
agg$feature <- paste(agg$hemisphere, agg$subfield, sep="_")

wide_base <- reshape(agg[agg$session == baseline_ses, c("subject_id","feature","volume")],
                      idvar="subject_id", timevar="feature", direction="wide")
names(wide_base) <- gsub("^volume\\.", "baseline_", names(wide_base))
wide_final <- reshape(agg[agg$session == final_ses, c("subject_id","feature","volume")],
                       idvar="subject_id", timevar="feature", direction="wide")
names(wide_final) <- gsub("^volume\\.", "final_", names(wide_final))

merged <- merge(wide_base, wide_final, by="subject_id")
feature_names <- unique(agg$feature)
for (f in feature_names) {
  bcol <- paste0("baseline_", f); fcol <- paste0("final_", f); ccol <- paste0("change_", f)
  if (bcol %in% names(merged) && fcol %in% names(merged)) {
    merged[[ccol]] <- merged[[fcol]] - merged[[bcol]]
  }
}
change_cols <- grep("^change_", names(merged), value=TRUE)
merged <- merged[, c("subject_id", change_cols)]

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
merged <- merge(merged, participants[, c("subject_id","group")], by="subject_id")
merged <- merged[merged$group %in% c(dance_groups, control_group), , drop=FALSE]
merged$dance <- factor(ifelse(merged$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
merged$group <- NULL

n_complete <- sum(complete.cases(merged))
msg("Subjects with complete baseline+final change scores: %d of %d\n", n_complete, nrow(merged))
merged <- merged[complete.cases(merged), , drop=FALSE]
if (nrow(merged) < 20) {
  warning(sprintf("Only %d complete subjects -- random forest results at this n should be treated as unreliable", nrow(merged)))
}

X <- merged[, change_cols, drop=FALSE]
y <- merged$dance

# ------------------------------------------------------------------------
# Nested cross-validation: outer loop estimates generalization accuracy,
# inner loop only tunes mtry (no algorithm search, no feature selection
# search -- both would inflate the estimate at this sample size).
# ------------------------------------------------------------------------
msg("Running nested cross-validation (random forest, single algorithm)...\n")
outer_folds <- caret::createFolds(y, k=5, list=TRUE)
fold_acc <- numeric(length(outer_folds))
importance_list <- list()

for (i in seq_along(outer_folds)) {
  test_idx <- outer_folds[[i]]
  train_X <- X[-test_idx, , drop=FALSE]; train_y <- y[-test_idx]
  test_X <- X[test_idx, , drop=FALSE]; test_y <- y[test_idx]

  ctrl <- caret::trainControl(method="cv", number=5)
  tuned <- caret::train(x=train_X, y=train_y, method="rf", trControl=ctrl,
                         tuneLength=3, ntree=500)
  pred <- predict(tuned, newdata=test_X)
  fold_acc[i] <- mean(pred == test_y)
  importance_list[[i]] <- caret::varImp(tuned)$importance
}

observed_accuracy <- mean(fold_acc)
msg("Observed 5-fold nested CV accuracy: %.3f (chance = %.3f, majority class)\n",
    observed_accuracy, max(prop.table(table(y))))

# ------------------------------------------------------------------------
# Permutation test: shuffle labels, refit, see how often permuted accuracy
# matches or exceeds the observed accuracy
# ------------------------------------------------------------------------
msg("Running %d label permutations for a null accuracy distribution...\n", opt$`n-permutations`)
perm_acc <- numeric(opt$`n-permutations`)
for (p in seq_len(opt$`n-permutations`)) {
  y_perm <- sample(y)
  folds_p <- caret::createFolds(y_perm, k=5, list=TRUE)
  acc_p <- numeric(length(folds_p))
  for (i in seq_along(folds_p)) {
    test_idx <- folds_p[[i]]
    rf_p <- randomForest::randomForest(x=X[-test_idx, , drop=FALSE], y=y_perm[-test_idx], ntree=200)
    acc_p[i] <- mean(predict(rf_p, X[test_idx, , drop=FALSE]) == y_perm[test_idx])
  }
  perm_acc[p] <- mean(acc_p)
}
perm_p_value <- mean(perm_acc >= observed_accuracy)

# ------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------
importance_df <- do.call(rbind, lapply(seq_along(importance_list), function(i) {
  df <- importance_list[[i]]; df$feature <- rownames(df); df$fold <- i; df
}))
write.csv(importance_df, file.path(opt$outdir, "feature_importance_by_fold.csv"), row.names=FALSE)
write.csv(data.frame(fold=seq_along(fold_acc), accuracy=fold_acc), file.path(opt$outdir, "fold_accuracy.csv"), row.names=FALSE)
write.csv(data.frame(permutation=seq_along(perm_acc), accuracy=perm_acc), file.path(opt$outdir, "permutation_null_accuracy.csv"), row.names=FALSE)

summary_txt <- sprintf(
  "Exploratory random forest classification (dance vs control), change scores %s -> %s\n\nSubjects: %d\nFeatures: %d\nObserved 5-fold nested CV accuracy: %.3f\nChance (majority class): %.3f\nPermutation test p-value (accuracy >= observed under label shuffling, %d permutations): %.4f\n\nInterpretation: this analysis is EXPLORATORY/hypothesis-generating.\nGiven the sample size, treat this result as unreliable unless the\npermutation p-value is clearly significant AND the primary confirmatory\nanalysis (01_primary_lmm.R) also shows a signal.\n",
  baseline_ses, final_ses, nrow(merged), length(change_cols),
  observed_accuracy, max(prop.table(table(y))), opt$`n-permutations`, perm_p_value
)
writeLines(summary_txt, file.path(opt$outdir, "summary.txt"))
cat(summary_txt)
