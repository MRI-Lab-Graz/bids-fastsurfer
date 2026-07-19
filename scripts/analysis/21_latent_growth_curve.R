#!/usr/bin/env Rscript
#
# Longitudinal SEM: a latent growth curve model (LGM) as a genuinely
# different modelling framework to cross-check the age x CA3 x intervention
# finding, using ALL THREE waves at once (intercept + slope factors) rather
# than reducing the trajectory to a single baseline->endpoint change score.
#
# Scope decision: only CA3 is modelled (the flagship ROI), not swept across
# the full subfield set -- this is a targeted cross-validation of an
# existing finding via a different method, not a new fishing expedition.
# With only 3 waves, a linear growth model (intercept + slope, no
# quadratic) is the most that's identified; residual variances are
# constrained equal across waves to leave 1 df for a model-fit test
# (freely estimating them would give a saturated, untestable model).
# Missing waves are handled via full-information maximum likelihood
# (lavaan's missing="fiml"), so subjects with 2 of 3 sessions still
# contribute rather than being dropped.
#
# The slope factor -- each subject's estimated linear rate of CA3 change --
# is regressed on age, group (intervention vs control), study (if
# combined), and their age x group interaction: this is the SEM-native
# analogue of "does rate of change depend on age x intervention", using
# the full repeated-measures structure instead of a two-point difference.

suppressPackageStartupMessages({
  library(optparse)
  library(lavaan)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to a tidy volumes TSV (hippo_subfields_tidy.tsv or combined_hippo_tidy.tsv)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV (subject_id, group, age, sex[, study])"),
  make_option(c("-o", "--outdir"), type="character", default="results/21_latent_growth_curve"),
  make_option(c("--roi"), type="character", default="CA3", help="Single ROI/subfield value to model [default %default]"),
  make_option(c("--roi-column"), type="character", default="subfield"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the intervention contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control"),
  make_option(c("--has-study-column"), action="store_true", default=FALSE,
              help="Set if --participants/--tidy include a 'study' column (combined-sample run)"),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_col <- opt$`roi-column`
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

msg("Loading tidy data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
if (!roi_col %in% names(tidy)) stop(sprintf("--roi-column '%s' not found in tidy file", roi_col))

roi_dat <- tidy[tidy[[roi_col]] == opt$roi, , drop=FALSE]
if (!nrow(roi_dat)) stop(sprintf("No rows matched --roi '%s'", opt$roi))

# Bilateral sum (lh + rh) per subject/session -- a single composite value
# per timepoint, as required by a univariate growth curve.
agg <- aggregate(volume ~ subject_id + session, data=roi_dat, FUN=sum)
agg$log_volume <- log(agg$volume)

sessions <- sort(unique(agg$session))
if (length(sessions) != 3) stop(sprintf("Expected exactly 3 sessions, found %d: %s", length(sessions), paste(sessions, collapse=", ")))
msg("Sessions found: %s\n", paste(sessions, collapse=", "))

wide <- reshape(agg[, c("subject_id","session","log_volume")], idvar="subject_id", timevar="session", direction="wide")
names(wide) <- sub("^log_volume\\.", "y_", names(wide))
y_cols <- paste0("y_", sessions)
names(wide)[match(y_cols, names(wide))] <- c("y1","y2","y3")

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
req_cols <- c("subject_id","group","age","sex")
if (opt$`has-study-column`) req_cols <- c(req_cols, "study")
missing_pcols <- setdiff(req_cols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

dat <- merge(wide, participants[, req_cols], by="subject_id")
dat <- dat[dat$group %in% c(dance_groups, control_group) | dat$group %in% c("intervention","control"), , drop=FALSE]
dat$dance01 <- ifelse(dat$group %in% dance_groups | dat$group == "intervention", 1, 0)
dat$age_z <- as.numeric(scale(dat$age))
dat$age_x_dance <- dat$age_z * dat$dance01
if (opt$`has-study-column`) dat$study01 <- ifelse(dat$study == sort(unique(dat$study))[1], 0, 1)

n_complete3 <- sum(complete.cases(dat[, c("y1","y2","y3")]))
n_any <- sum(rowSums(!is.na(dat[, c("y1","y2","y3")])) >= 1)
msg("Subjects with all 3 waves: %d. Subjects with >=1 wave (used via FIML): %d\n", n_complete3, n_any)

# ------------------------------------------------------------------------
# Unconditional growth model first (no predictors) -- establishes the
# baseline shape (mean intercept/slope, and whether slope VARIANCE is
# distinguishable from zero) before adding covariates.
# ------------------------------------------------------------------------
model_unconditional <- '
  i =~ 1*y1 + 1*y2 + 1*y3
  s =~ 0*y1 + 1*y2 + 2*y3
  y1 ~~ r*y1
  y2 ~~ r*y2
  y3 ~~ r*y3
'
fit0 <- lavaan::growth(model_unconditional, data=dat, missing="fiml")

sink(file.path(opt$outdir, paste0(opt$roi, "_unconditional_model.txt")))
cat("Unconditional growth model:", opt$roi, "\n\n")
print(summary(fit0, fit.measures=TRUE, standardized=TRUE))
sink()

# ------------------------------------------------------------------------
# Conditional model: slope (and intercept) regressed on age, dance,
# age x dance, and study (if combined)
# ------------------------------------------------------------------------
covariate_terms <- if (opt$`has-study-column`) "age_z + dance01 + age_x_dance + study01" else "age_z + dance01 + age_x_dance"
model_conditional <- sprintf('
  i =~ 1*y1 + 1*y2 + 1*y3
  s =~ 0*y1 + 1*y2 + 2*y3
  y1 ~~ r*y1
  y2 ~~ r*y2
  y3 ~~ r*y3
  i ~ %s
  s ~ %s
', covariate_terms, covariate_terms)
fit1 <- lavaan::growth(model_conditional, data=dat, missing="fiml")

sink(file.path(opt$outdir, paste0(opt$roi, "_conditional_model.txt")))
cat("Conditional growth model:", opt$roi, "\n\n")
print(summary(fit1, fit.measures=TRUE, standardized=TRUE))
sink()

pe <- parameterEstimates(fit1, standardized=TRUE)
slope_rows <- pe[pe$lhs == "s" & pe$op == "~", ]
intercept_rows <- pe[pe$lhs == "i" & pe$op == "~", ]
write.csv(pe, file.path(opt$outdir, paste0(opt$roi, "_parameter_estimates.csv")), row.names=FALSE)

fm0 <- fitMeasures(fit0, c("chisq","df","pvalue","cfi","rmsea"))
fm1 <- fitMeasures(fit1, c("chisq","df","pvalue","cfi","rmsea"))

con <- file(file.path(opt$outdir, paste0(opt$roi, "_summary.txt")), open="wt")
on.exit(close(con), add=TRUE)
cat(sprintf("Latent growth curve model: %s (bilateral sum, log volume)\n", opt$roi), file=con)
cat(sprintf("n with all 3 waves: %d, n with >=1 wave (FIML): %d\n\n", n_complete3, n_any), file=con)
cat("Unconditional model fit: chisq=", sprintf("%.2f", fm0["chisq"]), " df=", fm0["df"],
    " p=", sprintf("%.4f", fm0["pvalue"]), " CFI=", sprintf("%.3f", fm0["cfi"]),
    " RMSEA=", sprintf("%.3f", fm0["rmsea"]), "\n\n", sep="", file=con)
cat("Conditional model fit: chisq=", sprintf("%.2f", fm1["chisq"]), " df=", fm1["df"],
    " p=", sprintf("%.4f", fm1["pvalue"]), " CFI=", sprintf("%.3f", fm1["cfi"]),
    " RMSEA=", sprintf("%.3f", fm1["rmsea"]), "\n\n", sep="", file=con)
cat("Slope factor (rate of change) regressed on:\n", file=con)
for (i in seq_len(nrow(slope_rows))) {
  r <- slope_rows[i, ]
  cat(sprintf("  %s: est=%.4f, SE=%.4f, z=%.3f, p=%.4f\n", r$rhs, r$est, r$se, r$z, r$pvalue), file=con)
}
cat("\nIntercept factor (baseline level) regressed on:\n", file=con)
for (i in seq_len(nrow(intercept_rows))) {
  r <- intercept_rows[i, ]
  cat(sprintf("  %s: est=%.4f, SE=%.4f, z=%.3f, p=%.4f\n", r$rhs, r$est, r$se, r$z, r$pvalue), file=con)
}
msg("\nDone. See %s/%s_summary.txt\n", opt$outdir, opt$roi)
