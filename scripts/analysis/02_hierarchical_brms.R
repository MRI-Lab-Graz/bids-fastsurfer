#!/usr/bin/env Rscript
#
# Secondary EXPLORATORY analysis: full-subfield hierarchical (partial-pooling)
# model, as a principled alternative to the pre-specified ROI set + FDR
# correction used in 01_primary_lmm.R.
#
# Rather than picking a handful of subfields a priori (or testing all of them
# independently and correcting for multiplicity post hoc -- the approach that
# produced 54 near-duplicate scripts and post-selection p-values in the
# sibling 129run project), this fits ONE model across every subfield with
# subfield-varying dance:time slopes. Noisy subfield estimates shrink toward
# the group mean instead of being tested/corrected individually.
#
#   log(volume) ~ dance * time + hemisphere + age_z + sex + etiv_z
#                 + (1 + dance:time | subfield)
#                 + (1 | subject_id) + (1 | subject_id:subfield)
#
# This is exploratory / hypothesis-generating and reported as a complement to
# the confirmatory Step 3 (01_primary_lmm.R), not a replacement for it.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/02_hierarchical_brms",
              help="Output directory [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--include-composite"), action="store_true", default=FALSE,
              help="Include composite/summary measures (Whole_hippocampus etc.); excluded by default to avoid double-counting with their component subfields"),
  make_option(c("--chains"), type="integer", default=4, help="Number of MCMC chains [default %default]"),
  make_option(c("--iter"), type="integer", default=4000, help="Iterations per chain [default %default]"),
  make_option(c("--cores"), type="integer", default=4, help="Parallel cores for chains [default %default]"),
  make_option(c("--seed"), type="integer", default=134, help="Random seed [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))
if (!requireNamespace("brms", quietly=TRUE)) {
  stop("Package 'brms' is required for this script. Install it (and a working Stan backend,\n",
       "e.g. install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')));\n",
       "cmdstanr::install_cmdstan()) before running.")
}
suppressPackageStartupMessages(library(brms))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

# ------------------------------------------------------------------------
# Load and prepare data (full subfield set, head/body summed per subfield
# as in 01_primary_lmm.R, unless --include-composite is also requested)
# ------------------------------------------------------------------------
msg("Loading tidy subfield data from %s...\n", opt$tidy)
tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$is_composite <- as.logical(tidy$is_composite)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

model_dat_raw <- if (isTRUE(opt$`include-composite`)) tidy else tidy[!tidy$is_composite, , drop=FALSE]
# "other" region (parasubiculum/fimbria/HATA/hippocampal-fissure) has no
# head/body split and is anatomically distinct from CA/DG/subiculum; the
# fissure in particular is CSF, not parenchyma -- exclude fissure specifically.
model_dat_raw <- model_dat_raw[model_dat_raw$subfield_raw != "hippocampal-fissure", , drop=FALSE]

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=model_dat_raw, FUN=sum)
etiv_lookup <- unique(tidy[, c("subject_id","etiv")])
etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), ]
agg <- merge(agg, etiv_lookup, by="subject_id", all.x=TRUE)

msg("Loading participant metadata from %s...\n", opt$participants)
participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(agg, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(dance_groups, control_group), , drop=FALSE]

agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$dance <- factor(ifelse(agg$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$sex <- factor(agg$sex)
agg$age_z <- as.numeric(scale(agg$age))
agg$etiv_z <- as.numeric(scale(agg$etiv))
agg$log_volume <- log(agg$volume)
agg$subfield <- factor(agg$subfield)
agg$subject_id <- factor(agg$subject_id)

n_subfields <- nlevels(agg$subfield)
n_subjects <- length(unique(agg$subject_id))
msg("Modelling data: %d rows, %d subjects, %d subfields\n", nrow(agg), n_subjects, n_subfields)

# ------------------------------------------------------------------------
# Hierarchical model with subfield-varying dance:time slopes
# ------------------------------------------------------------------------
priors <- c(
  brms::prior(normal(0, 1), class="b"),
  brms::prior(normal(0, 1), class="sd"),
  brms::prior(normal(0, 1), class="sigma")
)

msg("Fitting hierarchical brms model (this can take a while)...\n")
fit <- brms::brm(
  log_volume ~ dance * time_f + hemisphere + age_z + sex + etiv_z +
    (1 + dance:time_f | subfield) + (1 | subject_id) + (1 | subject_id:subfield),
  data = agg,
  prior = priors,
  chains = opt$chains,
  iter = opt$iter,
  cores = opt$cores,
  seed = opt$seed,
  control = list(adapt_delta = 0.95)
)

saveRDS(fit, file.path(opt$outdir, "hierarchical_model.rds"))

# ------------------------------------------------------------------------
# Per-subfield posterior summaries for the dance effect at each timepoint
# ------------------------------------------------------------------------
# "subfield" is a group-level (random-effect) grouping factor, not a
# fixed-effect term, so its levels aren't part of brms::coef()'s
# population-level naming scheme (which fully dummy-codes group-level
# interactions and does not map cleanly back to a single "dance:time"
# coefficient per subfield). Instead, compute the dance-vs-control
# difference directly from the posterior predictive distribution: predict
# at both dance levels x all timepoints x each subfield, holding covariates
# at reference values, with re_formula restricted to the subfield-level
# grouping only (so predictions are marginalized over subject_id and
# subject_id:subfield, i.e. a population-average subfield effect).
msg("Extracting per-subfield posterior summaries...\n")

subfields <- levels(agg$subfield)
times <- levels(agg$time_f)
ref_sex <- names(sort(table(agg$sex), decreasing=TRUE))[1]

newdata <- expand.grid(
  subfield = subfields,
  dance = c("control", "dance"),
  time_f = times,
  hemisphere = "lh",
  age_z = 0,
  etiv_z = 0,
  sex = ref_sex,
  stringsAsFactors = FALSE
)
newdata$dance <- factor(newdata$dance, levels=levels(agg$dance))
newdata$time_f <- factor(newdata$time_f, levels=levels(agg$time_f))
newdata$hemisphere <- factor(newdata$hemisphere, levels=levels(agg$hemisphere))
newdata$sex <- factor(newdata$sex, levels=levels(agg$sex))
newdata$subfield <- factor(newdata$subfield, levels=levels(agg$subfield))

epred <- brms::posterior_epred(
  fit, newdata = newdata,
  re_formula = ~ (1 + dance:time_f | subfield),
  allow_new_levels = FALSE
)
# epred: draws x rows(newdata). For each subfield/time, compute the
# per-draw (dance - control) difference on the log-volume scale, then
# summarize across draws.
posterior_rows <- list()
for (sf in subfields) {
  for (tm in times) {
    dance_idx <- which(newdata$subfield == sf & newdata$time_f == tm & newdata$dance == "dance")
    control_idx <- which(newdata$subfield == sf & newdata$time_f == tm & newdata$dance == "control")
    diff_draws <- epred[, dance_idx] - epred[, control_idx]
    posterior_rows[[paste(sf, tm)]] <- data.frame(
      subfield = sf,
      time_f = tm,
      estimate_log_diff = mean(diff_draws),
      q2.5 = as.numeric(quantile(diff_draws, 0.025)),
      q97.5 = as.numeric(quantile(diff_draws, 0.975)),
      prob_dance_gt_control = mean(diff_draws > 0),
      stringsAsFactors = FALSE
    )
  }
}
posterior_summary <- do.call(rbind, posterior_rows)
write.csv(posterior_summary, file.path(opt$outdir, "subfield_dance_time_posteriors.csv"), row.names=FALSE)

msg("Model diagnostics (Rhat, ESS)...\n")
diag_df <- as.data.frame(brms::rhat(fit))
write.csv(diag_df, file.path(opt$outdir, "rhat_diagnostics.csv"), row.names=TRUE)
max_rhat <- suppressWarnings(max(brms::rhat(fit), na.rm=TRUE))
if (is.finite(max_rhat) && max_rhat > 1.01) {
  warning(sprintf("Max Rhat = %.3f (> 1.01) -- consider more iterations/chains before trusting these posteriors", max_rhat))
}

msg("Done. Model saved to hierarchical_model.rds; per-subfield posteriors in subfield_dance_time_posteriors.csv\n")
