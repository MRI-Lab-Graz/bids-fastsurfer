#!/usr/bin/env Rscript
#
# Reviewer-facing repeated-measures ANOVA on whole hippocampus and the
# head/body/tail regions. This is a SECONDARY analysis reported alongside
# the primary mixed model (01_primary_lmm.R) because reviewers often expect
# it -- it is NOT the primary analysis, since rm-ANOVA requires complete
# cases and drops any subject missing a timepoint (unlike the mixed model,
# which uses all available data).

suppressPackageStartupMessages({
  library(optparse)
  library(afex)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv (from extract_hippo_subfields.py)"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV with columns: subject_id, group, age, sex"),
  make_option(c("-o", "--outdir"), type="character", default="results/03_rm_anova",
              help="Output directory [default %default]"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,Whole_hippocampal_head,Whole_hippocampal_body,Hippocampal_tail",
              help="Comma-separated regions to test [default %default]"),
  make_option(c("--dance-groups"), type="character", default="ballet,contemporary",
              help="Comma-separated group values pooled into the 'dance' contrast [default %default]"),
  make_option(c("--control-group"), type="character", default="control",
              help="Group value treated as control [default %default]"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")
if (!file.exists(opt$tidy)) stop(sprintf("tidy file not found: %s", opt$tidy))
if (!file.exists(opt$participants)) stop(sprintf("participants file not found: %s", opt$participants))

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])
dance_groups <- trimws(strsplit(opt$`dance-groups`, ",")[[1]])
control_group <- trimws(opt$`control-group`)

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))

roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set; check region names against the tidy file (Whole_hippocampus / Whole_hippocampal_head / Whole_hippocampal_body / Hippocampal_tail)")

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=roi_dat, FUN=sum)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

agg <- merge(agg, participants[, required_pcols], by="subject_id")
agg <- agg[agg$group %in% c(dance_groups, control_group), , drop=FALSE]
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$dance <- factor(ifelse(agg$group %in% dance_groups, "dance", "control"), levels=c("control","dance"))
agg$subject_id <- factor(agg$subject_id)

results_list <- list()
for (roi_name in roi_set) {
  msg("Running rm-ANOVA for region: %s\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  if (!nrow(d_roi)) { msg("  -> no data, skipping\n"); next }

  n_before <- length(unique(d_roi$subject_id))
  fit <- tryCatch(
    afex::aov_ez(
      id = "subject_id",
      dv = "volume",
      data = d_roi,
      within = c("time_f", "hemisphere"),
      between = "dance",
      fun_aggregate = mean # collapse any accidental duplicate cells
    ),
    error = function(e) { message(sprintf("aov_ez failed for %s: %s", roi_name, e$message)); NULL }
  )
  if (is.null(fit)) next

  n_complete <- nrow(fit$data$long) / (nlevels(d_roi$time_f) * nlevels(d_roi$hemisphere))
  if (n_complete < n_before) {
    msg("  -> %d of %d subjects retained after listwise deletion for missing timepoints\n", n_complete, n_before)
  }

  anova_table <- as.data.frame(fit$anova_table)
  anova_table$term <- rownames(anova_table)
  write.csv(anova_table, file.path(opt$outdir, paste0(roi_name, "_anova_table.csv")), row.names=FALSE)
  saveRDS(fit, file.path(opt$outdir, paste0(roi_name, "_aov_fit.rds")))

  results_list[[roi_name]] <- data.frame(
    roi = roi_name,
    n_subjects_complete_cases = n_complete,
    n_subjects_available = n_before,
    dance_time_p = anova_table[grepl("dance.*time_f|time_f.*dance", anova_table$term), "Pr(>F)"][1],
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, results_list)
if (!is.null(summary_df)) {
  write.csv(summary_df, file.path(opt$outdir, "rm_anova_summary.csv"), row.names=FALSE)
  msg("Wrote summary to rm_anova_summary.csv\n")
} else {
  warning("No ANOVA models were successfully fit")
}
