#!/usr/bin/env Rscript
#
# Direct replication check of a colleague's reported finding on study 129:
# a significant HEMISPHERE x TIME x GROUP interaction in the hippocampal
# tail, F(2,115)=3.70, p=.028, np2=.060 (3-level GROUP: solo / small-group
# / control -- NOT our usual pooled dance-vs-control binary contrast).
#
# An earlier check on our FULL sample found this null (np2=0.0000279,
# p=0.408). The user hypothesized the discrepancy might be explained by
# the colleague using a smaller, differently-selected sample (excluding
# high-baseline-aerobic-activity participants) -- a hypothesis we could not
# verify at the time because we didn't have that selection.
#
# We now have it: a jamovi .omv export (data/JAMOVI_final.omv) contains
# exactly 127 subjects with group sizes 42/43/42 and matching demographics
# to the colleague's reported solo/group/control ns -- i.e. very likely
# their actual analysis sample. This script re-runs the SAME rm-ANOVA
# design (afex::aov_ez, within = time_f x hemisphere, between = 3-level
# group) restricted to exactly these 127 subjects, to see whether the
# finding replicates once the sample matches.
#
# LOSO validation follows immediately given how much rides on this one
# result -- this project's standing rule that no hit gets reported without
# checking single-subject fragility applies here more than anywhere else.

suppressPackageStartupMessages({
  library(optparse)
  library(afex)
})

option_list <- list(
  make_option(c("-t", "--tidy"), type="character",
              help="Path to hippo_subfields_tidy.tsv"),
  make_option(c("-p", "--participants"), type="character",
              help="Path to participants TSV restricted to the replication subsample (subject_id, group, age, sex)"),
  make_option(c("-o", "--outdir"), type="character", default="results/22_hippotail_colleague_replication"),
  make_option(c("--roi-set"), type="character",
              default="Whole_hippocampus,Whole_hippocampal_head,Whole_hippocampal_body,Hippocampal_tail"),
  make_option(c("--solo-pattern"), type="character", default="Single"),
  make_option(c("--group-pattern"), type="character", default="Group"),
  make_option(c("--control-value"), type="character", default="Control"),
  make_option(c("--alpha"), type="double", default=0.05),
  make_option(c("--quiet"), action="store_true", default=FALSE)
)
opt <- parse_args(OptionParser(option_list=option_list))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

if (is.null(opt$tidy)) stop("--tidy is required")
if (is.null(opt$participants)) stop("--participants is required")

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

roi_set <- trimws(strsplit(opt$`roi-set`, ",")[[1]])

tidy <- read.delim(opt$tidy, header=TRUE, sep="\t", stringsAsFactors=FALSE)
tidy$volume <- suppressWarnings(as.numeric(tidy$volume))
roi_dat <- tidy[tidy$subfield %in% roi_set, , drop=FALSE]
if (!nrow(roi_dat)) stop("No rows matched --roi-set")

agg <- aggregate(volume ~ subject_id + session + hemisphere + subfield, data=roi_dat, FUN=sum)

participants <- read.delim(opt$participants, header=TRUE, sep="\t", stringsAsFactors=FALSE)
required_pcols <- c("subject_id","group","age","sex")
missing_pcols <- setdiff(required_pcols, names(participants))
if (length(missing_pcols)) stop(sprintf("participants file missing required columns: %s", paste(missing_pcols, collapse=", ")))

participants$group3 <- ifelse(grepl(opt$`solo-pattern`, participants$group), "Solo",
                        ifelse(grepl(opt$`group-pattern`, participants$group), "Group",
                        ifelse(participants$group == opt$`control-value`, "Control", NA_character_)))
if (any(is.na(participants$group3))) stop("Some subjects in --participants didn't classify into Solo/Group/Control -- check --solo-pattern/--group-pattern/--control-value")

msg("3-level group sizes in this replication sample:\n")
print(table(participants$group3))

agg <- merge(agg, participants[, c("subject_id","group3","age","sex")], by="subject_id")
agg$time_f <- factor(agg$session, levels=sort(unique(agg$session)))
agg$hemisphere <- factor(agg$hemisphere, levels=c("lh","rh"))
agg$group3 <- factor(agg$group3, levels=c("Control","Solo","Group"))
agg$subject_id <- factor(agg$subject_id)

fit_roi <- function(d_roi) {
  tryCatch(
    afex::aov_ez(
      id = "subject_id", dv = "volume", data = d_roi,
      within = c("time_f", "hemisphere"), between = "group3",
      fun_aggregate = mean, anova_table = list(es = "pes")
    ),
    error = function(e) { message(sprintf("aov_ez failed: %s", e$message)); NULL }
  )
}

get_3way <- function(fit) {
  at <- as.data.frame(fit$anova_table)
  at$term <- rownames(at)
  row <- at[grepl("group3.*time_f.*hemisphere|time_f.*hemisphere.*group3|hemisphere.*time_f.*group3", at$term), ]
  row
}

results_list <- list()
fit_by_roi <- list()
for (roi_name in roi_set) {
  msg("\nFitting HEM x TIME x GROUP(3) rm-ANOVA for %s...\n", roi_name)
  d_roi <- agg[agg$subfield == roi_name, , drop=FALSE]
  n_before <- length(unique(d_roi$subject_id))
  fit <- fit_roi(d_roi)
  if (is.null(fit)) next
  fit_by_roi[[roi_name]] <- d_roi

  anova_table <- as.data.frame(fit$anova_table)
  anova_table$term <- rownames(anova_table)
  write.csv(anova_table, file.path(opt$outdir, paste0(roi_name, "_anova_table.csv")), row.names=FALSE)

  row3 <- get_3way(fit)
  n_complete <- nrow(fit$data$long) / (nlevels(d_roi$time_f) * nlevels(d_roi$hemisphere))
  msg("  n_subjects: %d available, %d complete cases\n", n_before, n_complete)
  if (nrow(row3)) {
    msg("  HEM x TIME x GROUP: F(%s)=%.3f, p=%.4f, np2=%.4f\n",
        row3$term[1], row3$F[1], row3$`Pr(>F)`[1], row3$pes[1])
  }

  results_list[[roi_name]] <- data.frame(
    roi = roi_name, n_available = n_before, n_complete = n_complete,
    F_value = if (nrow(row3)) row3$F[1] else NA_real_,
    df = if (nrow(row3)) row3$term[1] else NA_character_,
    p_value = if (nrow(row3)) row3$`Pr(>F)`[1] else NA_real_,
    np2 = if (nrow(row3)) row3$pes[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, results_list)
if (is.null(summary_df)) stop("No ROI models were successfully fit")
write.csv(summary_df, file.path(opt$outdir, "replication_summary.csv"), row.names=FALSE)
msg("\nWrote replication_summary.csv\n")

# ------------------------------------------------------------------------
# LOSO for the hippocampal tail's 3-way term specifically (and any other
# ROI landing below alpha) -- given how much rides on this single test
# ------------------------------------------------------------------------
screen <- summary_df[!is.na(summary_df$p_value) & summary_df$p_value < opt$alpha, ]
if (nrow(screen)) {
  msg("\nRunning LOSO for %d ROI(s) with p < %.2f...\n", nrow(screen), opt$alpha)
  con <- file(file.path(opt$outdir, "loso_summary.txt"), open="wt")
  on.exit(close(con), add=TRUE)
  cat("LOSO for the HEM x TIME x GROUP replication test\n\n", file=con)
  for (roi_name in screen$roi) {
    d_roi <- fit_by_roi[[roi_name]]
    full_p <- screen$p_value[screen$roi == roi_name]
    subjects <- unique(d_roi$subject_id)
    loso_p <- sapply(subjects, function(s) {
      d_sub <- d_roi[d_roi$subject_id != s, , drop=FALSE]
      d_sub$subject_id <- factor(d_sub$subject_id)
      if (length(unique(d_sub$group3)) < 3) return(NA_real_)
      fit <- fit_roi(d_sub)
      if (is.null(fit)) return(NA_real_)
      row3 <- get_3way(fit)
      if (!nrow(row3)) return(NA_real_)
      row3$`Pr(>F)`[1]
    })
    loso_p <- loso_p[!is.na(loso_p)]
    n_flips <- sum((loso_p < opt$alpha) != (full_p < opt$alpha))
    cat(sprintf("%s: full p=%.4f, LOSO range [%.4f, %.4f], %d/%d exclusions flip significance\n",
                roi_name, full_p, min(loso_p), max(loso_p), n_flips, length(loso_p)), file=con)
    write.csv(data.frame(excluded_subject=subjects[!is.na(sapply(subjects, function(s) 1))][seq_along(loso_p)], loso_p=loso_p),
              file.path(opt$outdir, paste0(roi_name, "_loso.csv")), row.names=FALSE)
  }
  msg("Done. See %s/loso_summary.txt\n", opt$outdir)
} else {
  msg("\nNo ROI reached p < %.2f; nothing to LOSO-check.\n", opt$alpha)
}
