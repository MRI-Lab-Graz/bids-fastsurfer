# core_data.R -- load tidy measure tables + participants/moderators, ROI
# selection, and the two data-prep steps every profile-based test needs:
# head/body ROI aggregation and the baseline-session eTIV lookup.

flex_load_tidy <- function(cfg, measure_name) {
  path <- flex_tidy_path(cfg, measure_name)
  if (!file.exists(path)) {
    stop(sprintf("tidy file not found for measure '%s': %s -- run the extract stage first (scripts/flex/extract_freesurfer.py --config ...)",
                  measure_name, path))
  }
  tidy <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  tidy$is_composite <- as.logical(tidy$is_composite)
  tidy$value <- suppressWarnings(as.numeric(tidy$value))
  required <- c("subject_id", "session", "hemisphere", "region", "value")
  missing_cols <- setdiff(required, names(tidy))
  if (length(missing_cols)) {
    stop(sprintf("tidy file for measure '%s' missing required column(s): %s",
                  measure_name, paste(missing_cols, collapse = ", ")))
  }
  tidy
}

# Participants file, with subject/group role columns aliased to canonical
# names (subject_id, group) so every downstream module can join on a fixed
# name regardless of what the study's own file calls them -- the original
# columns (incl. age/sex/whatever design$covariates references by name)
# are kept as-is alongside the aliases.
flex_load_participants <- function(cfg) {
  path <- flex_resolve_path(cfg$design$participants)
  if (!file.exists(path)) stop(sprintf("participants file not found: %s", path))
  participants <- read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  subj_col <- cfg$design$roles$subject
  group_col <- cfg$design$roles$group
  if (!subj_col %in% names(participants)) stop(sprintf("design.roles.subject '%s' not found in participants file", subj_col))
  if (!group_col %in% names(participants)) stop(sprintf("design.roles.group '%s' not found in participants file", group_col))

  participants$subject_id <- participants[[subj_col]]
  participants$group <- participants[[group_col]]
  participants
}

flex_load_moderators <- function(cfg) {
  mod_cfg <- cfg$design$moderators
  if (is.null(mod_cfg)) return(NULL)
  path <- flex_resolve_path(mod_cfg$file)
  if (!file.exists(path)) stop(sprintf("moderators file not found: %s", path))
  read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

# roi_set: character vector of region names, or "all" (keep everything).
flex_select_roi <- function(tidy, roi_set, roi_column = "region") {
  if (is.null(roi_set) || identical(roi_set, "all") || (length(roi_set) == 1 && roi_set == "all")) {
    return(tidy)
  }
  roi_set <- unlist(roi_set)
  out <- tidy[tidy[[roi_column]] %in% roi_set, , drop = FALSE]
  if (!nrow(out)) {
    stop(sprintf("no rows matched roi_set (%s); check region names against the tidy file's '%s' column",
                  paste(roi_set, collapse = ", "), roi_column))
  }
  out
}

# Head/body ROI aggregation: sums `value` within subject_id/session/
# hemisphere/roi_column. Identical in spirit to 01_primary_lmm.R's
# `aggregate(volume ~ subject_id + session + hemisphere + subfield, ...,
# FUN=sum)` -- subfield identity is already head/body-collapsed by
# extract_freesurfer.py's `region` column (hippo_base_name()), so rows
# differing only in the stripped head/body suffix land in the same group
# here and get summed into one subfield-total. For measures with no
# head/body split (aseg/thalamus/brainstem/cortex), every group already has
# exactly one row, so the sum is a no-op -- safe to apply unconditionally.
flex_aggregate_roi <- function(tidy, roi_column = "region") {
  f <- as.formula(paste("value ~ subject_id + session + hemisphere +", roi_column))
  stats::aggregate(f, data = tidy, FUN = sum)
}

# eTIV varies slightly session-to-session (FreeSurfer re-estimation noise,
# not real anatomical change) -- use each subject's BASELINE session eTIV as
# a fixed per-subject covariate, explicitly selected (not relying on
# incidental row order), exactly mirroring 01_primary_lmm.R lines 110-112.
# baseline_session: explicit session value, or NULL to use the earliest
# session present in the tidy table (alphabetical/lexicographic min, same
# convention as 01_primary_lmm.R's sort(unique(...))[1]).
flex_baseline_etiv <- function(tidy, baseline_session = NULL) {
  if (is.null(baseline_session)) {
    baseline_session <- sort(unique(tidy$session))[1]
  }
  etiv_lookup <- unique(tidy[tidy$session == baseline_session, c("subject_id", "etiv")])
  etiv_lookup <- etiv_lookup[!duplicated(etiv_lookup$subject_id), , drop = FALSE]
  etiv_lookup
}

flex_prepare_measure <- function(cfg, measure_name) {
  measure <- flex_get_measure(cfg, measure_name)
  profile <- flex_resolve_profile(cfg, measure$profile)
  tidy <- flex_load_tidy(cfg, measure_name)

  sessions_cfg <- cfg$sessions
  if (!is.null(sessions_cfg) && !is.null(sessions_cfg$include)) {
    include <- unlist(sessions_cfg$include)
    tidy <- tidy[tidy$session %in% include, , drop = FALSE]
  }

  roi_column <- measure$roi_column %||% "region"
  tidy <- flex_select_roi(tidy, measure$roi_set, roi_column)

  list(measure = measure, profile = profile, tidy = tidy, roi_column = roi_column)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
