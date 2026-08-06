# core_config.R -- load a flex study.json for the R analysis engine.
#
# Primary schema validation happens in Python (scripts/flex/config.py,
# used by ingest.py/extract_freesurfer.py/run_study.py) so there is exactly
# one definition of "valid" shared by every stage; run_study.py validates
# before invoking any Rscript. This file does light defensive checks only
# (a config could in principle be handed to an R module directly), and is
# the single place every test_*.R module gets a config, a profile, and a
# resolved path from -- no module parses JSON or expands a profile itself.

suppressPackageStartupMessages(library(jsonlite))

# Built-in measure profiles -- MUST mirror scripts/flex/config.py's
# BUILTIN_PROFILES exactly (Python is authoritative for what's a *valid*
# profile at config-validation time; this is the modelling-time mirror).
# See docs/FLEX_PIPELINE.md for the rationale: "volume" (log-transform +
# eTIV covariate) carries forward 01_primary_lmm.R's header rationale,
# "thickness" (raw scale, no eTIV) carries forward 29_hippo_thickness_lmm.R's.
flex_builtin_profiles <- list(
  volume          = list(transform = "log",      etiv_covariate = TRUE,  hemisphere = "pooled", aggregate = "sum"),
  thickness       = list(transform = "identity", etiv_covariate = FALSE, hemisphere = "pooled", aggregate = "mean"),
  volume_midline  = list(transform = "log",      etiv_covariate = TRUE,  hemisphere = "none",   aggregate = "sum")
)

flex_load_config <- function(path) {
  if (!file.exists(path)) stop(sprintf("flex study config not found: %s", path))
  cfg <- jsonlite::fromJSON(path, simplifyVector = FALSE)

  required_top <- c("study", "datasets", "measures", "design", "analyses", "output")
  missing_top <- setdiff(required_top, names(cfg))
  if (length(missing_top)) {
    stop(sprintf("config missing required top-level field(s): %s -- was this validated by run_study.py / scripts/flex/config.py first?",
                  paste(missing_top, collapse = ", ")))
  }
  if (is.null(cfg$design$roles$subject)) stop("config missing design.roles.subject")
  if (is.null(cfg$design$roles$group)) stop("config missing design.roles.group")

  cfg$.config_path <- normalizePath(path)
  cfg$.config_dir <- dirname(cfg$.config_path)
  cfg
}

# Relative paths in a config (participants, moderators.file, region_groups,
# output.root) are resolved against the repo root -- Rscript invocations in
# this repo are always run from the repo root (see e.g.
# tmp/pipeline_2tp_fs82/run_pipeline.sh's `cd` at the top), matching every
# other config in configs/.
flex_resolve_path <- function(p) {
  if (grepl("^(/|~|[A-Za-z]:)", p)) return(p)
  file.path(getwd(), p)
}

flex_resolve_profile <- function(cfg, profile_name) {
  custom <- cfg$profiles
  if (!is.null(custom) && !is.null(custom[[profile_name]])) {
    prof <- custom[[profile_name]]
  } else if (!is.null(flex_builtin_profiles[[profile_name]])) {
    prof <- flex_builtin_profiles[[profile_name]]
  } else {
    stop(sprintf("unknown profile '%s' (known: %s)", profile_name,
                  paste(union(names(flex_builtin_profiles), names(custom)), collapse = ", ")))
  }
  prof
}

flex_get_measure <- function(cfg, measure_name) {
  names_avail <- vapply(cfg$measures, function(m) m$name, character(1))
  idx <- match(measure_name, names_avail)
  if (is.na(idx)) stop(sprintf("unknown measure '%s' (known: %s)", measure_name, paste(names_avail, collapse = ", ")))
  cfg$measures[[idx]]
}

flex_output_dir <- function(cfg, ...) {
  root <- flex_resolve_path(cfg$output$root)
  d <- file.path(root, "results", ...)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

flex_tidy_path <- function(cfg, measure_name) {
  file.path(flex_resolve_path(cfg$output$root), "tidy", paste0(measure_name, ".tsv"))
}

flex_alpha <- function(cfg) {
  a <- cfg$output$alpha
  if (is.null(a)) 0.05 else a
}
