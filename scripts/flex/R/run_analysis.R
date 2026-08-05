#!/usr/bin/env Rscript
#
# Analysis-stage entry point for the flex pipeline: sources the core engine
# + every available test_*.R module, then runs cfg$analyses (optionally
# filtered by --only) against the tidy tables the extract stage wrote.
# Invoked by scripts/flex/run_study.py's "analyse" stage; also runnable
# directly for one-off reruns of a subset of the battery.
#
# Usage:
#   Rscript scripts/flex/R/run_analysis.R --config configs/flex/study.pk01.json
#   Rscript scripts/flex/R/run_analysis.R --config ... --only lmm,moderator
#   Rscript scripts/flex/R/run_analysis.R --config ... --only lmm --measures hippo,cortex

suppressPackageStartupMessages(library(optparse))
`%||%` <- function(a, b) if (is.null(a)) b else a

this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file_arg)) return(normalizePath(file_arg))
  normalizePath("scripts/flex/R/run_analysis.R")
}
R_DIR <- dirname(this_file())

option_list <- list(
  make_option("--config", type = "character", help = "Path to a flex study.json"),
  make_option("--only", type = "character", default = NULL,
              help = "Comma-separated test names to run [default: every analysis entry in the config]"),
  make_option("--measures", type = "character", default = NULL,
              help = "Comma-separated measure names to restrict to (within the selected analyses) [default: all]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config)) stop("--config is required")

for (f in c("core_config.R", "core_data.R", "core_design.R", "core_model.R", "core_report.R")) {
  source(file.path(R_DIR, f))
}
test_files <- list.files(R_DIR, pattern = "^test_.*\\.R$", full.names = TRUE)
for (f in test_files) source(f)

cfg <- flex_load_config(opt$config)

only <- if (!is.null(opt$only)) trimws(strsplit(opt$only, ",")[[1]]) else NULL
measures_filter <- if (!is.null(opt$measures)) trimws(strsplit(opt$measures, ",")[[1]]) else NULL

analyses <- cfg$analyses
if (!is.null(only)) {
  analyses <- Filter(function(a) a$test %in% only, analyses)
}
if (!length(analyses)) {
  stop(sprintf("no analyses selected (--only=%s against config test types: %s)",
               paste(only, collapse = ","),
               paste(unique(vapply(cfg$analyses, function(a) a$test, character(1))), collapse = ",")))
}

n_run <- 0L
n_skipped <- 0L

for (spec in analyses) {
  test_name <- spec$test
  fn_name <- paste0("flex_run_", test_name)
  if (!exists(fn_name, mode = "function")) {
    message(sprintf(paste("SKIP: test '%s' has no module yet (expected function %s() in scripts/flex/R/test_%s.R) --",
                           "use the frozen reference script under scripts/analysis/ for this test in the meantime"),
                     test_name, fn_name, test_name))
    n_skipped <- n_skipped + 1L
    next
  }
  fn <- get(fn_name, mode = "function")

  measures <- unlist(spec$measures)
  if (!is.null(measures_filter)) measures <- intersect(measures, measures_filter)

  spec_id <- spec$id %||% test_name
  for (measure_name in measures) {
    outdir <- flex_output_dir(cfg, paste0(spec_id, "_", measure_name))
    message(sprintf("=== %s / %s -> %s ===", test_name, measure_name, outdir))
    result <- tryCatch({
      fn(cfg, measure_name, spec, outdir)
      TRUE
    }, error = function(e) {
      message(sprintf("ERROR in %s/%s: %s", test_name, measure_name, conditionMessage(e)))
      FALSE
    })
    if (isTRUE(result)) n_run <- n_run + 1L else n_skipped <- n_skipped + 1L
  }
}

message(sprintf("\nDone: %d test/measure run(s) completed, %d skipped/failed.", n_run, n_skipped))
if (n_run == 0L) quit(status = 1L)
