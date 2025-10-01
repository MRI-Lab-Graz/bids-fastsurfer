#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  # GLM base
  library(stats)
  # GAM (optional); load lazily later to avoid hard dependency
})

# Flexible univariate LME analysis for aseg/aparc tables using fslmer

# Create comprehensive help text with examples
get_help_text <- function() {
  help_text <- "
DESCRIPTION:
    Flexible univariate linear mixed-effects (LME) analysis for subcortical/cortical 
    volume data from FreeSurfer/FastSurfer longitudinal processing. Supports single ROI 
    or multi-region analysis with various statistical models (fslmer, GLM, GAM).

USAGE:
    # Single ROI analysis
    Rscript fslmer_univariate.R --qdec qdec.table.dat --aseg aseg.long.table \\
                                --roi Left.Hippocampus --formula '~ tp*group' \\
                                --outdir results_hippo

    # Multiple regions matching pattern
    Rscript fslmer_univariate.R --qdec qdec.table.dat --aseg aseg.long.table \\
                                --region-pattern 'Hippocampus|Amygdala' \\
                                --formula '~ tp*group' --outdir results_limbic

    # All brain regions analysis
    Rscript fslmer_univariate.R --qdec qdec.table.dat --aseg aseg.long.table \\
                                --all-regions --formula '~ tp' \\
                                --outdir results_all

    # Using configuration file
    Rscript fslmer_univariate.R --config configs/fslmer_univariate.example.json

    # Via bash wrapper (recommended)
    bash scripts/run_fslmer_univariate.sh --help

EXAMPLES:
    # Basic longitudinal change
    --formula '~ tp'
    
    # Group by time interaction
    --formula '~ tp*group'
    
    # With baseline covariate
    --add-baseline --formula '~ baseline_value + tp*group'
    
    # Multiple time points with quadratic
    --formula '~ tp + I(tp^2)'

    # GAM with smooth time effect
    --engine gam --formula '~ s(tp)'

INPUT FILES:
    qdec.table.dat:     Subject metadata (TSV) with fsid, fsid_base, tp, group, etc.
    aseg.long.table:    FreeSurfer longitudinal volume table with .long. identifiers

OUTPUT:
    {outdir}/results_summary.csv:     Summary statistics for all ROIs
    {outdir}/models/{roi}_model.txt:  Individual model results per ROI
    {outdir}/merged_data.csv:         Combined input data (if --save-merged)
"
  return(help_text)
}

# Create argument groups for better organization
option_list <- list(
  # Required inputs
  make_option(c("--config"), type="character", default=NULL, 
             help="Path to JSON config file with arguments"),
  make_option(c("-q", "--qdec"), type="character", 
             help="Path to qdec.table.dat (TSV with subject metadata)"),
  make_option(c("-a", "--aseg"), type="character", 
             help="Path to aseg/aparc longitudinal table (TSV/whitespace format)"),
  
  # ROI selection (choose one approach)
  make_option(c("-r", "--roi"), type="character", 
             help="Single ROI column name (e.g., 'Left.Hippocampus')"),
  make_option(c("--all-regions"), action="store_true", default=FALSE, 
             help="Analyze all brain regions (excludes summary measures by default)"),
  make_option(c("--region-pattern"), type="character", default=NULL, 
             help="Regex pattern to match ROI names (e.g., 'Hippocampus|Amygdala')"),
  make_option(c("--include-summary"), action="store_true", default=FALSE, 
             help="Include global/summary volume measures in multi-region analysis"),
  
  # Statistical model configuration
  make_option(c("-f", "--formula"), type="character", default="~ tp", 
             help="Model formula for fixed effects (e.g., '~ tp*group') [default: %default]"),
  make_option(c("--random-effects"), type="character", default=NULL, 
             help="Random-effects structure: RI (random intercept), RS (random slope), RIRS (both), or indices '1,2'"),
  make_option(c("-z", "--zcols"), type="character", default="1,2", 
             help="[Deprecated] Use --random-effects instead. Random effects indices [default: %default]"),
  make_option(c("-c", "--contrast"), type="character", default=NULL, 
             help="Comma-separated contrast vector (length must match model coefficients)"),
  make_option(c("--engine"), type="character", default="fslmer", 
             help="Statistical engine: 'fslmer' (default), 'glm', or 'gam'"),
  make_option(c("--family"), type="character", default="gaussian", 
             help="GLM/GAM family: gaussian, binomial, poisson, etc. [ignored for fslmer]"),
  
  # Advanced model options
  make_option(c("--add-baseline"), action="store_true", default=FALSE, 
             help="Add per-ROI baseline covariate (value at first visit per subject)"),
  make_option(c("--derive-group5"), action="store_true", default=FALSE, 
             help="Derive intervention covariates from group_5: interv, smallgroup, weeks4, weeks"),
  make_option(c("--gam-k"), type="integer", default=NA, 
             help="GAM basis dimension k for s(time) [default: auto, min(5, unique timepoints)]"),
  make_option(c("--no-gam-re"), action="store_true", default=FALSE, 
             help="Disable random intercept s(fsid_base, bs='re') in GAM models"),
  
  # Input/output options
  make_option(c("-o", "--outdir"), type="character", default="fslmer_out", 
             help="Output directory [default: %default]"),
  make_option(c("--save-merged"), action="store_true", default=FALSE, 
             help="Save merged input data as CSV for inspection"),
  make_option(c("--time-col"), type="character", default=NULL, 
             help="Time variable column name in qdec [default: auto-detect 'tp']"),
  make_option(c("--id-col"), type="character", default=NULL, 
             help="ID column in aseg table [default: auto-detect 'Measure:volume']"),
  
  # Utility options
  make_option(c("--print-cols"), action="store_true", default=FALSE, 
             help="Print column names and preview data then exit"),
  make_option(c("--quiet"), action="store_true", default=FALSE, 
             help="Reduce output verbosity"),
  make_option(c("--debug"), action="store_true", default=FALSE, 
             help="Enable detailed optimizer logs (fslmer only)"),
  
  # Statistical summary and visualization options
  make_option(c("--alpha"), type="double", default=0.05, 
             help="Significance threshold (p < alpha) [default: %default]"),
  make_option(c("--trend"), type="double", default=0.1, 
             help="Trend threshold (alpha <= p < trend) [default: %default]"),
  make_option(c("--plot-results"), action="store_true", default=FALSE, 
             help="Generate plots for significant and trend results"),
  make_option(c("--no-summary"), action="store_true", default=FALSE, 
             help="Skip statistical significance summary in terminal output"),
  
  make_option(c("-h", "--help"), action="store_true", default=FALSE,
             help="Show this help message and exit")
)

# Custom argument parser to handle help
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && (args[1] %in% c("-h", "--help") || any(grepl("^--help", args)))) {
  cat(get_help_text())
  cat("\nARGUMENT DETAILS:\n")
  parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
  print_help(parser)
  quit(status=0)
}

opt <- parse_args(OptionParser(option_list=option_list, add_help_option=FALSE))
msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

# Config JSON merge (CLI overrides)
if (!is.null(opt$config)) {
  if (!file.exists(opt$config)) stop(sprintf("Config file not found: %s", opt$config))
  cfg <- jsonlite::fromJSON(opt$config, simplifyVector = TRUE)

  defaults <- list(
    qdec = NULL,
    aseg = NULL,
    roi = NULL,
    formula = "~ tp",
    zcols = "1,2",
    `random-effects` = NULL,
    contrast = NULL,
    outdir = "fslmer_out",
    time_col = NULL,
    id_col = NULL,
    region_pattern = NULL,
    all_regions = FALSE,
    add_baseline = FALSE,
    derive_group5 = FALSE
  )

  apply_cfg <- function(name) {
    default_val <- defaults[[name]]
    cli_val <- opt[[name]]
    cfg_val <- cfg[[name]]

    if (!is.null(cli_val) && !identical(cli_val, default_val)) {
      return(cli_val)
    }
    if (!is.null(cfg_val)) {
      return(cfg_val)
    }
    if (!is.null(cli_val)) {
      return(cli_val)
    }
    default_val
  }

  opt$qdec     <- apply_cfg("qdec")
  opt$aseg     <- apply_cfg("aseg")
  opt$roi      <- apply_cfg("roi")
  opt$formula  <- apply_cfg("formula")
  opt$zcols    <- apply_cfg("zcols")
  opt$`random-effects` <- apply_cfg("random-effects")
  opt$contrast <- apply_cfg("contrast")
  opt$outdir   <- apply_cfg("outdir")
  opt$time_col <- apply_cfg("time_col")
  opt$id_col   <- apply_cfg("id_col")
  # Region selection from JSON (supports region_pattern or region-pattern, all_regions)
  rp <- cfg[["region_pattern"]]; if (is.null(rp)) rp <- cfg[["region-pattern"]]
  if (!is.null(rp)) opt$`region-pattern` <- rp
  ar <- cfg[["all_regions"]]; if (is.null(ar)) ar <- cfg[["all-regions"]]
  if (!is.null(ar)) opt$`all-regions` <- isTRUE(ar)
  opt$add_baseline <- isTRUE(apply_cfg("add_baseline"))
  opt$derive_group5 <- isTRUE(apply_cfg("derive_group5"))

  # Accept random-effects from JSON using underscore or hyphen key
  re_json <- cfg[["random_effects"]]; if (is.null(re_json)) re_json <- cfg[["random-effects"]]
  if (!is.null(re_json)) opt$`random-effects` <- re_json
}

# Prefer --random-effects over --zcols if provided
if (!is.null(opt$`random-effects`) && nzchar(as.character(opt$`random-effects`))) {
  opt$zcols <- opt$`random-effects`
}

# Comprehensive argument validation
validate_arguments <- function(opt) {
  errors <- character(0)
  
  # Required input files
  if (is.null(opt$qdec)) {
    errors <- c(errors, "Missing required argument: --qdec (path to qdec.table.dat)")
  } else if (!file.exists(opt$qdec)) {
    errors <- c(errors, sprintf("qdec file not found: %s", opt$qdec))
  }
  
  if (is.null(opt$aseg)) {
    errors <- c(errors, "Missing required argument: --aseg (path to aseg/aparc table)")
  } else if (!file.exists(opt$aseg)) {
    errors <- c(errors, sprintf("aseg/aparc table not found: %s", opt$aseg))
  }
  
  # ROI selection validation - exactly one method required
  roi_methods <- sum(c(
    !is.null(opt$roi),
    isTRUE(opt$`all-regions`),
    !is.null(opt$`region-pattern`)
  ))
  
  if (roi_methods == 0) {
    errors <- c(errors, "ROI selection required: specify --roi, --all-regions, or --region-pattern")
  } else if (roi_methods > 1) {
    errors <- c(errors, "Multiple ROI selection methods specified. Choose only one: --roi, --all-regions, or --region-pattern")
  }
  
  # Engine validation
  valid_engines <- c("fslmer", "glm", "gam")
  if (!opt$engine %in% valid_engines) {
    errors <- c(errors, sprintf("Invalid engine '%s'. Valid options: %s", opt$engine, paste(valid_engines, collapse=", ")))
  }
  
  # Formula validation
  if (is.null(opt$formula) || nchar(trimws(opt$formula)) == 0) {
    errors <- c(errors, "Empty or missing formula. Use --formula '~ tp' for basic time effect")
  } else {
    tryCatch(as.formula(opt$formula), error = function(e) {
      errors <<- c(errors, sprintf("Invalid formula syntax: %s", e$message))
    })
  }
  
  # GAM-specific validation
  if (opt$engine == "gam") {
    if (!is.na(opt$`gam-k`) && opt$`gam-k` < 2) {
      errors <- c(errors, "GAM basis dimension k (--gam-k) must be >= 2")
    }
  }
  
  # Output directory creation
  if (!dir.exists(opt$outdir)) {
    tryCatch({
      dir.create(opt$outdir, recursive = TRUE)
      msg("Created output directory: %s\n", opt$outdir)
    }, error = function(e) {
      errors <<- c(errors, sprintf("Cannot create output directory '%s': %s", opt$outdir, e$message))
    })
  }
  
  if (length(errors) > 0) {
    cat("VALIDATION ERRORS:\n", file = stderr())
    for (i in seq_along(errors)) {
      cat(sprintf("  %d. %s\n", i, errors[i]), file = stderr())
    }
    cat("\nUse --help for usage examples and argument details.\n", file = stderr())
    quit(status = 1)
  }
  
  # Warnings for deprecated options
  if (!is.null(opt$zcols) && is.null(opt$`random-effects`)) {
    msg("WARNING: --zcols is deprecated. Use --random-effects instead.\n")
  }
}

# Run validation (skip for --print-cols which has its own validation)
if (!isTRUE(opt$`print-cols`)) {
  validate_arguments(opt)
}

# For --print-cols we only need to read files
if (isTRUE(opt$`print-cols`)) {
  if (is.null(opt$qdec) || is.null(opt$aseg)) {
    cat("ERROR: --print-cols requires both --qdec and --aseg arguments\n", file = stderr())
    cat("Use --help for usage examples.\n", file = stderr())
    quit(status = 1)
  }
  if (!file.exists(opt$qdec)) {
    cat(sprintf("ERROR: qdec file not found: %s\n", opt$qdec), file = stderr())
    quit(status = 1)
  }
  if (!file.exists(opt$aseg)) {
    cat(sprintf("ERROR: aseg/aparc table not found: %s\n", opt$aseg), file = stderr())
    quit(status = 1)
  }
  
  msg("Reading files to display column information...\n")
  qdec <- tryCatch(read.delim(opt$qdec, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(qdec)) qdec <- read.table(opt$qdec, header=TRUE, stringsAsFactors=FALSE)
  aseg <- read.table(opt$aseg, header=TRUE, stringsAsFactors=FALSE)
  cat("qdec columns:\n"); print(names(qdec))
  cat("aseg/aparc columns:\n"); print(names(aseg))
  cat("\nFirst few rows of qdec:\n"); print(head(qdec))
  cat("\nFirst few rows of aseg:\n"); print(head(aseg))
  cat(sprintf("\nqdec dimensions: %d rows x %d cols\n", nrow(qdec), ncol(qdec)))
  cat(sprintf("aseg dimensions: %d rows x %d cols\n", nrow(aseg), ncol(aseg)))
  quit(status=0)
}

# Read inputs (validation already ensures files exist)
msg("Reading input files...\n")
qdec <- tryCatch(read.delim(opt$qdec, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
if (is.null(qdec)) qdec <- read.table(opt$qdec, header=TRUE, stringsAsFactors=FALSE)

coerce_numeric_like <- function(x) {
  if (!is.character(x)) return(x)
  stripped <- trimws(tolower(x))
  stripped[stripped %in% c("", "na", "n/a", "nan", "null")] <- NA_character_
  nums <- suppressWarnings(as.numeric(stripped))
  if (any(!is.na(stripped) & is.na(nums))) {
    x[is.na(stripped)] <- NA_character_
    return(x)
  }
  x[is.na(stripped)] <- NA
  nums
}

if (nrow(qdec) > 0) {
  for (nm in names(qdec)) {
    qdec[[nm]] <- coerce_numeric_like(qdec[[nm]])
  }
}

aseg <- read.table(opt$aseg, header=TRUE, stringsAsFactors=FALSE)

# Standardize qdec id column
if ("fsid-base" %in% names(qdec)) names(qdec)[names(qdec)=="fsid-base"] <- "fsid_base"
if ("fsid.base" %in% names(qdec)) names(qdec)[names(qdec)=="fsid.base"] <- "fsid_base"

# Determine id column in aseg and derive fsid/fsid_base
id_col <- if (!is.null(opt$id_col) && opt$id_col %in% names(aseg)) opt$id_col else {
  if ("Measure.volume" %in% names(aseg)) "Measure.volume" else if ("Measure:volume" %in% names(aseg)) "Measure:volume" else names(aseg)[1]
}
ids <- as.character(aseg[[id_col]])
has_long <- grepl("\\.long\\.", ids, perl=TRUE)
aseg$fsid <- ifelse(has_long, sub("^(.*)\\.long\\..*$", "\\1", ids, perl=TRUE), ids)
aseg$fsid_base <- ifelse(has_long, sub("^.*\\.long\\.(.*)$", "\\1", ids, perl=TRUE), NA_character_)

# Time column
time_col <- if (!is.null(opt$time_col)) opt$time_col else if ("tp" %in% names(qdec)) "tp" else stop("time_col not specified and 'tp' not in qdec")

# Keys check
if (!all(c("fsid","fsid_base") %in% names(qdec))) stop(sprintf("qdec must have fsid and fsid_base; has: %s", paste(names(qdec), collapse=", ")))
if (!all(c("fsid","fsid_base") %in% names(aseg))) stop(sprintf("aseg missing derived fsid/fsid_base; has: %s", paste(names(aseg), collapse=", ")))

# Merge
dat <- merge(qdec, aseg, by=c("fsid","fsid_base"))
if (nrow(dat) == 0) stop("Merged data is empty; check IDs")

# Order and ni
if (!time_col %in% names(dat)) stop(sprintf("time_col '%s' not found after merge", time_col))
dat <- dat[order(dat$fsid_base, dat[[time_col]]), ]

# Optionally derive interpretable covariates from group_5 for intervention design decomposition
if (isTRUE(opt$derive_group5) && ("group_5" %in% names(dat))) {
  g5 <- as.character(dat$group_5)
  is_na <- is.na(g5)
  dat$interv <- as.integer(!is_na & g5 != "control")
  dat$smallgroup <- as.integer(!is_na & grepl("^smallgroup", g5))
  dat$weeks4 <- as.integer(!is_na & grepl("_4w$", g5))
  dat$weeks <- ifelse(!is_na & grepl("_4w$", g5), 4L,
                      ifelse(!is_na & grepl("_2w$", g5), 2L, 0L))
}

# ROI selection
multi_region <- isTRUE(opt$`all-regions`) || !is.null(opt$`region-pattern`)
rois_to_analyze <- NULL
if (multi_region) {
  exclude <- unique(c("fsid","fsid_base", time_col, names(qdec), id_col))
  brain_cols <- setdiff(names(dat), exclude)
  # Drop summary/global measures unless explicitly included
  if (!isTRUE(opt$`include-summary`)) {
    summary_patterns <- c(
      "^BrainSegVol$", "^BrainSegVolNotVent$", "^CortexVol$",
      "^lhCortexVol$", "^rhCortexVol$", "^CerebralWhiteMatterVol$",
      "^lhCerebralWhiteMatterVol$", "^rhCerebralWhiteMatterVol$",
      "^TotalGrayVol$", "^SubCortGrayVol$", "^SupraTentorialVol$",
      "^SupraTentorialVolNotVent$", "^MaskVol$", "eTIV", "to\\.eTIV$",
      "^EstimatedTotalIntraCranialVol$", "^WM\\.hypointensities$",
      "hypointensities$", "^CC_", "^Optic\\.Chiasm$", "^CSF$",
      "^X3rd\\.Ventricle$", "^X4th\\.Ventricle$", "^X5th\\.Ventricle$"
    )
    keep <- !Reduce(`|`, lapply(summary_patterns, function(p) grepl(p, brain_cols, perl=TRUE)))
    brain_cols <- brain_cols[keep]
  }
  if (isTRUE(opt$`all-regions`)) {
    rois_to_analyze <- brain_cols
  } else {
    rois_to_analyze <- brain_cols[grepl(opt$`region-pattern`, brain_cols, ignore.case=TRUE, perl=TRUE)]
  }
  # Always exclude analysis/design helper columns that are not brain ROIs
  rois_to_analyze <- setdiff(rois_to_analyze, c(
    "baseline_value", "interv", "smallgroup", "weeks4", "weeks"
  ))
  if (!length(rois_to_analyze)) stop("No ROIs found to analyze")
} else {
  if (is.null(opt$roi)) stop("ROI is required (or use --all-regions/--region-pattern)")
  rois_to_analyze <- opt$roi
}

analyze_roi <- function(roi_name, dat, opt, time_col) {
  if (!(roi_name %in% names(dat))) {
    roi2 <- gsub("-", ".", roi_name, fixed=TRUE)
    if (!(roi2 %in% names(dat))) return(list(error=sprintf("ROI '%s' not found", roi_name)))
    roi_name <- roi2
  }
  form <- as.formula(opt$formula)
  # If requested or referenced in the formula, add per-subject baseline covariate for this ROI at the full-data level
  if (isTRUE(opt$add_baseline) || ("baseline_value" %in% all.vars(form))) {
    if (!("baseline_value" %in% names(dat))) {
      tp_vec_all <- dat[[time_col]]
      yvec_all <- suppressWarnings(as.numeric(dat[[roi_name]]))
      bases_all <- unique(dat$fsid_base)
      base_time_all <- tapply(tp_vec_all, dat$fsid_base, function(x) suppressWarnings(min(x, na.rm=TRUE)))
      base_val_all <- vapply(bases_all, function(sb) {
        idx <- which(dat$fsid_base == sb)
        if (!length(idx)) return(NA_real_)
        tmin <- base_time_all[[sb]]
        if (is.infinite(tmin) || is.na(tmin)) return(NA_real_)
        i_pick <- idx[which.min(abs(tp_vec_all[idx] - tmin))]
        as.numeric(yvec_all[i_pick])
      }, numeric(1))
      dat$baseline_value <- as.numeric(setNames(base_val_all, bases_all)[dat$fsid_base])
    }
  }
  needed_vars <- unique(c(all.vars(form), roi_name))
  missing_vars <- setdiff(needed_vars, names(dat))
  if (length(missing_vars)) {
    return(list(error=sprintf("Missing columns for ROI '%s': %s", roi_name, paste(missing_vars, collapse=", "))))
  }

  keep_mask <- stats::complete.cases(dat[, needed_vars, drop=FALSE])
  if (!any(keep_mask)) {
    return(list(error=sprintf("ROI '%s' has no complete cases for required covariates", roi_name)))
  }

  dat_roi <- dat[keep_mask, , drop=FALSE]
  dat_roi <- dat_roi[order(dat_roi$fsid_base, dat_roi[[time_col]]), , drop=FALSE]

  # dat_roi already inherits baseline_value from dat if computed

  ni_vec <- table(dat_roi$fsid_base)
  if (any(ni_vec < 2)) {
    keep_bases <- names(ni_vec)[ni_vec >= 2]
    dat_roi <- dat_roi[dat_roi$fsid_base %in% keep_bases, , drop=FALSE]
    if (!length(keep_bases)) {
      return(list(error=sprintf("ROI '%s' lacks subjects with >=2 timepoints after filtering", roi_name)))
    }
    ni_vec <- table(dat_roi$fsid_base)
  }

  yvec <- as.numeric(dat_roi[[roi_name]])
  if (all(is.na(yvec)) || sd(yvec, na.rm=TRUE) == 0) {
    return(list(error=sprintf("ROI '%s' has zero variance or all NA after filtering; skipped", roi_name)))
  }
  Y <- matrix(yvec, ncol=1)

  X <- model.matrix(form, dat_roi)
  if (nrow(X) != nrow(Y)) {
    return(list(error=sprintf("Design/response mismatch for ROI '%s': %d rows in X vs %d in Y", roi_name, nrow(X), nrow(Y))))
  }

  ni_local <- matrix(as.integer(ni_vec), ncol=1)
  # Compact progress indicator
  if (!isTRUE(opt$quiet)) cat(sprintf("%s: ", roi_name))
  flush.console()

  engine <- tolower(opt$engine)
  # Wrapper to optionally silence optimizer output and choose engine
  fit_fn <- switch(engine,
    fslmer = {
      if (!requireNamespace("fslmer", quietly=TRUE)) return(list(error="Package 'fslmer' is required for --engine fslmer"))
      # Determine zcols from indices or human-friendly spec
      parse_zcols <- function(spec, X, time_col) {
        cn <- colnames(X)
        # helper to find time slope columns (main effect only)
        time_idx <- integer(0)
        # factor-coded time columns like factor(tp)2, factor(tp)3
        fac_pat <- paste0("^factor\\(", gsub("[.\\[\\]\\^$]", "\\\\\\0", time_col), "\\)\\d+$")
        m_fac <- grepl(fac_pat, cn, perl=TRUE)
        if (any(m_fac)) {
          time_idx <- which(m_fac)
        } else {
          # numeric/continuous time column named exactly time_col
          m_num <- which(cn == time_col)
          if (length(m_num)) time_idx <- m_num
        }
        spec_u <- toupper(trimws(as.character(spec)))
        # If looks like indices, return parsed ints
        if (grepl("^[0-9, ]+$", spec_u)) {
          vals <- suppressWarnings(as.integer(strsplit(spec_u, ",")[[1]]))
          return(vals)
        }
        # Map synonyms to canonical
        spec_u <- switch(spec_u,
          "RIFS" = "RI",
          "FIFS" = "RI",   # no random effects in fslmer makes little sense; map to RI
          "FIRS" = "RS",
          spec_u
        )
        if (spec_u %in% c("RI", "RANDOM_INTERCEPT")) {
          return(1L)
        } else if (spec_u %in% c("RS", "RANDOM_SLOPE")) {
          return(time_idx)
        } else if (spec_u %in% c("RIRS", "RI+RS", "RI_RS", "RANDOM_INTERCEPT_RANDOM_SLOPE")) {
          return(unique(c(1L, time_idx)))
        } else {
          # unknown token -> try to parse as integers anyway
          vals <- suppressWarnings(as.integer(strsplit(spec_u, ",")[[1]]))
          return(vals)
        }
      }
      zcols <- if (is.numeric(opt$zcols)) as.integer(opt$zcols) else parse_zcols(opt$zcols, X, time_col)
      if (length(zcols) == 0 || any(is.na(zcols))) return(list(error="--random-effects (or --zcols) must be indices or one of RI/RIRS/RS (synonyms: RIFS->RI, FIRS->RS, FIFS->RI)"))
      function() fslmer::lme_fit_FS(X, zcols, Y, ni_local)
    },
    glm = {
      dat_roi$y <- as.numeric(dat_roi[[roi_name]])
      fam <- tryCatch(get(opt$family, mode="function"), error=function(e) NULL)
      if (is.null(fam)) fam <- gaussian
      fam <- fam()
      function() stats::glm(update(form, y ~ .), data=dat_roi, family=fam)
    },
    gam = {
      if (!requireNamespace("mgcv", quietly=TRUE)) return(list(error="Package 'mgcv' is required for --engine gam"))
      dat_roi$y <- as.numeric(dat_roi[[roi_name]])
      # Ensure time variable is numeric for s(time). If not, try to coerce sensibly.
      tp_raw <- dat_roi[[time_col]]
      if (!is.numeric(tp_raw)) {
        # Try to extract numeric part (e.g., ses-1 -> 1)
        tp_num <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", as.character(tp_raw))))
        if (all(is.na(tp_num))) {
          # Fallback: ordered factor mapping to 1..k by sorted unique levels
          lev <- sort(unique(as.character(tp_raw)))
          f <- factor(tp_raw, levels=lev, ordered=TRUE)
          tp_num <- as.numeric(f)
        }
        dat_roi[[time_col]] <- tp_num
      }
      # After coercion, validate there are at least 2 unique time points
      uniq_t <- length(unique(dat_roi[[time_col]][!is.na(dat_roi[[time_col]])]))
      if (uniq_t < 2) return(list(error=sprintf("Not enough unique values in '%s' for GAM (need >=2)", time_col)))
      fam <- tryCatch(get(opt$family, mode="function"), error=function(e) NULL)
      if (is.null(fam)) fam <- gaussian
      fam <- fam()
      # Choose k adaptively unless user specified --gam-k
      k_user <- if (!is.na(opt$`gam-k`)) as.integer(opt$`gam-k`) else NA_integer_
      k_time <- if (!is.na(k_user)) k_user else max(2L, min(5L, uniq_t))
      if (!isTRUE(opt$quiet) && isTRUE(opt$debug)) cat(sprintf("[GAM] Using k=%d for s(%s) with %d unique values\n", k_time, time_col, uniq_t))
      base_form <- paste0("y ~ s(", time_col, ", k=", k_time, ")")
      if (!isTRUE(opt$`no-gam-re`)) base_form <- paste0(base_form, "+ s(fsid_base, bs='re')")
      gform <- as.formula(base_form)
      function() {
        # Try with chosen k; if the common mgcv error about unique covariate combinations occurs, reduce k and retry
        fit <- try(mgcv::gam(gform, data=dat_roi, family=fam, method="REML"), silent=TRUE)
        if (inherits(fit, "try-error")) {
          msg_txt <- as.character(fit)
          if (grepl("fewer unique covariate combinations", msg_txt, fixed=TRUE)) {
            k_fallback <- max(2L, min(k_time, uniq_t))
            if (k_fallback < k_time) {
              gform_fb <- as.formula(paste0("y ~ s(", time_col, ", k=", k_fallback, ")", if (!isTRUE(opt$`no-gam-re`)) "+ s(fsid_base, bs='re')" else ""))
              if (!isTRUE(opt$quiet) && isTRUE(opt$debug)) cat(sprintf("[GAM] Retrying with k=%d due to limited unique values\n", k_fallback))
              return(mgcv::gam(gform_fb, data=dat_roi, family=fam, method="REML"))
            }
          }
          stop(attr(fit, "condition")$message)
        }
        fit
      }
    },
    {
      return(list(error=sprintf("Unknown --engine '%s' (use fslmer|glm|gam)", opt$engine)))
    }
  )
  run_fit <- function() {
    if (isTRUE(opt$debug)) {
      fit_fn()
    } else {
      # Capture and discard optimizer chatter
      zz <- file(tempfile(), open="wt"); zz2 <- file(tempfile(), open="wt")
      on.exit({ try(close(zz), silent=TRUE); try(close(zz2), silent=TRUE) }, add=TRUE)
      sink(zz); sink(zz2, type="message")
      on.exit({ try(sink(NULL), silent=TRUE); try(sink(NULL, type="message"), silent=TRUE) }, add=TRUE)
      fit <- fit_fn()
      fit
    }
  }

  res <- tryCatch({
    fit <- run_fit()
    out <- list(roi=roi_name, engine=engine, X=X, Y=Y, ni=ni_local)
    if (engine == "fslmer") {
      F_C <- NULL; Cvec <- NULL
      # Robust contrast parsing: accept NULL, character, numeric; ignore empty lists/objects
      if (!is.null(opt$contrast)) {
        if (is.character(opt$contrast)) {
          parts <- strsplit(opt$contrast, ",")[[1]]
          Cvec <- suppressWarnings(as.numeric(parts))
        } else if (is.numeric(opt$contrast)) {
          Cvec <- as.numeric(opt$contrast)
        } else {
          Cvec <- NULL
        }
        if (!is.null(Cvec)) {
          if (length(Cvec) != ncol(X)) return(list(error=sprintf("Contrast length %d != ncol(X) %d", length(Cvec), ncol(X))))
          C <- matrix(Cvec, nrow=1)
          F_C <- fslmer::lme_F(fit, C)
        }
      }
      out$stats <- fit; out$F_C <- F_C; out$Cvec <- Cvec; out$zcols <- zcols
    } else if (engine == "glm") {
      out$glm <- fit
    } else if (engine == "gam") {
      out$gam <- fit
    }
    if (!isTRUE(opt$quiet)) cat("######## done\n")
    out
  }, error=function(e) list(error=sprintf("Model failed for %s: %s", roi_name, e$message)))
  res
}

results <- list(); failed <- character(0)
for (roi in rois_to_analyze) {
  r <- analyze_roi(roi, dat, opt, time_col)
  if (!is.null(r$error)) { cat("ERROR:", r$error, "\n"); failed <- c(failed, roi); next }
  results[[roi]] <- r
}
if (!length(results)) stop("No ROIs were successfully analyzed")

outdir <- opt$outdir; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
if (isTRUE(opt$save_merged)) write.csv(dat, file=file.path(outdir, "merged_data.csv"), row.names=FALSE)

if (multi_region) {
  coef_summary <- do.call(rbind, lapply(names(results), function(roi) {
    r <- results[[roi]]
    engine <- if (!is.null(r$engine)) r$engine else "fslmer"
    if (engine == "fslmer") {
      bhat <- r$stats$Bhat
      beta <- as.numeric(bhat)
      cn <- rownames(bhat)
      if (is.null(cn)) {
        cn <- colnames(r$X)
        if (is.null(cn) || length(cn) != length(beta)) cn <- paste0("beta", seq_along(beta))
      }
      data.frame(roi=roi, coef=cn, beta=beta, engine=engine, stringsAsFactors=FALSE)
    } else if (engine == "glm") {
      co <- tryCatch(stats::coef(r$glm), error=function(e) NULL)
      if (is.null(co)) return(NULL)
      data.frame(roi=roi, coef=names(co), beta=as.numeric(co), engine=engine, stringsAsFactors=FALSE)
    } else if (engine == "gam") {
      co <- tryCatch(stats::coef(r$gam), error=function(e) NULL)
      if (is.null(co)) return(NULL)
      data.frame(roi=roi, coef=names(co), beta=as.numeric(co), engine=engine, stringsAsFactors=FALSE)
    } else {
      NULL
    }
  }))
  write.csv(coef_summary, file=file.path(outdir, "lme_coefficients.csv"), row.names=FALSE)
  if (!is.null(opt$contrast)) {
    f_summary <- do.call(rbind, lapply(names(results), function(roi) {
      r <- results[[roi]]; if (!is.null(r$F_C)) data.frame(roi=roi, F=r$F_C$F, pval=r$F_C$pval, sgn=r$F_C$sgn, df=r$F_C$df)
    }))
    if (!is.null(f_summary)) write.csv(f_summary, file=file.path(outdir, "lme_F_test.csv"), row.names=FALSE)
  }
  roi_dir <- file.path(outdir, "individual_rois"); dir.create(roi_dir, showWarnings=FALSE, recursive=TRUE)
  for (roi in names(results)) {
    r <- results[[roi]]; safe <- gsub("[^A-Za-z0-9_-]", "_", roi)
    saveRDS(r$stats, file=file.path(roi_dir, paste0(safe, "_fit.rds")))
  }
} else {
  r <- results[[1]]; saveRDS(r$stats, file=file.path(outdir, "fit.rds"))
  if (!is.null(r$F_C)) write.table(data.frame(F=r$F_C$F, pval=r$F_C$pval, sgn=r$F_C$sgn, df=r$F_C$df), file=file.path(outdir, "F_test.txt"), quote=FALSE, sep="\t", row.names=FALSE)
}

resolved <- list(
  qdec=opt$qdec, aseg=opt$aseg,
  roi=if(multi_region) paste(rois_to_analyze, collapse=",") else opt$roi,
  formula=opt$formula,
  random_effects=if (!is.null(opt$`random-effects`)) opt$`random-effects` else NULL,
  zcols=if(length(results)) results[[1]]$zcols else NULL,
  contrast=if(length(results) && !is.null(results[[1]]$Cvec)) results[[1]]$Cvec else NULL,
  outdir=outdir, time_col=time_col, id_col=id_col, multi_region=multi_region,
  n_rois_analyzed=length(results), n_rois_failed=length(failed)
)
jsonlite::write_json(resolved, file.path(outdir, "used_config.json"), pretty=TRUE, auto_unbox=TRUE)

cat(sprintf("Done. Outputs in %s. Success: %d, Failed: %d\n", outdir, length(results), length(failed)))

# Optional note about optimizer logs (only if --debug)
if (isTRUE(opt$debug) && !isTRUE(opt$quiet)) {
  cat("\nNote about optimizer progress:\n")
  cat("- Lines like 'Likelihood at FS iteration k : <value>' come from fslmer's Fisher Scoring optimizer.\n")
  cat("- 'Likelihood' is the model log-likelihood; it should increase (become less negative) and stabilize as the fit converges.\n")
  cat("- 'Gradient norm' is the size of the score vector; small values (e.g., < 1e-3) indicate convergence.\n\n")
}
