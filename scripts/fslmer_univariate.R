#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  # GLM base
  library(stats)
  # GAM (optional); load lazily later to avoid hard dependency
})

# Flexible univariate LME analysis for aseg/aparc tables using fslmer

option_list <- list(
  make_option(c("--config"), type="character", default=NULL, help="Path to JSON config file with arguments"),
  make_option(c("-q", "--qdec"), type="character", help="Path to qdec.table.dat (TSV)"),
  make_option(c("-a", "--aseg"), type="character", help="Path to aseg/aparc long table (TSV/whitespace)"),
  make_option(c("-r", "--roi"), type="character", help="ROI column name in aseg/aparc table (e.g., Left.Hippocampus)"),
  make_option(c("-f", "--formula"), type="character", default="~ tp", help="Model formula for fixed effects (e.g., '~ tp*group') [default %default]"),
  make_option(c("-z", "--zcols"), type="character", default="1,2", help="Comma-separated columns for random effects (e.g., '1,2' for intercept+time) [default %default]"),
  make_option(c("-c", "--contrast"), type="character", default=NULL, help="Comma-separated contrast vector (length must match ncol(X))"),
  make_option(c("-o", "--outdir"), type="character", default="fslmer_out", help="Output directory [default %default]"),
  make_option(c("--save-merged"), action="store_true", default=FALSE, help="Save merged dat as CSV"),
  make_option(c("--time-col"), type="character", default=NULL, help="Name of time variable in qdec; if NULL and 'tp' exists, uses 'tp'"),
  make_option(c("--id-col"), type="character", default=NULL, help="Override ID column in aseg (default autodetect 'Measure:volume')"),
  make_option(c("--print-cols"), action="store_true", default=FALSE, help="Print column names (and heads) of qdec and aseg then exit"),
  make_option(c("--all-regions"), action="store_true", default=FALSE, help="Run analysis on all brain regions (subcortical volumes)"),
  make_option(c("--region-pattern"), type="character", default=NULL, help="Regex to match ROI names (e.g., 'Hippocampus|Amygdala')"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity"),
  make_option(c("--include-summary"), action="store_true", default=FALSE, help="Include global/summary volume measures in multi-region analysis"),
  make_option(c("--engine"), type="character", default="fslmer", help="Model engine: 'fslmer' (default), 'glm', or 'gam'"),
  make_option(c("--family"), type="character", default="gaussian", help="GLM/GAM family (gaussian, binomial, poisson, Gamma, etc.) [ignored for fslmer]"),
  make_option(c("--gam-k"), type="integer", default=NA, help="Basis dimension k for s(time) in GAM [default auto: min(5, unique time levels); must be >=2]"),
  make_option(c("--no-gam-re"), action="store_true", default=FALSE, help="Disable random intercept s(fsid_base, bs='re') in GAM"),
  make_option(c("--debug"), action="store_true", default=FALSE, help="Print detailed optimizer logs from fslmer")
)

opt <- parse_args(OptionParser(option_list=option_list))
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
    contrast = NULL,
    outdir = "fslmer_out",
    time_col = NULL,
    id_col = NULL
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
  opt$contrast <- apply_cfg("contrast")
  opt$outdir   <- apply_cfg("outdir")
  opt$time_col <- apply_cfg("time_col")
  opt$id_col   <- apply_cfg("id_col")
}

# For --print-cols we only need to read files
if (isTRUE(opt$`print-cols`)) {
  if (is.null(opt$qdec) || is.null(opt$aseg)) stop("qdec and aseg paths are required for --print-cols")
  if (!file.exists(opt$qdec)) stop(sprintf("qdec file not found: %s", opt$qdec))
  if (!file.exists(opt$aseg)) stop(sprintf("aseg/aparc table not found: %s", opt$aseg))
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

if (is.null(opt$qdec) || is.null(opt$aseg)) stop("qdec and aseg paths are required (via flags or config)")
if (!file.exists(opt$qdec)) stop(sprintf("qdec file not found: %s", opt$qdec))
if (!file.exists(opt$aseg)) stop(sprintf("aseg/aparc table not found: %s", opt$aseg))

# Read inputs
qdec <- tryCatch(read.delim(opt$qdec, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
if (is.null(qdec)) qdec <- read.table(opt$qdec, header=TRUE, stringsAsFactors=FALSE)

coerce_numeric_like <- function(x) {
  if (!is.character(x)) return(x)
  stripped <- trimws(tolower(x))
  stripped[stripped %in% c("", "na", "n/a", "nan", "null")] <- NA_character_
  nums <- suppressWarnings(as.numeric(stripped))
  if (any(!is.na(stripped) & is.na(nums))) {
    return(x)
  }
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
ni <- matrix(unname(table(dat$fsid_base)), ncol=1)

# ROI selection
multi_region <- isTRUE(opt$`all-regions`) || !is.null(opt$`region-pattern`)
rois_to_analyze <- NULL
if (multi_region) {
  exclude <- unique(c("fsid","fsid_base", time_col, names(qdec)))
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
  if (!length(rois_to_analyze)) stop("No ROIs found to analyze")
} else {
  if (is.null(opt$roi)) stop("ROI is required (or use --all-regions/--region-pattern)")
  rois_to_analyze <- opt$roi
}

analyze_roi <- function(roi_name, dat, ni, opt, time_col) {
  if (!(roi_name %in% names(dat))) {
    roi2 <- gsub("-", ".", roi_name, fixed=TRUE)
    if (!(roi2 %in% names(dat))) return(list(error=sprintf("ROI '%s' not found", roi_name)))
    roi_name <- roi2
  }
  Y <- matrix(dat[[roi_name]], ncol=1)
  # Skip degenerate ROIs (no variance) or mostly missing
  yvec <- as.numeric(dat[[roi_name]])
  if (all(is.na(yvec)) || sd(yvec, na.rm=TRUE) == 0) {
    return(list(error=sprintf("ROI '%s' has zero variance or all NA; skipped", roi_name)))
  }
  form <- as.formula(opt$formula)
  X <- model.matrix(form, dat)
  # Compact progress indicator
  if (!isTRUE(opt$quiet)) cat(sprintf("%s: ", roi_name))
  flush.console()

  engine <- tolower(opt$engine)
  # Wrapper to optionally silence optimizer output and choose engine
  fit_fn <- switch(engine,
    fslmer = {
      if (!requireNamespace("fslmer", quietly=TRUE)) return(list(error="Package 'fslmer' is required for --engine fslmer"))
      zcols <- if (is.numeric(opt$zcols)) as.integer(opt$zcols) else as.integer(strsplit(opt$zcols, ",")[[1]])
      if (any(is.na(zcols))) return(list(error="--zcols must be integers"))
      function() fslmer::lme_fit_FS(X, zcols, Y, ni)
    },
    glm = {
      dat$y <- as.numeric(dat[[roi_name]])
      fam <- tryCatch(get(opt$family, mode="function"), error=function(e) NULL)
      if (is.null(fam)) fam <- gaussian
      fam <- fam()
      function() stats::glm(update(form, y ~ .), data=dat, family=fam)
    },
    gam = {
      if (!requireNamespace("mgcv", quietly=TRUE)) return(list(error="Package 'mgcv' is required for --engine gam"))
      dat$y <- as.numeric(dat[[roi_name]])
      # Ensure time variable is numeric for s(time). If not, try to coerce sensibly.
      tp_raw <- dat[[time_col]]
      if (!is.numeric(tp_raw)) {
        # Try to extract numeric part (e.g., ses-1 -> 1)
        tp_num <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", as.character(tp_raw))))
        if (all(is.na(tp_num))) {
          # Fallback: ordered factor mapping to 1..k by sorted unique levels
          lev <- sort(unique(as.character(tp_raw)))
          f <- factor(tp_raw, levels=lev, ordered=TRUE)
          tp_num <- as.numeric(f)
        }
        dat[[time_col]] <- tp_num
      }
      # After coercion, validate there are at least 2 unique time points
      uniq_t <- length(unique(dat[[time_col]][!is.na(dat[[time_col]])]))
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
        fit <- try(mgcv::gam(gform, data=dat, family=fam, method="REML"), silent=TRUE)
        if (inherits(fit, "try-error")) {
          msg_txt <- as.character(fit)
          if (grepl("fewer unique covariate combinations", msg_txt, fixed=TRUE)) {
            k_fallback <- max(2L, min(k_time, uniq_t))
            if (k_fallback < k_time) {
              gform_fb <- as.formula(paste0("y ~ s(", time_col, ", k=", k_fallback, ")", if (!isTRUE(opt$`no-gam-re`)) "+ s(fsid_base, bs='re')" else ""))
              if (!isTRUE(opt$quiet) && isTRUE(opt$debug)) cat(sprintf("[GAM] Retrying with k=%d due to limited unique values\n", k_fallback))
              return(mgcv::gam(gform_fb, data=dat, family=fam, method="REML"))
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
    out <- list(roi=roi_name, engine=engine, X=X, Y=Y)
    if (engine == "fslmer") {
      F_C <- NULL; Cvec <- NULL
      if (!is.null(opt$contrast)) {
        Cvec <- if (is.numeric(opt$contrast)) as.numeric(opt$contrast) else as.numeric(strsplit(opt$contrast, ",")[[1]])
        if (length(Cvec) != ncol(X)) return(list(error=sprintf("Contrast length %d != ncol(X) %d", length(Cvec), ncol(X))))
        C <- matrix(Cvec, nrow=1)
        F_C <- fslmer::lme_F(fit, C)
      }
      out$stats <- fit; out$F_C <- F_C; out$Cvec <- Cvec
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
  r <- analyze_roi(roi, dat, ni, opt, time_col)
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
  formula=opt$formula, zcols=if(length(results)) results[[1]]$zcols else NULL,
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
