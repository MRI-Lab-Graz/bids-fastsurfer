#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(fslmer)
})

# Flexible univariate LME analysis for aseg/aparc tables using fslmer
# Inputs:
#  - Qdec TSV (from scripts/generate_qdec.py)
#  - aseg/aparc table (from asegstats2table/aparcstats2table with --qdec-long)
#  - Model specification via R formula string (e.g., "~ tp*group")
#  - Random effects columns (Zcols); default: c(1,2) for intercept and first time-like term
#  - ROI column name to analyze (e.g., "Left.Hippocampus")
#  - Outputs: RDS of fit, text summary, optional CSV of merged data

option_list <- list(
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
  make_option(c("--print-cols"), action="store_true", default=FALSE, help="Print column names of qdec and aseg then exit"),
  make_option(c("--quiet"), action="store_true", default=FALSE, help="Reduce output verbosity")
)

opt <- parse_args(OptionParser(option_list=option_list))

msg <- function(...) { if (!isTRUE(opt$quiet)) cat(sprintf(...), sep="") }

stopifnot(!is.null(opt$qdec), !is.null(opt$aseg))
if (is.null(opt$roi)) stop("--roi is required (e.g., Left.Hippocampus)")

if (!file.exists(opt$qdec)) stop(sprintf("qdec file not found: %s", opt$qdec))
if (!file.exists(opt$aseg)) stop(sprintf("aseg/aparc table not found: %s", opt$aseg))

# Read inputs
qdec <- tryCatch(read.delim(opt$qdec, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
if (is.null(qdec)) qdec <- read.table(opt$qdec, header=TRUE, stringsAsFactors=FALSE)
aseg <- read.table(opt$aseg, header=TRUE, stringsAsFactors=FALSE)

if (isTRUE(opt$print_cols)) {
  cat("qdec columns:\n"); print(names(qdec))
  cat("aseg/aparc columns:\n"); print(names(aseg))
  quit(status=0)
}

# Identify ID column in aseg
id_col <- opt$id_col
if (is.null(id_col)) {
  cand <- grep("^Measure[.:]volume$", names(aseg), ignore.case=TRUE, value=TRUE)
  if (length(cand) == 0) stop("Could not find 'Measure:volume' column in aseg/aparc table; use --id-col to override")
  id_col <- cand[[1]]
}

# Extract fsid and fsid.base
aseg$fsid      <- sub("\\.long\\..*", "",  aseg[[id_col]])
aseg$fsid.base <- sub(".*\\.long\\.",  "", aseg[[id_col]])

# Ensure time column
time_col <- opt$time_col
if (is.null(time_col)) {
  if ("Time.From.Baseline" %in% names(qdec)) time_col <- "Time.From.Baseline"
  else if ("tp" %in% names(qdec)) time_col <- "tp"
  else stop("Could not infer time column. Set --time-col (e.g., tp or Time.From.Baseline)")
}
if (!(time_col %in% names(qdec))) stop(sprintf("Time column '%s' not found in qdec", time_col))

# Merge
dat <- merge(qdec, aseg, by=c("fsid.base", "fsid"))
if (nrow(dat) == 0) stop("Merged data is empty; check IDs and inputs")

# Order by subject then time
dat <- dat[order(dat$fsid.base, dat[[time_col]]), ]

# Build response Y and ni
if (!(opt$roi %in% names(dat))) {
  # Try to adapt: replace '-' with '.' in ROI name
  roi2 <- gsub("-", ".", opt$roi, fixed=TRUE)
  if (!(roi2 %in% names(dat))) stop(sprintf("ROI column '%s' not found in merged data", opt$roi))
  opt$roi <- roi2
}
Y <- matrix(dat[[opt$roi]], ncol=1)
ni <- matrix(unname(table(dat$fsid.base)), ncol=1)

# Design matrix from formula (fixed effects)
form <- as.formula(opt$formula)
X <- model.matrix(form, dat)

# Random effects columns
zcols <- as.integer(strsplit(opt$zcols, ",")[[1]])
if (any(is.na(zcols))) stop("--zcols must be a comma-separated list of integers")

msg("Design matrix: %d x %d\n", nrow(X), ncol(X))
msg("Random effects columns (Zcols): %s\n", paste(zcols, collapse=","))
msg("ROI: %s; observations: %d; subjects: %d\n", opt$roi, nrow(Y), length(ni))

# Fit model
stats <- lme_fit_FS(X, zcols, Y, ni)

# Optional contrast
F_C <- NULL
if (!is.null(opt$contrast)) {
  Cvec <- as.numeric(strsplit(opt$contrast, ",")[[1]])
  if (length(Cvec) != ncol(X)) stop(sprintf("Contrast length (%d) must match ncol(X) (%d)", length(Cvec), ncol(X)))
  C <- matrix(Cvec, nrow=1)
  F_C <- lme_F(stats, C)
}

# Outputs
outdir <- opt$outdir
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
saveRDS(stats, file=file.path(outdir, "fit.rds"))
write.table(data.frame(Bhat=I(stats$Bhat)), file=file.path(outdir, "Bhat.txt"), quote=FALSE, sep="\t")
if (!is.null(F_C)) {
  write.table(data.frame(F=F_C$F, pval=F_C$pval, sgn=F_C$sgn, df=F_C$df),
              file=file.path(outdir, "F_test.txt"), quote=FALSE, sep="\t", row.names=FALSE)
}
if (isTRUE(opt$save_merged)) {
  write.csv(dat, file=file.path(outdir, "merged_data.csv"), row.names=FALSE)
}

msg("Done. Outputs written to %s\n", outdir)
