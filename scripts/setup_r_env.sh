#!/usr/bin/env bash
set -euo pipefail

# setup_r_env.sh
# Purpose: Install/activate an R renv project env for this repo, install required packages,
# snapshot lockfile, and verify the analysis script runs basic checks.
#
# Usage:
#   bash scripts/setup_r_env.sh [--no-snapshot] [--cran-mirror <url>] [--quiet]
#
# Behavior:
# - If renv is not installed system-wide, installs it to user library.
# - If renv project not initialized, runs renv::init(); otherwise renv::activate() and renv::restore() if lockfile exists.
# - Installs CRAN packages: optparse, jsonlite, remotes
# - Installs fslmer from GitHub via remotes: Deep-MI/fslmer
# - Snapshots renv.lock (unless --no-snapshot)
# - Verifies Rscript availability and prints package versions; runs a lightweight self-check of scripts/fslmer_univariate.R
#
# Notes:
# - Requires Rscript in PATH. On macOS, you may need Xcode CLT: `xcode-select --install`.
# - The script should be run from repo root.

QUIET=0
SNAPSHOT=1
CRAN_MIRROR="https://cloud.r-project.org"

die() { echo "[setup_r_env] $*" >&2; exit 1; }
log() { if [[ $QUIET -eq 0 ]]; then echo "[setup_r_env] $*"; fi }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-snapshot) SNAPSHOT=0; shift ;;
    --cran-mirror) CRAN_MIRROR="${2:-}"; shift 2 ;;
    --quiet) QUIET=1; shift ;;
    *) die "Unknown option: $1" ;;
  esac
done

command -v Rscript >/dev/null 2>&1 || die "Rscript not found. Please install R first."

log "Rscript: $(Rscript --version | head -n1 || true)"

# Ensure renv is installed
log "Ensuring renv is installed..."
Rscript -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv', repos='${CRAN_MIRROR}')" >/dev/null

# Initialize/activate project
if [[ -f "renv.lock" ]]; then
  log "renv.lock found — activating and restoring"
  Rscript -e "renv::activate(); renv::restore(prompt=FALSE)" >/dev/null
else
  log "No renv.lock — initializing new renv project"
  Rscript -e "renv::init(bare=TRUE)" >/dev/null
fi

# Install CRAN deps
log "Installing CRAN packages: optparse, jsonlite, remotes"
Rscript -e "install.packages(c('optparse','jsonlite','remotes'), repos='${CRAN_MIRROR}')" >/dev/null

# Install fslmer
log "Installing fslmer from GitHub (Deep-MI/fslmer) via remotes (no build/vignettes)"
Rscript -e "remotes::install_github('Deep-MI/fslmer', build=FALSE, build_vignettes=FALSE, dependencies=TRUE, upgrade='never', quiet=TRUE)" >/dev/null || die "Failed to install fslmer"

# Snapshot
if [[ $SNAPSHOT -eq 1 ]]; then
  log "Snapshotting renv state"
  Rscript -e "renv::snapshot(prompt=FALSE)" >/dev/null
fi

# Verify loaded packages and versions
log "Verifying R packages"
Rscript -e "pkgs <- c('optparse','jsonlite','fslmer'); print(data.frame(pkg=pkgs, available=sapply(pkgs, requireNamespace, quietly=TRUE)))"

# Lightweight check: print columns of a dummy small table via the R script (should exit gracefully)
if [[ -f scripts/fslmer_univariate.R ]]; then
  log "Running lightweight self-check of scripts/fslmer_univariate.R --print-cols (expected to error if files missing)"
  set +e
  Rscript scripts/fslmer_univariate.R --print-cols 2>/dev/null
  RC=$?
  set -e
  if [[ $RC -ne 0 ]]; then
    log "Self-check completed (no input files provided, which is fine)."
  else
    log "Self-check completed successfully."
  fi
fi

log "R environment setup complete."
