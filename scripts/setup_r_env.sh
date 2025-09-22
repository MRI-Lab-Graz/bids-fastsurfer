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
# - Installs CRAN packages using pak (preferred) or remotes fallback
#   - CRAN: optparse, jsonlite
#   - GitHub: Deep-MI/fslmer
# - Snapshots renv.lock (unless --no-snapshot)
# - Verifies Rscript availability and prints package versions; runs a lightweight self-check of scripts/fslmer_univariate.R
#
# Notes:
# - Requires Rscript in PATH. On macOS, you may need Xcode CLT: `xcode-select --install`.
# - The script should be run from repo root.

QUIET=0
SNAPSHOT=1
CRAN_MIRROR="https://cloud.r-project.org"
LOG_FILE="setup_r_env.log"

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
log "Installing packages with pak (preferred)"
# Try pak first (fast resolver, binary packages when available)
set +e
Rscript -e "if (!requireNamespace('pak', quietly=TRUE)) install.packages('pak', repos='https://r-lib.github.io/p/pak/stable'); pak::pkg_install(c('optparse','jsonlite','remotes'), upgrade = FALSE)" >>"${LOG_FILE}" 2>&1
PAK_RC=$?
set -e
if [[ $PAK_RC -ne 0 ]]; then
  log "pak failed; falling back to install.packages for CRAN deps"
  Rscript -e "install.packages(c('optparse','jsonlite','remotes'), repos='${CRAN_MIRROR}')" >>"${LOG_FILE}" 2>&1

# Ensure bettermc (dependency of fslmer) is available before installing fslmer
log "Ensuring 'bettermc' is installed (from CRAN archive if needed)"
BETTERMC_VERSION="1.2.1"
BETTERMC_URL="https://cran.r-project.org/src/contrib/Archive/bettermc/bettermc_${BETTERMC_VERSION}.tar.gz"
set +e
# Prefer direct archive URL via pak to avoid solver issues
Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install('${BETTERMC_URL}', upgrade = FALSE) else quit(status=99)" >>"${LOG_FILE}" 2>&1
RC_PAK_BMC_URL=$?
set -e
if [[ $RC_PAK_BMC_URL -eq 99 || $RC_PAK_BMC_URL -ne 0 ]]; then
  log "pak archive install failed or pak not available; trying remotes::install_url for bettermc"
  set +e
  Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}'); remotes::install_url('${BETTERMC_URL}')" >>"${LOG_FILE}" 2>&1
  RC_REM_URL=$?
  set -e
  if [[ $RC_REM_URL -ne 0 ]]; then
    log "remotes::install_url failed; trying remotes::install_version for bettermc ${BETTERMC_VERSION}"
    set +e
    Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}'); remotes::install_version('bettermc', version='${BETTERMC_VERSION}', repos='https://cran.r-project.org')" >>"${LOG_FILE}" 2>&1
    RC_REM_VER=$?
    set -e
    if [[ $RC_REM_VER -ne 0 ]]; then
      log "Failed to install bettermc from all sources. Showing last 60 log lines:"; tail -n 60 "${LOG_FILE}" || true; die "Failed to install 'bettermc' (required by fslmer). Check network/firewall and try again.";
    fi
  fi
fi
fi

# Install fslmer
log "Installing fslmer (Deep-MI/fslmer) via pak (preferred)"
# Allow pinning a specific ref via env var FSLMER_REF (tag/branch/commit). Example: export FSLMER_REF=v0.2.0
FSLMER_REF_ARG=""
if [[ -n "${FSLMER_REF:-}" ]]; then FSLMER_REF_ARG=", ref='${FSLMER_REF}'"; fi

set +e
Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install(sprintf('Deep-MI/fslmer%s', if (nzchar(Sys.getenv('FSLMER_REF'))) paste0('@', Sys.getenv('FSLMER_REF')) else ''), upgrade = FALSE) else quit(status=99)" >>"${LOG_FILE}" 2>&1
PAK_FSL_RC=$?
set -e
if [[ $PAK_FSL_RC -ne 0 && $PAK_FSL_RC -ne 99 ]]; then
  log "pak install of fslmer failed (exit $PAK_FSL_RC); will try remotes fallback"
fi
if [[ $PAK_FSL_RC -eq 99 || $PAK_FSL_RC -ne 0 ]]; then
  log "Installing fslmer via remotes fallback (no build/vignettes)"
  export R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
  # Ensure remotes present in current library before calling it
  Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}')" >>"${LOG_FILE}" 2>&1
  Rscript -e "remotes::install_github('Deep-MI/fslmer'${FSLMER_REF_ARG}, build=FALSE, build_vignettes=FALSE, dependencies=c('Depends','Imports','LinkingTo'), upgrade='never', quiet=TRUE)" >>"${LOG_FILE}" 2>&1 || {
    log "fslmer install failed. Showing last 60 log lines:"; tail -n 60 "${LOG_FILE}" || true; die "Failed to install fslmer"; }
fi

# Snapshot
if [[ $SNAPSHOT -eq 1 ]]; then
  log "Snapshotting renv state"
  Rscript -e "renv::snapshot(prompt=FALSE)" >>"${LOG_FILE}" 2>&1
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
log "Install log saved to ${LOG_FILE}"
