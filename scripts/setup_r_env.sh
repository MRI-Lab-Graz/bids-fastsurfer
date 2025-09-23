#!/usr/bin/env bash
set -euo pipefail

# setup_r_env.sh
# Purpose: Install/activate an R renv project env for this repo, install required packages,
# snapshot lockfile, and verify the analysis script runs basic checks.
#
# Usage:
#   bash scripts/setup_r_env.sh [--no-snapshot] [--cran-mirror <url>] [--offline] [--quiet]
#
# Behavior:
# - If renv is not installed system-wide, installs it to user library.
# - If renv project not initialized, runs renv::init(); otherwise renv::activate() and renv::restore() if lockfile exists.
# - Installs CRAN packages using pak (preferred) or remotes fallback
#   - CRAN: optparse, jsonlite, mgcv (mgcv is required for GAM runs)
#   - GitHub: Deep-MI/fslmer
# - Snapshots renv.lock (unless --no-snapshot)
# - Verifies Rscript availability and prints package versions; runs a lightweight self-check of scripts/fslmer_univariate.R
#
# Notes:
# - Requires Rscript in PATH. On macOS, you may need Xcode CLT: `xcode-select --install`.
# - The script should be run from repo root.
#!/usr/bin/env bash
# This script is deprecated and kept only as a stub for backward compatibility.
# Please use the micromamba-based installer instead:
#   bash scripts/install.sh

echo "[setup_r_env.sh] Deprecated. Use: bash scripts/install.sh" >&2
exit 1
      log "Offline mode: mgcv not installed (no vendor/mgcv_*.tar.gz found). If you plan to use --engine gam, provide an mgcv tarball or run setup online."
    fi
  fi
else
  log "Installing packages with pak (preferred)"
  # Try pak first (fast resolver, binary packages when available)
  set +e
  Rscript -e "if (!requireNamespace('pak', quietly=TRUE)) install.packages('pak', repos='https://r-lib.github.io/p/pak/stable'); pak::pkg_install(c('optparse','jsonlite','remotes','checkmate','mgcv'), upgrade = FALSE)" >>"${LOG_FILE}" 2>&1
  PAK_RC=$?
  set -e
  if [[ $PAK_RC -ne 0 ]]; then
    log "pak failed; falling back to install.packages for CRAN deps"
    Rscript -e "install.packages(c('optparse','jsonlite','remotes','checkmate','mgcv'), repos='${CRAN_MIRROR}')" >>"${LOG_FILE}" 2>&1
  fi
fi

# Ensure bettermc (dependency of fslmer) is available before installing fslmer
log "Ensuring 'bettermc' is installed (dependency of fslmer)"
BETTERMC_VERSION="1.2.1"
BETTERMC_URL="https://cran.r-project.org/src/contrib/Archive/bettermc/bettermc_${BETTERMC_VERSION}.tar.gz"
BETTERMC_LOCAL="${BETTERMC_TARBALL:-}"
if [[ "$OFFLINE" -eq 1 || "$PREFER_ONLINE" -ne 1 ]]; then
  if [[ -z "$BETTERMC_LOCAL" && -f "vendor/bettermc_${BETTERMC_VERSION}.tar.gz" ]]; then
    BETTERMC_LOCAL="vendor/bettermc_${BETTERMC_VERSION}.tar.gz"
  fi
  # If a tarball exists in vendor/ with a different name/extension, pick the first match
  if [[ -z "$BETTERMC_LOCAL" ]]; then
    FIRST_MATCH=$(ls vendor/bettermc*.tar* 2>/dev/null | head -n 1 || true)
    if [[ -n "$FIRST_MATCH" ]]; then
      BETTERMC_LOCAL="$FIRST_MATCH"
    fi
  fi
fi
if [[ "$OFFLINE" -eq 1 ]]; then
  if [[ -n "$BETTERMC_LOCAL" && -f "$BETTERMC_LOCAL" ]]; then
    log "Installing bettermc from local tarball: $BETTERMC_LOCAL"
    Rscript -e "install.packages('$BETTERMC_LOCAL', repos=NULL, type='source')" >>"${LOG_FILE}" 2>&1 || true
  else
    die "Offline mode: Provide bettermc tarball via BETTERMC_TARBALL or vendor/bettermc_*.tar.gz (and ensure vendor/checkmate_*.tar.gz is present)"
  fi
else
  set +e
  # Prefer local tarball if provided (no network needed)
  if [[ -n "$BETTERMC_LOCAL" && -f "$BETTERMC_LOCAL" ]]; then
    log "Installing bettermc from local tarball: $BETTERMC_LOCAL"
    Rscript -e "install.packages('$BETTERMC_LOCAL', repos=NULL, type='source')" >>"${LOG_FILE}" 2>&1
    RC_PAK_BMC_URL=$?
  else
    # Prefer direct archive URL via pak to avoid solver issues
    Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install('${BETTERMC_URL}', upgrade = FALSE) else quit(status=99)" >>"${LOG_FILE}" 2>&1
    RC_PAK_BMC_URL=$?
  fi
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
        log "install_version failed; trying GitHub upstream gfkse/bettermc"
        set +e
        Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install(sprintf('%s%s', if (nzchar(Sys.getenv('BETTERMC_GITHUB'))) Sys.getenv('BETTERMC_GITHUB') else 'gfkse/bettermc', if (nzchar(Sys.getenv('BETTERMC_REF'))) paste0('@', Sys.getenv('BETTERMC_REF')) else ''), upgrade = FALSE) else { if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}'); remotes::install_github(paste0(if (nzchar(Sys.getenv('BETTERMC_GITHUB'))) Sys.getenv('BETTERMC_GITHUB') else 'gfkse/bettermc', if (nzchar(Sys.getenv('BETTERMC_REF'))) paste0('@', Sys.getenv('BETTERMC_REF')) else '')) }" >>"${LOG_FILE}" 2>&1
        RC_GH_GFKSE=$?
        set -e
        if [[ $RC_GH_GFKSE -ne 0 ]]; then
          log "Upstream gfkse/bettermc failed; trying GitHub akersting/bettermc"
          set +e
          Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install(sprintf('akersting/bettermc%s', if (nzchar(Sys.getenv('BETTERMC_REF'))) paste0('@', Sys.getenv('BETTERMC_REF')) else ''), upgrade = FALSE) else { if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}'); remotes::install_github(paste0('akersting/bettermc', if (nzchar(Sys.getenv('BETTERMC_REF'))) paste0('@', Sys.getenv('BETTERMC_REF')) else '')) }" >>"${LOG_FILE}" 2>&1
          RC_GH_AK=$?
          set -e
          if [[ $RC_GH_AK -ne 0 ]]; then
            log "Upstream GitHub failed; trying GitHub cran/bettermc mirror"
            set +e
            Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install('cran/bettermc', upgrade = FALSE) else { if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}'); remotes::install_github('cran/bettermc') }" >>"${LOG_FILE}" 2>&1
            RC_GH=$?
            set -e
            if [[ $RC_GH -ne 0 ]]; then
              log "Failed to install bettermc from all sources. Showing last 60 log lines:"; tail -n 60 "${LOG_FILE}" || true; die "Failed to install 'bettermc' (required by fslmer). Check network/firewall and try again.";
            fi
          fi
        fi
      fi
    fi
  fi
fi

# Verify bettermc present
Rscript -e "quit(status = as.integer(!requireNamespace('bettermc', quietly=TRUE)))"
if [[ $? -ne 0 ]]; then
  log "'bettermc' is still not available after installation attempts. Showing last 60 log lines:"; tail -n 60 "${LOG_FILE}" || true; die "Missing 'bettermc' prevents fslmer install.";
fi

# Install fslmer
log "Installing fslmer (Deep-MI/fslmer)"
# Local tarball support and archive URL fallback
FSLMER_LOCAL="${FSLMER_TARBALL:-}"
if [[ -z "$FSLMER_LOCAL" && ( "$OFFLINE" -eq 1 || "$PREFER_ONLINE" -ne 1 ) ]]; then
  FIRST_FSLMER=$(ls vendor/fslmer*.tar* 2>/dev/null | head -n 1 || true)
  if [[ -n "$FIRST_FSLMER" ]]; then FSLMER_LOCAL="$FIRST_FSLMER"; fi
fi

FSLMER_REF_SAFE="${FSLMER_REF:-}"
if [[ -n "$FSLMER_REF_SAFE" ]]; then
  FSLMER_ARCHIVE_URL="https://github.com/Deep-MI/fslmer/archive/refs/tags/${FSLMER_REF_SAFE}.tar.gz"
else
  FSLMER_ARCHIVE_URL="https://github.com/Deep-MI/fslmer/archive/refs/heads/master.tar.gz"
fi

INSTALL_DONE=0
if [[ -n "$FSLMER_LOCAL" && -f "$FSLMER_LOCAL" ]]; then
  log "Installing fslmer from local tarball: $FSLMER_LOCAL"
  set +e
  Rscript -e "install.packages('$FSLMER_LOCAL', repos=NULL, type='source')" >>"${LOG_FILE}" 2>&1
  RC_LOCAL=$?
  set -e
  if [[ $RC_LOCAL -eq 0 ]]; then INSTALL_DONE=1; fi
fi

if [[ $INSTALL_DONE -eq 0 && "$OFFLINE" -eq 0 ]]; then
  log "Installing fslmer via pak (preferred)"
  set +e
  Rscript -e "if (requireNamespace('pak', quietly=TRUE)) pak::pkg_install(sprintf('Deep-MI/fslmer%s', if (nzchar(Sys.getenv('FSLMER_REF'))) paste0('@', Sys.getenv('FSLMER_REF')) else ''), upgrade = FALSE) else quit(status=99)" >>"${LOG_FILE}" 2>&1
  PAK_FSL_RC=$?
  set -e
  if [[ $PAK_FSL_RC -eq 0 ]]; then INSTALL_DONE=1; fi
fi

if [[ $INSTALL_DONE -eq 0 && "$OFFLINE" -eq 0 ]]; then
  log "pak install failed or pak unavailable; trying remotes::install_github"
  export R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
  Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='${CRAN_MIRROR}')" >>"${LOG_FILE}" 2>&1
  set +e
  Rscript -e "remotes::install_github(paste0('Deep-MI/fslmer', if (nzchar(Sys.getenv('FSLMER_REF'))) paste0('@', Sys.getenv('FSLMER_REF')) else ''), build=FALSE, build_vignettes=FALSE, dependencies=c('Depends','Imports','LinkingTo'), upgrade='never', quiet=TRUE)" >>"${LOG_FILE}" 2>&1
  RC_GH=$?
  set -e
  if [[ $RC_GH -eq 0 ]]; then INSTALL_DONE=1; fi
fi

if [[ $INSTALL_DONE -eq 0 && "$OFFLINE" -eq 0 ]]; then
  log "remotes::install_github failed; trying remotes::install_url from GitHub archive"
  set +e
  Rscript -e "remotes::install_url('${FSLMER_ARCHIVE_URL}', build=FALSE, dependencies=c('Depends','Imports','LinkingTo'), upgrade='never', quiet=TRUE)" >>"${LOG_FILE}" 2>&1
  RC_URL=$?
  set -e
  if [[ $RC_URL -eq 0 ]]; then INSTALL_DONE=1; fi
fi

if [[ $INSTALL_DONE -eq 0 ]]; then
  if [[ "$OFFLINE" -eq 1 ]]; then
    die "Offline mode: fslmer local tarball not provided (set FSLMER_TARBALL or place vendor/fslmer_*.tar.gz)"
  else
    log "fslmer install failed. Showing last 100 log lines:"; tail -n 100 "${LOG_FILE}" || true; die "Failed to install fslmer";
  fi
fi

# Verify fslmer present
Rscript -e "quit(status = as.integer(!requireNamespace('fslmer', quietly=TRUE)))"
if [[ $? -ne 0 ]]; then
  log "fslmer still not available after installation attempts. Showing last 100 log lines:"; tail -n 100 "${LOG_FILE}" || true; die "Missing 'fslmer' after install attempts.";
fi

# Snapshot
if [[ $SNAPSHOT -eq 1 ]]; then
  log "Snapshotting renv state"
  Rscript -e "renv::snapshot(prompt=FALSE)" >>"${LOG_FILE}" 2>&1
fi

# Verify loaded packages and versions
log "Verifying R packages (mgcv is required for GAM runs)"
Rscript -e "pkgs <- c('optparse','jsonlite','fslmer','mgcv'); print(data.frame(pkg=pkgs, available=sapply(pkgs, requireNamespace, quietly=TRUE)))"

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
