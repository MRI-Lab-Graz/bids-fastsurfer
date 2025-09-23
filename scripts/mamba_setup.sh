#!/usr/bin/env bash
set -euo pipefail

# mamba_setup.sh
# Purpose: Bootstrap a local micromamba, create/update the 'fastsurfer-r' env from env/environment.yml,
# and install R packages not on conda-forge (bettermc specific version, Deep-MI/fslmer) inside the env.
#
# Usage:
#   bash scripts/mamba_setup.sh [--prefix <dir>] [--env fastsurfer-r] [--r 4.5]
#
# Notes:
# - No sudo required. Installs micromamba under ~/.local/micromamba by default.
# - Adds a small activation note at the end. This script does NOT modify your shell rc files.

PREFIX="$HOME/.local/micromamba"
ENV_NAME="fastsurfer-r"
R_VERSION=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="${2:-}"; shift 2 ;;
    --env) ENV_NAME="${2:-}"; shift 2 ;;
    --r) R_VERSION="${2:-}"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

mkdir -p "$PREFIX"

# Determine platform triplet expected by micromamba API
UNAME_S=$(uname -s)
UNAME_M=$(uname -m)
case "$UNAME_S" in
  Linux) OS="linux" ;;
  Darwin) OS="osx" ;;
  *) OS="linux" ;;
esac
case "$UNAME_M" in
  x86_64|amd64) ARCH_TAG="64" ;;
  aarch64|arm64)
    if [[ "$OS" == "osx" ]]; then ARCH_TAG="arm64"; else ARCH_TAG="aarch64"; fi ;;
  *) ARCH_TAG="64" ;;
esac
PLATFORM_TAG="$OS-$ARCH_TAG"

MAMBA_BIN="$PREFIX/bin/micromamba"
if [[ ! -x "$MAMBA_BIN" ]]; then
  echo "[mamba_setup] Installing micromamba ($PLATFORM_TAG) to $PREFIX"
  URL="https://micro.mamba.pm/api/micromamba/${PLATFORM_TAG}/latest"
  mkdir -p "$PREFIX/bin"
  # Stream-extract; requires bzip2 support in tar
  if ! curl -fsSL "$URL" | tar -xj -C "$PREFIX" --strip-components=1 bin/micromamba 2>/dev/null; then
    echo "[mamba_setup] Download/extract failed from $URL" >&2
    echo "[mamba_setup] Trying to download to file for inspection..." >&2
    TMP_ARCHIVE="$PREFIX/micromamba_${PLATFORM_TAG}.tar.bz2"
    if curl -fsSL "$URL" -o "$TMP_ARCHIVE"; then
      echo "[mamba_setup] Saved archive to $TMP_ARCHIVE (first 3 lines):"
      head -n 3 "$TMP_ARCHIVE" 2>/dev/null || true
      echo "[mamba_setup] Retrying tar extraction with verbose output..."
      tar -xjvf "$TMP_ARCHIVE" -C "$PREFIX" --strip-components=1 bin/micromamba || {
        echo "[mamba_setup] Unable to extract micromamba. Ensure 'bzip2' and 'tar' are available and that the platform tag is correct ($PLATFORM_TAG)." >&2
        exit 4
      }
    else
      echo "[mamba_setup] Failed to download micromamba from $URL" >&2
      exit 4
    fi
  fi
fi

export MAMBA_ROOT_PREFIX="$PREFIX"
"$MAMBA_BIN" shell hook -s bash >/dev/null 2>&1 || true

# Prefer an environment.yml located alongside this script; fallback to env/environment.yml
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YAML="$SCRIPT_DIR/environment.yml"
FALLBACK_YAML="env/environment.yml"

if [[ -f "$FALLBACK_YAML" && ! -f "$YAML" ]]; then
  echo "[mamba_setup] Using fallback $FALLBACK_YAML (no scripts/environment.yml found)."
  YAML="$FALLBACK_YAML"
fi

if [[ ! -f "$YAML" ]]; then
  echo "[mamba_setup] No environment.yml found — creating scripts/environment.yml with defaults."
  cat >"$SCRIPT_DIR/environment.yml" <<'YML'
name: fastsurfer-r
channels:
  - conda-forge
dependencies:
  - r-base=4.5
  - make
  - compilers
  - r-optparse
  - r-jsonlite
  - r-mgcv
  - r-remotes
  - r-bh
  - r-matrix
  - r-rcpp
  - r-rcpparmadillo
YML
  YAML="$SCRIPT_DIR/environment.yml"
fi

# Optionally override R version by editing a temp YAML on the fly
TMP_YAML="${YAML}.tmp"
if [[ -n "$R_VERSION" ]]; then
  echo "[mamba_setup] Using R version r-base=${R_VERSION}"
  # Replace any line starting with '  - r-base' with the chosen version
  awk -v ver="$R_VERSION" '{ if ($0 ~ /^\s*-\s*r-base/) { print "  - r-base=" ver } else { print } }' "$YAML" > "$TMP_YAML"
else
  cp "$YAML" "$TMP_YAML"
fi

echo "[mamba_setup] Creating/updating env from $TMP_YAML"
"$MAMBA_BIN" create -y -f "$TMP_YAML" || "$MAMBA_BIN" env update -f "$TMP_YAML"

echo "[mamba_setup] Installing R packages inside env: bettermc (robust), Deep-MI/fslmer"
"$MAMBA_BIN" run -n "$ENV_NAME" Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")'

# Try local tarball if provided, else CRAN archive version, else GitHub mirrors
BETTERMC_TARBALL="${BETTERMC_TARBALL:-}"
if [[ -n "$BETTERMC_TARBALL" && -f "$BETTERMC_TARBALL" ]]; then
  echo "[mamba_setup] Installing bettermc from local tarball: $BETTERMC_TARBALL"
  "$MAMBA_BIN" run -n "$ENV_NAME" Rscript -e "install.packages('${BETTERMC_TARBALL//"'"/"'\''"}', repos=NULL, type='source')"
else
  echo "[mamba_setup] Installing bettermc via remotes (archive/version or GitHub mirrors)"
  set +e
  "$MAMBA_BIN" run -n "$ENV_NAME" Rscript -e 'ok <- TRUE; tryCatch({ if (!requireNamespace("bettermc", quietly=TRUE)) remotes::install_version("bettermc", version="1.2.1", repos="https://cran.r-project.org") }, error=function(e) ok<<-FALSE); quit(status = as.integer(!ok))'
  RC=$?
  set -e
  if [[ $RC -ne 0 ]]; then
    echo "[mamba_setup] bettermc archived version failed; trying GitHub mirrors (cran/bettermc, gfkse/bettermc, akersting/bettermc)"
    set +e
    "$MAMBA_BIN" run -n "$ENV_NAME" Rscript -e 'ok <- FALSE; for (repo in c("cran/bettermc","gfkse/bettermc","akersting/bettermc")) { try({ remotes::install_github(repo, quiet=TRUE); if (requireNamespace("bettermc", quietly=TRUE)) { ok <- TRUE; break } }, silent=TRUE) }; quit(status = as.integer(!ok))'
    RC2=$?
    set -e
    if [[ $RC2 -ne 0 ]]; then
      echo "[mamba_setup] Failed to install bettermc automatically. You can set BETTERMC_TARBALL to a local .tar.gz and re-run."
      exit 1
    fi
  fi
fi

# fslmer: try local tarball first, else GitHub
FSLMER_TARBALL="${FSLMER_TARBALL:-}"
if [[ -n "$FSLMER_TARBALL" && -f "$FSLMER_TARBALL" ]]; then
  echo "[mamba_setup] Installing fslmer from local tarball: $FSLMER_TARBALL"
  "$MAMBA_BIN" run -n "$ENV_NAME" Rscript -e "install.packages('${FSLMER_TARBALL//"'"/"'\''"}', repos=NULL, type='source')"
else
  echo "[mamba_setup] Installing fslmer from GitHub Deep-MI/fslmer"
  "$MAMBA_BIN" run -n "$ENV_NAME" Rscript -e 'remotes::install_github("Deep-MI/fslmer", upgrade="never", quiet=TRUE)'
fi

echo "[mamba_setup] Done. Activate with:"
echo "  eval \"$($MAMBA_BIN shell hook -s bash)\""
echo "  micromamba activate $ENV_NAME"
echo "Then run:"
echo "  Rscript scripts/fslmer_univariate.R --config /abs/path/to/config.json"
