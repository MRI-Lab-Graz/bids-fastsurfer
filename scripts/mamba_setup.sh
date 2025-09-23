#!/usr/bin/env bash
set -euo pipefail

# mamba_setup.sh
# Purpose: Bootstrap a local micromamba, create/update the 'fastsurfer-r' env from env/environment.yml,
# and install R packages not on conda-forge (bettermc specific version, Deep-MI/fslmer) inside the env.
#
# Usage:
#   bash scripts/mamba_setup.sh [--prefix <dir>] [--env fastsurfer-r] [--r 4.5]
#                                [--no-compilers] [--tmpdir /big/tmp] [--pkgs-dir /big/mamba-pkgs]
#
# Notes:
# - No sudo required. Installs micromamba under ~/.local/micromamba by default.
# - For low disk space: use --no-compilers to skip heavy toolchains; redirect caches with --pkgs-dir and temp with --tmpdir.
# - Adds a small activation note at the end. This script does NOT modify your shell rc files.

PREFIX="$HOME/.local/micromamba"
ENV_NAME="fastsurfer-r"
R_VERSION=""
NO_COMPILERS=0
USER_TMPDIR=""
USER_PKGS_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="${2:-}"; shift 2 ;;
    --env) ENV_NAME="${2:-}"; shift 2 ;;
    --r) R_VERSION="${2:-}"; shift 2 ;;
    --no-compilers) NO_COMPILERS=1; shift ;;
    --tmpdir) USER_TMPDIR="${2:-}"; shift 2 ;;
    --pkgs-dir) USER_PKGS_DIR="${2:-}"; shift 2 ;;
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
  # Ensure tar supports bzip2; fallback to decompression via bunzip2 if available
  if ! tar --help 2>&1 | grep -qi bzip2; then
    echo "[mamba_setup] Warning: tar may lack bzip2 support; attempting fallback path" >&2
  fi
  # Stream-extract directly to PREFIX; some builds place bin/micromamba, others micromamba at top-level
  if ! curl -fsSL "$URL" | tar -xj -C "$PREFIX" 2>/dev/null; then
    echo "[mamba_setup] Download/extract failed from $URL" >&2
    echo "[mamba_setup] Trying to download to file for inspection..." >&2
    TMP_ARCHIVE="$PREFIX/micromamba_${PLATFORM_TAG}.tar.bz2"
    if curl -fsSL "$URL" -o "$TMP_ARCHIVE"; then
      echo "[mamba_setup] Saved archive to $TMP_ARCHIVE (file size): $(wc -c < "$TMP_ARCHIVE" 2>/dev/null || echo unknown) bytes"
      echo "[mamba_setup] Retrying tar extraction with verbose output..."
      tar -xjvf "$TMP_ARCHIVE" -C "$PREFIX" || {
        echo "[mamba_setup] Unable to extract micromamba. Ensure 'bzip2' and 'tar' are available and that the platform tag is correct ($PLATFORM_TAG)." >&2
        exit 4
      }
    else
      echo "[mamba_setup] Failed to download micromamba from $URL" >&2
      exit 4
    fi
  fi
  # Relocate binary to PREFIX/bin if needed
  if [[ -x "$PREFIX/bin/micromamba" ]]; then
    : # already correct
  elif [[ -x "$PREFIX/micromamba" ]]; then
    mv "$PREFIX/micromamba" "$PREFIX/bin/micromamba"
  elif [[ -x "$PREFIX/bin/micromamba.exe" ]]; then
    mv "$PREFIX/bin/micromamba.exe" "$PREFIX/bin/micromamba" || true
  elif [[ -x "$PREFIX/micromamba/micromamba" ]]; then
    mv "$PREFIX/micromamba/micromamba" "$PREFIX/bin/micromamba"
  fi
fi

export MAMBA_ROOT_PREFIX="$PREFIX"
if [[ -n "$USER_PKGS_DIR" ]]; then
  mkdir -p "$USER_PKGS_DIR"
  export CONDA_PKGS_DIRS="$USER_PKGS_DIR"
  export MAMBA_PKGS_DIRS="$USER_PKGS_DIR"
  echo "[mamba_setup] Using custom pkgs dir: $USER_PKGS_DIR"
fi
if [[ -n "$USER_TMPDIR" ]]; then
  mkdir -p "$USER_TMPDIR"
  export TMPDIR="$USER_TMPDIR"
  echo "[mamba_setup] Using custom TMPDIR: $TMPDIR"
fi
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

# Optionally override R version; also sanitize YAML (strip comments and blank lines)
TMP_RAW_YAML="${YAML}.raw.tmp"
TMP_YAML="${YAML}.tmp"
if [[ -n "$R_VERSION" ]]; then
  echo "[mamba_setup] Using R version r-base=${R_VERSION}"
  awk -v ver="$R_VERSION" '{ if ($0 ~ /^\s*-\s*r-base/) { print "  - r-base=" ver } else { print } }' "$YAML" > "$TMP_RAW_YAML"
else
  cp "$YAML" "$TMP_RAW_YAML"
fi
# Optionally remove heavy compiler toolchains to save disk space
if [[ $NO_COMPILERS -eq 1 ]]; then
  echo "[mamba_setup] --no-compilers enabled: filtering compiler packages from YAML"
  awk '{
    line=$0
    if (line ~ /^\s*-\s*(c-compiler|cxx-compiler|fortran-compiler)\s*$/) next
    print line
  }' "$TMP_RAW_YAML" > "${TMP_RAW_YAML}.noc"
  mv "${TMP_RAW_YAML}.noc" "$TMP_RAW_YAML"
fi
# Strip comments and empty lines safely, but keep section headers
awk '{
  line=$0
  # Remove inline comments that start with # and have space before or after
  sub(/\s*#.*$/, "", line)
  # Trim trailing spaces
  sub(/[ \t]+$/, "", line)
  # Skip empty lines
  if (line ~ /^\s*$/) next
  print line
}' "$TMP_RAW_YAML" > "$TMP_YAML"

echo "[mamba_setup] Creating/updating env '$ENV_NAME' from $TMP_YAML"
# Try creating from YAML first
set +e
"$MAMBA_BIN" create -y -n "$ENV_NAME" -f "$TMP_YAML"
CREATE_RC=$?
set -e
if [[ $CREATE_RC -ne 0 ]]; then
  echo "[mamba_setup] YAML-based create failed (rc=$CREATE_RC). Showing sanitized YAML (first 50 lines):"
  head -n 50 "$TMP_YAML" || true
  echo "[mamba_setup] Falling back to spec-based create by extracting dependencies..."
  # Extract dependency specs from sanitized YAML
  mapfile -t SPECS < <(awk '/^dependencies:/ {in_dep=1; next} /^\w/ {in_dep=0} in_dep && /^\s*-\s*/ { sub(/^\s*-\s*/, ""); if (length($0)>0) print $0 }' "$TMP_YAML")
  if [[ $NO_COMPILERS -eq 1 ]]; then
    NEW_SPECS=()
    for s in "${SPECS[@]}"; do
      case "$s" in
        c-compiler|cxx-compiler|fortran-compiler) continue ;;
        *) NEW_SPECS+=("$s") ;;
      esac
    done
    SPECS=("${NEW_SPECS[@]}")
  fi
  if [[ ${#SPECS[@]} -eq 0 ]]; then
    echo "[mamba_setup] No dependency specs could be parsed from $TMP_YAML" >&2
    exit 5
  fi
  echo "[mamba_setup] Specs: ${SPECS[*]}"
  # Use explicit channel from YAML or default to conda-forge
  CHANNELS=("-c" "conda-forge")
  set +e
  "$MAMBA_BIN" create -y -n "$ENV_NAME" "${CHANNELS[@]}" "${SPECS[@]}"
  CREATE_RC=$?
  set -e
  if [[ $CREATE_RC -ne 0 ]]; then
    echo "[mamba_setup] Spec-based create also failed (rc=$CREATE_RC). Aborting." >&2
    exit $CREATE_RC
  fi
else
  # On success, try to update from YAML to ensure sync (harmless if no changes)
  set +e
  "$MAMBA_BIN" env update -n "$ENV_NAME" -f "$TMP_YAML" --prune
  set -e
fi

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
