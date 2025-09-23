#!/usr/bin/env bash
set -euo pipefail

# mamba_setup.sh
# Purpose: Bootstrap a local micromamba, create/update the 'fastsurfer-r' env from env/environment.yml,
# and install R packages not on conda-forge (bettermc specific version, Deep-MI/fslmer) inside the env.
#
# Usage:
#   bash scripts/mamba_setup.sh [--prefix <dir>] [--env fastsurfer-r] [--r 4.5]
#                                [--no-compilers] [--tmpdir /big/tmp] [--pkgs-dir /big/mamba-pkgs]
#                                [--skip-extras] [--bettermc-version 1.2.1]
#                                [--yaml /abs/path/to/environment.yml]
#                                [--print-yaml] [--print-specs] [--use-specs]
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
USER_PREFIX_SET=0
SKIP_EXTRAS=0
BETTERMC_VERSION="1.2.1"
YAML_OVERRIDE=""
PRINT_YAML=0
PRINT_SPECS=0
USE_SPECS=0
STRICT_YAML=0

while [[ $# -gt 0 ]]; do
  case "$1" in
  --prefix) PREFIX="${2:-}"; USER_PREFIX_SET=1; shift 2 ;;
    --env) ENV_NAME="${2:-}"; shift 2 ;;
    --r) R_VERSION="${2:-}"; shift 2 ;;
    --no-compilers) NO_COMPILERS=1; shift ;;
    --tmpdir) USER_TMPDIR="${2:-}"; shift 2 ;;
    --pkgs-dir) USER_PKGS_DIR="${2:-}"; shift 2 ;;
    --skip-extras) SKIP_EXTRAS=1; shift ;;
    --bettermc-version) BETTERMC_VERSION="${2:-}"; shift 2 ;;
    --yaml) YAML_OVERRIDE="${2:-}"; shift 2 ;;
    --print-yaml) PRINT_YAML=1; shift ;;
    --print-specs) PRINT_SPECS=1; shift ;;
    --use-specs) USE_SPECS=1; shift ;;
    --strict-yaml) STRICT_YAML=1; shift ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

# If pkgs-dir is on a large volume and no explicit --prefix was given, place the envs there too
if [[ -n "$USER_PKGS_DIR" && $USER_PREFIX_SET -eq 0 ]]; then
  PKG_PARENT="$(cd "$(dirname "$USER_PKGS_DIR")" && pwd)"
  PREFIX="$PKG_PARENT/micromamba-root"
  echo "[mamba_setup] No --prefix provided; setting micromamba root to: $PREFIX"
fi

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
echo "[mamba_setup] Effective paths:" 
echo "  MAMBA_ROOT_PREFIX=$MAMBA_ROOT_PREFIX"
echo "  CONDA_PKGS_DIRS=${CONDA_PKGS_DIRS:-}"
echo "  MAMBA_PKGS_DIRS=${MAMBA_PKGS_DIRS:-}"
echo "  TMPDIR=${TMPDIR:-}"

# Prefer an environment.yml located alongside this script; fallback to env/environment.yml
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YAML="$SCRIPT_DIR/environment.yml"
FALLBACK_YAML="env/environment.yml"

# Allow explicit override
if [[ -n "$YAML_OVERRIDE" ]]; then
  if [[ ! -f "$YAML_OVERRIDE" ]]; then
    echo "[mamba_setup] --yaml provided but file not found: $YAML_OVERRIDE" >&2
    exit 2
  fi
  YAML="$YAML_OVERRIDE"
elif [[ -f "$FALLBACK_YAML" && ! -f "$YAML" ]]; then
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
  - c-compiler
  - cxx-compiler
  - fortran-compiler
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

# Build effective YAML path
STRICT_PASSTHRU=0
if [[ $STRICT_YAML -eq 1 && -z "$R_VERSION" ]]; then
  # Use the provided YAML as-is without sanitization
  TMP_YAML="$YAML"
  STRICT_PASSTHRU=1
else
  # Optionally override R version; also sanitize YAML (strip comments and blank lines)
  # Use temp files to avoid stale/partial repo files between runs
  _TMPDIR_FOR_YAML="${USER_TMPDIR:-${TMPDIR:-/tmp}}"
  mkdir -p "$_TMPDIR_FOR_YAML" 2>/dev/null || true
  TMP_RAW_YAML="$(mktemp -t env_raw_XXXXXXXX.yml 2>/dev/null || mktemp "$_TMPDIR_FOR_YAML/env_raw_XXXXXXXX.yml")"
  TMP_YAML="$(mktemp -t env_san_XXXXXXXX.yml 2>/dev/null || mktemp "$_TMPDIR_FOR_YAML/env_san_XXXXXXXX.yml")"
  trap 'rm -f "$TMP_RAW_YAML" "$TMP_YAML" 2>/dev/null || true' EXIT
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
    sub(/\s*#.*$/, "", line)
    sub(/[ \t]+$/, "", line)
    if (line ~ /^\s*$/) next
    print line
  }' "$TMP_RAW_YAML" > "$TMP_YAML"
  # Ensure conda-forge channel exists
  if ! grep -qE '^[[:space:]]*-[[:space:]]*conda-forge[[:space:]]*$' "$TMP_YAML"; then
    echo "[mamba_setup] Adding missing conda-forge channel to YAML"
    printf "channels:\n  - conda-forge\n%s" "$(cat "$TMP_YAML")" > "$TMP_YAML"
  fi
  # Guard: dependencies block must exist
  if ! grep -qE '^dependencies:' "$TMP_YAML"; then
    echo "[mamba_setup] Error: YAML has no dependencies block: $YAML" >&2
    exit 5
  fi
fi

# If user only wants to inspect the effective YAML/specs, print and exit
if [[ $PRINT_YAML -eq 1 ]]; then
  echo "[mamba_setup] Effective sanitized YAML ($TMP_YAML):"
  cat "$TMP_YAML"
  exit 0
fi

if [[ $USE_SPECS -eq 0 ]]; then
  echo "[mamba_setup] Creating/updating env '$ENV_NAME' from $TMP_YAML"
  # Try creating from YAML first
  set +e
  "$MAMBA_BIN" create -y -n "$ENV_NAME" -f "$TMP_YAML"
  CREATE_RC=$?
  set -e
else
  # Force spec-based path
  CREATE_RC=1
fi
if [[ $CREATE_RC -ne 0 ]]; then
  echo "[mamba_setup] YAML-based create failed (rc=$CREATE_RC). Showing sanitized YAML (first 50 lines):"
  head -n 50 "$TMP_YAML" || true
  if [[ $STRICT_YAML -eq 1 ]]; then
    echo "[mamba_setup] STRICT YAML mode: not falling back to specs. Fix the YAML and retry." >&2
    exit $CREATE_RC
  fi
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
  # Ensure required R packages are present even if YAML is minimal
  REQUIRED=(make r-optparse r-jsonlite r-mgcv r-remotes r-bh r-matrix r-rcpp r-rcpparmadillo)
  for req in "${REQUIRED[@]}"; do
    present=0
    for s in "${SPECS[@]}"; do
      if [[ "$s" == "$req"* ]]; then present=1; break; fi
    done
    if [[ $present -eq 0 ]]; then SPECS+=("$req"); fi
  done
  if [[ ${#SPECS[@]} -eq 0 ]]; then
    echo "[mamba_setup] No dependency specs parsed; synthesizing REQUIRED specs"
    SPECS=(make r-optparse r-jsonlite r-mgcv r-remotes r-bh r-matrix r-rcpp r-rcpparmadillo)
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

if [[ $SKIP_EXTRAS -eq 1 ]]; then
  echo "[mamba_setup] --skip-extras set: skipping bettermc and fslmer installs"
else
  echo "[mamba_setup] Installing R extras inside env: bettermc (robust), Deep-MI/fslmer"
  "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")'

  # Ensure prerequisite CRAN libs are present via conda-forge when possible (binary installs)
  echo "[mamba_setup] Ensuring base R deps via conda-forge: r-checkmate r-backports"
  set +e
  "$MAMBA_BIN" install -y -n "$ENV_NAME" -c conda-forge r-checkmate r-backports >/dev/null 2>&1
  RC_CF_DEPS=$?
  set -e
  if [[ $RC_CF_DEPS -ne 0 ]]; then
    echo "[mamba_setup] conda-forge install failed; falling back to install.packages for checkmate/backports"
    "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e 'if (!requireNamespace("backports", quietly=TRUE)) install.packages("backports", repos="https://cloud.r-project.org"); if (!requireNamespace("checkmate", quietly=TRUE)) install.packages("checkmate", repos="https://cloud.r-project.org")'
  fi

  # Helper to download bettermc tarball to vendor dir
  VENDOR_DIR="$SCRIPT_DIR/vendor"; mkdir -p "$VENDOR_DIR"
  BETTERMC_TARBALL="${BETTERMC_TARBALL:-}"
  if [[ -z "$BETTERMC_TARBALL" ]]; then
    CRAN_URL="https://cran.r-project.org/src/contrib/Archive/bettermc/bettermc_${BETTERMC_VERSION}.tar.gz"
    GH_URL="https://github.com/akersting/bettermc/releases/download/v${BETTERMC_VERSION}/bettermc_${BETTERMC_VERSION}.tar.gz"
    TGT="$VENDOR_DIR/bettermc_${BETTERMC_VERSION}.tar.gz"
    echo "[mamba_setup] Attempting to download bettermc ${BETTERMC_VERSION} tarball to $TGT"
    set +e
    curl -fL "$CRAN_URL" -o "$TGT"
    DL_RC=$?
    if [[ $DL_RC -ne 0 ]]; then
      echo "[mamba_setup] CRAN download failed (rc=$DL_RC), trying GitHub release asset"
      curl -fL "$GH_URL" -o "$TGT"
      DL_RC=$?
    fi
    set -e
    if [[ $DL_RC -eq 0 ]]; then
      BETTERMC_TARBALL="$TGT"
      echo "[mamba_setup] Downloaded bettermc tarball: $BETTERMC_TARBALL"
    else
      echo "[mamba_setup] Warning: could not download bettermc tarball from CRAN or GitHub. Will try remotes next."
    fi
  fi

  # Prefer local tarball install
  if [[ -n "$BETTERMC_TARBALL" && -f "$BETTERMC_TARBALL" ]]; then
    echo "[mamba_setup] Installing bettermc from tarball: $BETTERMC_TARBALL"
    set +e
    "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e "install.packages('${BETTERMC_TARBALL//"'"/"'\\''"}', repos=NULL, type='source')"
    BM_RC=$?
    set -e
  else
    BM_RC=1
  fi

  # Fallback to remotes if needed
  if [[ $BM_RC -ne 0 ]]; then
    echo "[mamba_setup] Installing bettermc via remotes (archive or GitHub)"
    set +e
    "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e "ok <- TRUE; tryCatch({ if (!requireNamespace('bettermc', quietly=TRUE)) remotes::install_version('bettermc', version='${BETTERMC_VERSION}', repos='https://cran.r-project.org') }, error=function(e) ok<<-FALSE); quit(status = as.integer(!ok))"
    RC=$?
    set -e
    if [[ $RC -ne 0 ]]; then
      echo "[mamba_setup] bettermc install via CRAN version failed; trying GitHub mirrors"
      set +e
      "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e "ok <- FALSE; for (repo in c('cran/bettermc','gfkse/bettermc','akersting/bettermc')) { try({ remotes::install_github(repo, quiet=TRUE); if (requireNamespace('bettermc', quietly=TRUE)) { ok <- TRUE; break } }, silent=TRUE) }; quit(status = as.integer(!ok))"
      RC2=$?
      set -e
      if [[ $RC2 -ne 0 ]]; then
        echo "[mamba_setup] Warning: bettermc installation failed. You can set BETTERMC_TARBALL to a local .tar.gz and re-run. Continuing without bettermc."
      fi
    fi
  fi

  # Verify bettermc is available before attempting fslmer (hard dependency)
  set +e
  "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e 'quit(status = as.integer(!requireNamespace("bettermc", quietly=TRUE)))'
  BM_OK=$?
  set -e
  if [[ $BM_OK -ne 0 ]]; then
    echo "[mamba_setup] Error: 'bettermc' is not installed; fslmer depends on it. Skipping fslmer installation."
    echo "[mamba_setup] Hint: ensure checkmate/backports are present, provide BETTERMC_TARBALL, or allow internet for CRAN installs."
  else
    # fslmer: try local tarball first, else GitHub
    FSLMER_TARBALL="${FSLMER_TARBALL:-}"
    if [[ -n "$FSLMER_TARBALL" && -f "$FSLMER_TARBALL" ]]; then
      echo "[mamba_setup] Installing fslmer from local tarball: $FSLMER_TARBALL"
      set +e
      "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e "install.packages('${FSLMER_TARBALL//"'"/"'\\''"}', repos=NULL, type='source')"
      set -e
    else
      echo "[mamba_setup] Installing fslmer from GitHub Deep-MI/fslmer"
      set +e
      "$MAMBA_BIN" run -n "$ENV_NAME" Rscript --vanilla -e 'remotes::install_github("Deep-MI/fslmer", upgrade="never", quiet=TRUE)'
      FSL_RC=$?
      set -e
      if [[ $FSL_RC -ne 0 ]]; then
        echo "[mamba_setup] Warning: fslmer installation failed. The R helper can still run with --engine glm/gam."
      fi
    fi
  fi
fi

echo "[mamba_setup] Done. Activate with:"
# Persist activation info for convenience
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
{
  echo "MAMBA_ROOT_PREFIX=\"$MAMBA_ROOT_PREFIX\""
  echo "ENV_NAME=\"$ENV_NAME\""
} > "$ENV_RECORD"

echo "  source scripts/mamba_activate.sh   # activates like Python venv"
echo "Then run inside the env:"
echo "  Rscript scripts/fslmer_univariate.R --config /abs/path/to/config.json"
