#!/usr/bin/env bash
set -euo pipefail

# install.sh (now under scripts/)
# Purpose: One-command installer for the FastSurfer R environment.
# - Verifies basic CLI prerequisites
# - Bootstraps micromamba and creates the fastsurfer-r env from scripts/environment.yml
# - Installs R extras (bettermc, fslmer) in the right order
# - Prints activation and smoke-test instructions
#
# Usage:
#   bash scripts/install.sh [--use-specs] [--no-compilers] [--pkgs-dir DIR] [--tmpdir DIR]
#                           [--auto-apt] [--env fastsurfer-r] [--r 4.5] [--no-fsqc]
#
# Notes:
# - --auto-apt will attempt to apt-get install missing CLI tools on Debian/Ubuntu with sudo

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

USE_SPECS=0
STRICT_YAML=1
NO_COMPILERS=0
PKGS_DIR=""
TMP_DIR=""
AUTO_APT=0
ENV_NAME="fastsurfer-r"
R_VERSION=""
NO_FSQC=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --use-specs) USE_SPECS=1; shift ;;
    --no-compilers) NO_COMPILERS=1; shift ;;
    --pkgs-dir) PKGS_DIR="${2:-}"; shift 2 ;;
    --tmpdir) TMP_DIR="${2:-}"; shift 2 ;;
    --auto-apt) AUTO_APT=1; shift ;;
    --env) ENV_NAME="${2:-}"; shift 2 ;;
    --r) R_VERSION="${2:-}"; shift 2 ;;
    --no-fsqc) NO_FSQC=1; shift ;;
    -h|--help)
      echo "Usage: bash scripts/install.sh [--use-specs] [--no-compilers] [--pkgs-dir DIR] [--tmpdir DIR] [--auto-apt] [--env NAME] [--r VERSION] [--no-fsqc]"
      echo
      echo "Options:"
      echo "  --use-specs       Use explicit conda specs (env lock) instead of environment.yml"
      echo "  --no-compilers    Skip adding C/C++/Fortran compilers"
      echo "  --pkgs-dir DIR    Point micromamba to an existing packages cache"
      echo "  --tmpdir DIR      Use a specific temp directory for downloads"
      echo "  --auto-apt        Attempt to apt-get missing shell tools (Linux only; requires sudo)"
      echo "  --env NAME        Micromamba environment name (default: fastsurfer-r)"
      echo "  --r VERSION       R version override (default from environment.yml)"
  echo "  --no-fsqc         Skip installing fsqc into the micromamba env (default: install)"
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

echo "[install] Preflight checks"
NEEDED=(curl tar grep awk sed head)
MISSING=()
for c in "${NEEDED[@]}"; do
  command -v "$c" >/dev/null 2>&1 || MISSING+=("$c")
done

OS="$(uname -s)"
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "[install] Missing tools: ${MISSING[*]}"
  if [[ "$AUTO_APT" -eq 1 && "$OS" == "Linux" && -x "/usr/bin/apt-get" ]]; then
    echo "[install] Attempting to install prerequisites via apt (requires sudo)"
    sudo apt-get update -y
    sudo apt-get install -y "${MISSING[@]}" bzip2 ca-certificates
  else
    echo "[install] Please install the missing tools and re-run. Suggestions:"
    if [[ "$OS" == "Darwin" ]]; then
      echo "  - macOS: xcode-select --install (for basic dev tools) or brew install ${MISSING[*]}"
    else
      echo "  - Ubuntu/Debian: sudo apt-get update && sudo apt-get install -y ${MISSING[*]} bzip2"
    fi
    exit 3
  fi
fi

# Default to strict YAML mode (use environment.yml as-is). Users can opt-out with --use-specs.

echo "[install] Creating/Updating micromamba env: $ENV_NAME"
CMD=(bash mamba_setup.sh --yaml environment.yml --strict-yaml)
[[ "$USE_SPECS" -eq 1 ]] && CMD+=(--use-specs)
[[ "$NO_COMPILERS" -eq 1 ]] && CMD+=(--no-compilers)
[[ -n "$PKGS_DIR" ]] && CMD+=(--pkgs-dir "$PKGS_DIR")
[[ -n "$TMP_DIR" ]] && CMD+=(--tmpdir "$TMP_DIR")
[[ -n "$ENV_NAME" ]] && CMD+=(--env "$ENV_NAME")
[[ -n "$R_VERSION" ]] && CMD+=(--r "$R_VERSION")

echo "[install] Running: ${CMD[*]}"
"${CMD[@]}"

echo
if [[ "$NO_FSQC" -eq 0 ]]; then
  echo "[install] Installing fsqc (Deep-MI) into the micromamba environment ($ENV_NAME)"
  # Use python/pip inside the activated env; the setup script printed activation guidance above.
  # We'll temporarily activate to install fsqc.
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/mamba_activate.sh"
  python -m pip install --upgrade pip wheel setuptools
  python -m pip install fsqc
  echo "[install] fsqc installed into $ENV_NAME. 'run_fsqc' should now be available after activation."
else
  echo "[install] Skipping fsqc installation (--no-fsqc). You can later install with:"
  echo "  source scripts/mamba_activate.sh && python -m pip install fsqc"
fi

echo
echo "[install] Done. Activate the environment with:"
echo "  source scripts/mamba_activate.sh"
echo "Then verify R in the env:"
echo "  which Rscript"
echo "  Rscript -e 'pkgs <- c(\"optparse\",\"jsonlite\",\"mgcv\",\"checkmate\",\"bettermc\",\"fslmer\"); print(sapply(pkgs, requireNamespace, quietly=TRUE))'"
echo "Optionally verify fsqc is available:"
echo "  run_fsqc --help"
echo "Quick help for the analysis script:"
echo "  Rscript scripts/fslmer_univariate.R --help"
echo "Deactivate later with:"
echo "  source scripts/mamba_deactivate.sh"
