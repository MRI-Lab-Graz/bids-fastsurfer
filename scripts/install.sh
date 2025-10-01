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
#                           [--auto-apt] [--env fastsurfer-r] [--r 4.5]
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
# Precheck controls
MIN_TMP_GB=3
MIN_CACHE_GB=5
REQUIRE_FS=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --use-specs) USE_SPECS=1; shift ;;
    --no-compilers) NO_COMPILERS=1; shift ;;
    --pkgs-dir) PKGS_DIR="${2:-}"; shift 2 ;;
    --tmpdir) TMP_DIR="${2:-}"; shift 2 ;;
    --auto-apt) AUTO_APT=1; shift ;;
    --env) ENV_NAME="${2:-}"; shift 2 ;;
    --r) R_VERSION="${2:-}"; shift 2 ;;
    --min-tmp-gb) MIN_TMP_GB="${2:-3}"; shift 2 ;;
    --min-cache-gb) MIN_CACHE_GB="${2:-5}"; shift 2 ;;
  --require-fs) REQUIRE_FS=1; shift ;;
  --allow-no-fs) REQUIRE_FS=0; shift ;;
    -h|--help)
      echo "FastSurfer R environment installer"
      echo
      echo "Basic usage"
      echo "  # Standard install (requires FreeSurfer present)"
      echo "  bash scripts/install.sh"
      echo
      echo "  # Low-storage install: place caches on a large disk and skip compilers"
      echo "  bash scripts/install.sh \\
    --no-compilers \\
    --pkgs-dir /path/to/bigdisk/.mamba-cache \\
    --tmpdir   /path/to/bigdisk/tmp"
      echo
      echo "Expert usage"
      echo "  # Custom env name and R version"
      echo "  bash scripts/install.sh --env my-r-env --r 4.5"
      echo
      echo "  # Proceed without FreeSurfer (advanced; features depending on FS won’t work)"
      echo "  bash scripts/install.sh --allow-no-fs"
      echo
      echo "  # Auto-install missing CLI prerequisites (Debian/Ubuntu)"
      echo "  bash scripts/install.sh --auto-apt"
      echo
      echo "Options"
      echo "  --env NAME         Micromamba environment name (default: fastsurfer-r)"
      echo "  --r VERSION        R version override (default from environment.yml)"
      echo "  --no-compilers     Skip C/C++/Fortran compilers to save space"
      echo "  --pkgs-dir DIR     Package cache directory (use a large, writable path)"
      echo "  --tmpdir DIR       Temporary directory for downloads (use a large, writable path)"
      echo "  --min-tmp-gb N     Minimum free space in TMPDIR in GB (default: 3)"
      echo "  --min-cache-gb N   Minimum free space in package cache in GB (default: 5)"
      echo "  --require-fs       Fail if FreeSurfer is not detected locally (default)"
      echo "  --allow-no-fs      Proceed without FreeSurfer (advanced; not recommended)"
      echo "  --auto-apt         apt-get missing tools (Linux Debian/Ubuntu; requires sudo)"
      echo "  --use-specs        Use explicit conda specs instead of environment.yml"
      echo
      echo "Examples"
      echo "  # Require FreeSurfer (default) and set stricter space thresholds"
      echo "  bash scripts/install.sh --min-tmp-gb 5 --min-cache-gb 8"
      exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

# Friendly header
echo ""
echo "========================================"
echo "🧠  MRI Lab Graz  •  Karl Koschutnig"
echo "🚀  FastSurfer R Environment Installer"
echo "========================================"
echo ""

echo "🧪 Preflight checks"
NEEDED=(curl tar grep awk sed head)
MISSING=()
for c in "${NEEDED[@]}"; do
  command -v "$c" >/dev/null 2>&1 || MISSING+=("$c")
done

OS="$(uname -s)"
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "❗ Missing tools: ${MISSING[*]}"
  if [[ "$AUTO_APT" -eq 1 && "$OS" == "Linux" && -x "/usr/bin/apt-get" ]]; then
    echo "🛠️  Attempting to install prerequisites via apt (requires sudo)"
    sudo apt-get update -y
    sudo apt-get install -y "${MISSING[@]}" bzip2 ca-certificates
  else
    echo "ℹ️  Please install the missing tools and re-run. Suggestions:"
    if [[ "$OS" == "Darwin" ]]; then
      echo "  - macOS: xcode-select --install (for basic dev tools) or brew install ${MISSING[*]}"
    else
      echo "  - Ubuntu/Debian: sudo apt-get update && sudo apt-get install -y ${MISSING[*]} bzip2"
    fi
    exit 3
  fi
fi

# Disk space prechecks for tmp and package cache
check_space_gb() {
  local path="$1"; local need_gb="$2"; local label="$3"
  local probe="$path"
  if [[ ! -e "$probe" ]]; then
    probe="$(dirname "$probe" 2>/dev/null || echo /)"
  fi
  # df -Pk: POSIX portable, KB units
  local avail_kb
  avail_kb=$(df -Pk "$probe" 2>/dev/null | awk 'NR==2 {print $4}')
  if [[ -z "$avail_kb" ]]; then
    echo "⚠️  Could not determine free space for $label at $path"
    return 0
  fi
  local avail_gb
  avail_gb=$(( avail_kb / 1024 / 1024 ))
  if (( avail_gb < need_gb )); then
    echo "❌ Insufficient free space for $label at $path: ${avail_gb}GB available, need ≥ ${need_gb}GB" >&2
    return 1
  fi
  echo "✅ $label: ${avail_gb}GB free at $path (min ${need_gb}GB)"
  return 0
}

# Determine effective TMP and cache paths for checking
EFFECTIVE_TMP="${TMP_DIR:-${TMPDIR:-/tmp}}"
EFFECTIVE_CACHE="${PKGS_DIR:-$HOME/.local/micromamba}"

if ! check_space_gb "$EFFECTIVE_TMP" "$MIN_TMP_GB" "TMPDIR"; then
  echo "💡 Hint: re-run with --tmpdir pointing to a larger volume (and ensure write permissions)." >&2
  exit 3
fi
if ! check_space_gb "$EFFECTIVE_CACHE" "$MIN_CACHE_GB" "package cache"; then
  echo "💡 Hint: re-run with --pkgs-dir pointing to a larger volume, or free up space." >&2
  exit 3
fi

# FreeSurfer presence check (warn by default, or fail with --require-fs)
FS_HOME="${FREESURFER_HOME:-/usr/local/freesurfer}"
FS_SETUP="$FS_HOME/SetUpFreeSurfer.sh"
FS_LICENSE_PATH="${FS_LICENSE:-$FS_HOME/license.txt}"
if [[ -f "$FS_SETUP" ]]; then
  echo "🧠 Found FreeSurfer: $FS_HOME"
else
  if [[ "$REQUIRE_FS" -eq 1 ]]; then
    echo "❌ FreeSurfer not found. Please set up FreeSurfer before installation." >&2
    echo "  - Install FreeSurfer locally and locate SetUpFreeSurfer.sh (e.g., /usr/local/freesurfer)" >&2
    echo "  - Export FREESURFER_HOME to that path, e.g.:" >&2
    echo "      export FREESURFER_HOME=/usr/local/freesurfer" >&2
    echo "      source \"$FREESURFER_HOME/SetUpFreeSurfer.sh\"" >&2
    echo "  - Ensure your license is available at \$FS_LICENSE or \$FREESURFER_HOME/license.txt" >&2
    echo "If you intentionally want to proceed without FreeSurfer, re-run with --allow-no-fs." >&2
    exit 4
  fi
  echo "⚠️  FreeSurfer not detected. Set FREESURFER_HOME (or install at /usr/local/freesurfer)." >&2
fi
if [[ -f "$FS_LICENSE_PATH" ]]; then
  echo "🔑 License: $FS_LICENSE_PATH"
else
  echo "ℹ️  FreeSurfer license not found (expected at FS_LICENSE or $FS_HOME/license.txt)." >&2
fi

# Default to strict YAML mode (use environment.yml as-is). Users can opt-out with --use-specs.

echo "📦 Environment: $ENV_NAME"
CMD=(bash mamba_setup.sh --yaml environment.yml)
if [[ "$NO_COMPILERS" -eq 0 ]]; then
  CMD+=(--strict-yaml)
else
  echo "[install] --no-compilers requested: disabling strict YAML so compilers can be filtered"
fi
[[ "$USE_SPECS" -eq 1 ]] && CMD+=(--use-specs)
[[ "$NO_COMPILERS" -eq 1 ]] && CMD+=(--no-compilers)
[[ -n "$PKGS_DIR" ]] && CMD+=(--pkgs-dir "$PKGS_DIR")
[[ -n "$TMP_DIR" ]] && CMD+=(--tmpdir "$TMP_DIR")
[[ -n "$ENV_NAME" ]] && CMD+=(--env "$ENV_NAME")
[[ -n "$R_VERSION" ]] && CMD+=(--r "$R_VERSION")

export MAMBA_NO_PROGRESS_BARS=1
export MAMBA_NO_BANNER=1
echo "🚀 Setting up environment (this may take a while)..."
"${CMD[@]}"

echo
echo
echo "🧰 Installing fsqc into environment: $ENV_NAME"
# Use the micromamba environment Python directly to avoid PATH conflicts with FreeSurfer/FSL
MAMBA_PYTHON="${MAMBA_ROOT_PREFIX}/envs/${ENV_NAME}/bin/python"

if [[ ! -f "$MAMBA_PYTHON" ]]; then
  echo "❌ Python not found in micromamba environment: $MAMBA_PYTHON"
  echo "Environment creation may have failed."
  exit 1
fi

echo "  → Using Python: $MAMBA_PYTHON"
echo "  → Upgrading pip, wheel, setuptools..."
"$MAMBA_PYTHON" -m pip install --upgrade pip wheel setuptools

# Install fsqc with recommended dependencies
echo "  → Installing fsqc from PyPI..."
"$MAMBA_PYTHON" -m pip install fsqc

# Verify installation using the same Python
echo "  → Verifying fsqc installation..."
if "$MAMBA_PYTHON" -c "import fsqc; print('fsqc module imported successfully')" >/dev/null 2>&1; then
  echo "✅ fsqc Python module installed successfully."
  
  # Check if run_fsqc command is available in the environment
  RUN_FSQC_PATH="${MAMBA_ROOT_PREFIX}/envs/${ENV_NAME}/bin/run_fsqc"
  if [[ -f "$RUN_FSQC_PATH" ]]; then
    echo "✅ run_fsqc command available at: $RUN_FSQC_PATH"
    echo "  → Testing run_fsqc help..."
    if "$RUN_FSQC_PATH" --help >/dev/null 2>&1; then
      echo "  → run_fsqc help command works correctly."
    else
      echo "  ⚠️  run_fsqc command exists but help failed."
    fi
  else
    echo "  ⚠️  run_fsqc command not found at expected location."
    echo "  → This may be normal - fsqc might use a different entry point."
  fi
else
  echo "❌ fsqc installation failed - Python module not importable."
  echo
  echo "Installation failed. Please report this issue."
  exit 1
fi

echo
echo "🎉 Done"
echo "👉 Activate the environment:"
echo "   source scripts/mamba_activate.sh"
echo "👉 Verify tools:"
echo "   which Rscript"
echo "   python -c 'import fsqc; print(\"fsqc available\")'"
echo "   # Note: run_fsqc may require full path due to FreeSurfer/FSL PATH conflicts"
echo "👉 R helper usage:"
echo "   Rscript scripts/fslmer_univariate.R --help"
echo "👉 Deactivate later:"
echo "   source scripts/mamba_deactivate.sh"
