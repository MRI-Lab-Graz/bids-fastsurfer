#!/usr/bin/env bash
set -euo pipefail

# run_fslmer_univariate.sh
# Purpose: Run scripts/fslmer_univariate.R inside the project's micromamba env if present,
# otherwise fall back to system Rscript.
#
# Usage:
#   bash scripts/run_fslmer_univariate.sh --config /path/to/config.json [other args]
#   # arguments are passed through to the R script

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

# Try micromamba first
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
if [[ -f "$ENV_RECORD" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_RECORD"
fi
ENV_NAME="${ENV_NAME:-fastsurfer-r}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/micromamba}"
MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"

if [[ -x "$MAMBA_BIN" && -d "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME" ]]; then
  exec "$MAMBA_BIN" run -p "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME" Rscript scripts/fslmer_univariate.R "$@"
fi

# Fallback: system Rscript
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found in PATH and micromamba env unavailable" >&2; exit 127; }
exec Rscript scripts/fslmer_univariate.R "$@"
