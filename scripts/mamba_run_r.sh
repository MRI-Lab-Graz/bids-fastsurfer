#!/usr/bin/env bash
set -euo pipefail

# mamba_run_r.sh
# Purpose: Run scripts/fslmer_univariate.R inside the micromamba env without manual activation.
#
# Usage:
#   bash scripts/mamba_run_r.sh --config /path/to/config.json [other args]

ENV_NAME="fastsurfer-r"
PREFIX="$HOME/.local/micromamba"
MAMBA_BIN="$PREFIX/bin/micromamba"

if [[ ! -x "$MAMBA_BIN" ]]; then
  echo "micromamba not found at $MAMBA_BIN. Run: bash scripts/mamba_setup.sh" >&2
  exit 1
fi

exec "$MAMBA_BIN" run -n "$ENV_NAME" Rscript scripts/fslmer_univariate.R "$@"
