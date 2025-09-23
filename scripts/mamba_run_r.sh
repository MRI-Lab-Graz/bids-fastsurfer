#!/usr/bin/env bash
set -euo pipefail

# mamba_run_r.sh
# Purpose: Run scripts/fslmer_univariate.R inside the micromamba env without manual activation.
#
# Usage:
#   bash scripts/mamba_run_r.sh --config /path/to/config.json [other args]

ENV_NAME="fastsurfer-r"

# Resolve micromamba binary in this order:
# 1) MAMBA_ROOT_PREFIX/bin/micromamba if env var is set
# 2) micromamba on PATH
# 3) ~/.local/micromamba/bin/micromamba (default)

if [[ -n "${MAMBA_ROOT_PREFIX:-}" && -x "$MAMBA_ROOT_PREFIX/bin/micromamba" ]]; then
  MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
elif command -v micromamba >/dev/null 2>&1; then
  MAMBA_BIN="$(command -v micromamba)"
elif [[ -x "$HOME/.local/micromamba/bin/micromamba" ]]; then
  MAMBA_BIN="$HOME/.local/micromamba/bin/micromamba"
else
  echo "micromamba not found. Run setup: bash scripts/mamba_setup.sh" >&2
  exit 1
fi

exec "$MAMBA_BIN" run -n "$ENV_NAME" Rscript scripts/fslmer_univariate.R "$@"
