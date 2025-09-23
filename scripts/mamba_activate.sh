#!/usr/bin/env bash
set -euo pipefail

# mamba_activate.sh
# Purpose: Activate the fastsurfer-r micromamba environment like Python venv.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
if [[ -f "$ENV_RECORD" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_RECORD"
fi

# Fallbacks if .mamba_env not present
: "${MAMBA_ROOT_PREFIX:=${HOME}/.local/micromamba}"
: "${ENV_NAME:=fastsurfer-r}"

MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
if ! command -v micromamba >/dev/null 2>&1; then
  if [[ -x "$MAMBA_BIN" ]]; then
    eval "$("$MAMBA_BIN" shell hook -s bash)"
  else
    echo "micromamba not found. Run: bash scripts/mamba_setup.sh" >&2
    return 1 2>/dev/null || exit 1
  fi
fi

# If micromamba is available on PATH, still eval the shell hook to expose activate
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook -s bash)"
fi

micromamba activate "$ENV_NAME"
echo "Activated micromamba env: $ENV_NAME"