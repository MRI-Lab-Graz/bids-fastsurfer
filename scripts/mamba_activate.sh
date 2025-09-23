#!/usr/bin/env bash
set -euo pipefail

# mamba_activate.sh
# Purpose: Activate the fastsurfer-r micromamba environment like Python venv.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-${(%):-%N}}")" && pwd)"
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
if [[ -f "$ENV_RECORD" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_RECORD"
fi

# Fallbacks if .mamba_env not present
: "${MAMBA_ROOT_PREFIX:=${HOME}/.local/micromamba}"
: "${ENV_NAME:=fastsurfer-r}"

MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
CURRENT_SHELL="bash"
if [[ -n "${ZSH_VERSION:-}" ]] || [[ "${SHELL:-}" == *"zsh"* ]]; then
  CURRENT_SHELL="zsh"
fi

if ! command -v micromamba >/dev/null 2>&1; then
  if [[ -x "$MAMBA_BIN" ]]; then
    eval "$("$MAMBA_BIN" shell hook -s "$CURRENT_SHELL")"
  else
    echo "micromamba not found. Run: bash scripts/mamba_setup.sh" >&2
    return 1 2>/dev/null || exit 1
  fi
fi

# If micromamba is available on PATH, still eval the shell hook to expose activate
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook -s "$CURRENT_SHELL")"
fi

micromamba activate "$ENV_NAME" >/dev/null
RS_BIN="$(command -v Rscript || true)"
ENV_BIN_DIR="$MAMBA_ROOT_PREFIX/envs/$ENV_NAME/bin"
if [[ -z "$RS_BIN" || "$RS_BIN" != "$ENV_BIN_DIR/Rscript" ]]; then
  echo "Warning: Rscript not resolving to env bin. Got: ${RS_BIN:-not found}"
  echo "If PATH didn't update, try manually: eval \"$($MAMBA_BIN shell hook -s $CURRENT_SHELL)\" && micromamba activate $ENV_NAME"
fi
echo "Activated micromamba env: $ENV_NAME"

# Avoid renv autoloader and user/site profiles interfering in this env
export RENV_CONFIG_AUTOLOADER_ENABLED=false
# Disable reading of user and site profiles; users can still run R without wrapper if needed
export R_PROFILE_USER=
export R_PROFILE=

# Provide a convenience wrapper so plain 'Rscript' runs with --vanilla in this shell
# This limits surprises from project or user profiles trying to source renv/activate.R
Rscript() {
  command Rscript --vanilla "$@"
}
# 'export -f' is bash-only; safe to skip in zsh
if [[ -n "${BASH_VERSION:-}" ]]; then
  export -f Rscript 2>/dev/null || true
fi