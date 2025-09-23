#!/usr/bin/env bash
# Note: source this file to affect your current shell
#   zsh/bash:  source scripts/mamba_deactivate.sh

# Avoid exiting the parent shell on errors
set +e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-${(%):-%N}}")" && pwd)"
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
if [[ -f "$ENV_RECORD" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_RECORD"
fi

: "${MAMBA_ROOT_PREFIX:=${HOME}/.local/micromamba}"
MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
# If recorded root is stale, fall back to default HOME path
if [[ ! -x "$MAMBA_BIN" && "$MAMBA_ROOT_PREFIX" != "$HOME/.local/micromamba" ]]; then
  MAMBA_ROOT_PREFIX="$HOME/.local/micromamba"
  MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
fi

# Detect shell for hook
CURRENT_SHELL="bash"
if [[ -n "${ZSH_VERSION:-}" ]] || [[ "${SHELL:-}" == *"zsh"* ]]; then
  CURRENT_SHELL="zsh"
fi

if [[ ! -x "$MAMBA_BIN" ]]; then
  echo "[mamba_deactivate] micromamba not found at $MAMBA_ROOT_PREFIX/bin/micromamba"
  echo "[mamba_deactivate] If the env is active, you can try: micromamba deactivate"
else
  # Ensure shell integration is available in this session, then deactivate
  eval "$("$MAMBA_BIN" shell hook -s "$CURRENT_SHELL")"
  micromamba deactivate >/dev/null 2>&1 || true
fi

# Clean up Rscript wrapper defined by mamba_activate.sh
if [[ -n "${BASH_VERSION:-}" ]]; then
  unset -f Rscript 2>/dev/null || true
elif [[ -n "${ZSH_VERSION:-}" ]]; then
  unfunction Rscript 2>/dev/null || true
fi

# Remove env vars set by activator to avoid interference
unset RENV_CONFIG_AUTOLOADER_ENABLED
unset R_PROFILE_USER
unset R_PROFILE

echo "Deactivated micromamba (if active) and cleaned up Rscript wrapper."
echo "If your prompt still shows (fastsurfer-r), run: micromamba deactivate"
