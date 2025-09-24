#!/usr/bin/env bash

# mamba_activate.sh
# Purpose: Activate the fastsurfer-r micromamba environment like Python venv.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-${(%):-%N}}")" && pwd)"
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
if [[ -f "$ENV_RECORD" ]]; then
  # shellcheck disable=SC1090
  source "$ENV_RECORD"
fi

# Fallbacks if .mamba_env not present or invalid
: "${MAMBA_ROOT_PREFIX:=${HOME}/.local/micromamba}"
: "${ENV_NAME:=fastsurfer-r}"

# Export to ensure downstream hook/commands see the intended root
export MAMBA_ROOT_PREFIX

MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
# If recorded root is stale, fall back to default HOME path
if [[ ! -x "$MAMBA_BIN" && "$MAMBA_ROOT_PREFIX" != "$HOME/.local/micromamba" ]]; then
  MAMBA_ROOT_PREFIX="$HOME/.local/micromamba"
  MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
fi
CURRENT_SHELL="bash"
if [[ -n "${ZSH_VERSION:-}" ]] || [[ "${SHELL:-}" == *"zsh"* ]]; then
  CURRENT_SHELL="zsh"
fi

if [[ -x "$MAMBA_BIN" ]]; then
  # Always prefer the recorded micromamba for consistent root
  eval "$("$MAMBA_BIN" shell hook -s "$CURRENT_SHELL")"
else
  echo "micromamba not found at $MAMBA_ROOT_PREFIX/bin/micromamba. Run: bash scripts/mamba_setup.sh" >&2
  return 1 2>/dev/null || exit 1
fi

# Activate by env name using the recorded root (hook points to the right micromamba)
ENV_PREFIX="$MAMBA_ROOT_PREFIX/envs/$ENV_NAME"
# Use the shell-integrated activation so environment variables are applied to current shell
if ! micromamba activate "$ENV_NAME" >/dev/null 2>&1; then
  echo "Failed to activate micromamba env: $ENV_NAME" >&2
  echo "Try running: eval \"$($MAMBA_BIN shell hook -s $CURRENT_SHELL)\" && micromamba activate $ENV_NAME" >&2
  return 1 2>/dev/null || exit 1
fi
RS_BIN="$(command -v Rscript || true)"
ENV_BIN_DIR="$ENV_PREFIX/bin"
if [[ -z "$RS_BIN" || "$RS_BIN" != "$ENV_BIN_DIR/Rscript" ]]; then
  echo "Warning: Rscript not resolving to env bin. Got: ${RS_BIN:-not found}"
  echo "If PATH didn't update, try manually: eval \"$($MAMBA_BIN shell hook -s $CURRENT_SHELL)\" && micromamba activate $ENV_NAME"
fi
echo "Activated micromamba env: $ENV_NAME ($ENV_PREFIX)"

# Source local FreeSurfer setup so tools like asegstats2table are available
PREFERRED_FS_HOME="${FREESURFER_HOME:-}"
if [[ -n "$PREFERRED_FS_HOME" && -f "$PREFERRED_FS_HOME/SetUpFreeSurfer.sh" ]]; then
  export FREESURFER_HOME="$PREFERRED_FS_HOME"
elif [[ -f "/usr/local/freesurfer/SetUpFreeSurfer.sh" ]]; then
  export FREESURFER_HOME="/usr/local/freesurfer"
else
  echo "Warning: SetUpFreeSurfer.sh not found; asegstats2table may be unavailable." >&2
  FREESURFER_HOME=""
fi

if [[ -n "$FREESURFER_HOME" ]]; then
  # shellcheck disable=SC1091
  source "$FREESURFER_HOME/SetUpFreeSurfer.sh"
  if [[ -n "${FASTSURFER_SUBJECTS_DIR:-}" ]]; then
    export SUBJECTS_DIR="$FASTSURFER_SUBJECTS_DIR"
  fi
fi

# Avoid renv autoloader and user/site profiles interfering in this env
export RENV_CONFIG_AUTOLOADER_ENABLED=false
# Disable reading of user and site profiles; users can still run R without wrapper if needed
export R_PROFILE_USER=
export R_PROFILE=

# Provide a convenience wrapper so plain 'Rscript' runs with --vanilla for scripts
# but does not interfere with informational flags like --help/--version.
# Use a simple alias when possible to avoid function export issues across shells.
if [[ -n "${BASH_VERSION:-}" ]]; then
  Rscript() {
    if [[ $# -eq 0 ]]; then
      command Rscript "$@"
      return $?
    fi
    case "$1" in
      --help|--version|-h|-v)
        command Rscript "$@" ;;
      *)
        command Rscript --vanilla "$@" ;;
    esac
  }
  export -f Rscript 2>/dev/null || true
elif [[ -n "${ZSH_VERSION:-}" ]]; then
  # zsh: define function (no export needed) to avoid alias side-effects in tools like `which`
  Rscript() {
    if [[ $# -eq 0 ]]; then
      command Rscript "$@"
      return $?
    fi
    case "$1" in
      --help|--version|-h|-v)
        command Rscript "$@" ;;
      *)
        command Rscript --vanilla "$@" ;;
    esac
  }
fi
