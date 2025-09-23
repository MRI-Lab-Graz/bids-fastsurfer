#!/usr/bin/env bash
# Deprecated shim. Please use the new path:
#   bash scripts/install.sh [options]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
if [[ -x "scripts/install.sh" ]]; then
  echo "[install.sh] Note: this entrypoint is deprecated. Forwarding to scripts/install.sh ..." >&2
  exec bash scripts/install.sh "$@"
else
  echo "[install.sh] Error: scripts/install.sh not found. Please update your commands to use bash scripts/install.sh" >&2
  exit 1
fi
