#!/usr/bin/env bash
set -euo pipefail

# Deprecated: install_fsqc.sh
# fsqc is now installed directly into the micromamba environment by scripts/install.sh.
# This stub remains only for backward compatibility and to provide guidance.

echo "[DEPRECATED] scripts/install_fsqc.sh is no longer used."
echo "Please run: bash scripts/install.sh (optionally with --no-fsqc to skip)."
echo "After installation, activate the env and verify:"
echo "  source scripts/mamba_activate.sh && run_fsqc --help"
exit 1
