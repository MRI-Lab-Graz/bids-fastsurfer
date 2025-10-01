#!/usr/bin/env bash
set -euo pipefail

# Deprecated: install_fsqc.sh
# fsqc is now automatically installed as part of the standard environment setup.
# This stub remains only for backward compatibility and to provide guidance.

echo "[DEPRECATED] scripts/install_fsqc.sh is no longer used."
echo "fsqc is now installed automatically with the main environment."
echo
echo "To set up the complete environment including fsqc:"
echo "  bash scripts/install.sh"
echo
echo "After installation, activate the environment and verify:"
echo "  source scripts/mamba_activate.sh"
echo "  run_fsqc --help"
echo
echo "If fsqc is missing from an existing environment, reinstall:"
echo "  bash scripts/install.sh  # This will recreate the environment with fsqc"
exit 1
