#!/usr/bin/env bash
set -euo pipefail

# run_fslmer_univariate.sh
# Purpose: Run scripts/fslmer_univariate.R within the project's renv so that
# required packages (optparse, jsonlite, fslmer, mgcv, etc.) are available.
#
# Usage:
#   bash scripts/run_fslmer_univariate.sh --config /path/to/config.json [other args]
#   # arguments are passed through to the R script

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found in PATH" >&2; exit 127; }

if [[ ! -f "renv/activate.R" ]]; then
  echo "renv not initialized in this project (missing renv/activate.R). Run: bash scripts/setup_r_env.sh" >&2
  exit 1
fi

# Run the analysis script inside renv, forwarding all CLI args
Rscript -e 'args <- commandArgs(trailingOnly=TRUE); renv::run("scripts/fslmer_univariate.R", args=args)' -- "$@"
