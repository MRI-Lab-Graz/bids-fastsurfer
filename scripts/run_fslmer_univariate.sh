#!/usr/bin/env bash
set -euo pipefail

# run_fslmer_univariate.sh
# Purpose: Run univariate linear mixed-effects analysis for FreeSurfer/FastSurfer longitudinal data.
# This wrapper script manages the R environment and provides user-friendly interfaces.
#
# DESCRIPTION:
#   Flexible univariate linear mixed-effects (LME) analysis for subcortical/cortical 
#   volume data from FreeSurfer/FastSurfer longitudinal processing. Supports single ROI 
#   or multi-region analysis with various statistical models (fslmer, GLM, GAM).
#
# USAGE:
#   # Show help
#   bash scripts/run_fslmer_univariate.sh --help
#
#   # Single ROI analysis
#   bash scripts/run_fslmer_univariate.sh \\
#     --qdec results/prep_long/qdec.table.dat \\
#     --aseg results/prep_long/aseg.long.table \\
#     --roi Left.Hippocampus \\
#     --formula '~ tp*group' \\
#     --outdir results/hippo_analysis
#
#   # Multiple regions matching pattern
#   bash scripts/run_fslmer_univariate.sh \\
#     --qdec results/prep_long/qdec.table.dat \\
#     --aseg results/prep_long/aseg.long.table \\
#     --region-pattern 'Hippocampus|Amygdala' \\
#     --formula '~ tp*group' \\
#     --outdir results/limbic_analysis
#
#   # All brain regions analysis
#   bash scripts/run_fslmer_univariate.sh \\
#     --qdec results/prep_long/qdec.table.dat \\
#     --aseg results/prep_long/aseg.long.table \\
#     --all-regions \\
#     --formula '~ tp' \\
#     --outdir results/full_analysis
#
#   # Using configuration file
#   bash scripts/run_fslmer_univariate.sh --config configs/fslmer_univariate.example.json
#
# EXAMPLES:
#   # Basic longitudinal change
#   --formula '~ tp'
#   
#   # Group by time interaction
#   --formula '~ tp*group'
#   
#   # With baseline covariate
#   --add-baseline --formula '~ baseline_value + tp*group'
#   
#   # GAM with smooth time effect
#   --engine gam --formula '~ s(tp)'
#
# INPUT FILES:
#   qdec.table.dat:     Subject metadata (TSV) with fsid, fsid_base, tp, group, etc.
#   aseg.long.table:    FreeSurfer longitudinal volume table with .long. identifiers
#
# OUTPUT:
#   {outdir}/results_summary.csv:     Summary statistics for all ROIs
#   {outdir}/models/{roi}_model.txt:  Individual model results per ROI
#   {outdir}/merged_data.csv:         Combined input data (if --save-merged)

# Load shared functions
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/common_functions.sh"

# Handle help flag early
if [[ $# -eq 0 ]] || [[ "${1:-}" == "--help" ]] || [[ "${1:-}" == "-h" ]]; then
    # Extract help text from this script
    sed -n '/^# DESCRIPTION:/,/^[^#]/p' "$0" | sed 's/^# \?//' | head -n -1
    echo ""
    echo "R SCRIPT ARGUMENT DETAILS:"
    echo "Run with R script directly for detailed argument documentation:"
    echo "  Rscript scripts/fslmer_univariate.R --help"
    echo ""
    echo "ENVIRONMENT:"
    echo "  This script automatically manages R environment using micromamba if available,"
    echo "  otherwise falls back to system Rscript. Run scripts/install.sh to set up"
    echo "  the recommended environment with all required R packages."
    exit 0
fi

# Change to project root directory
ROOT_DIR="$(get_project_root "$SCRIPT_DIR")"
cd "$ROOT_DIR" || show_error "Cannot change to project directory: $ROOT_DIR"

# Validate that R script exists
R_SCRIPT="scripts/fslmer_univariate.R"
validate_path "$R_SCRIPT" file "R script" || exit 1

# Setup R environment and get Rscript command
RSCRIPT_CMD="$(setup_r_environment "$SCRIPT_DIR")"

echo "Executing R script..." >&2
exec $RSCRIPT_CMD "$R_SCRIPT" "$@"
