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
#   {outdir}/summary_tp_sexM.csv:     Significance summary with categories (sig/trend/non-sig)
#   {outdir}/summary.html:            HTML report with plots for significant results
#   {outdir}/plots/:                  Directory containing plots for significant ROIs

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

# Parse arguments to extract output directory for summarization
OUTDIR=""
CONFIG_FILE=""
ARGS=("$@")

for ((i=0; i<${#ARGS[@]}; i++)); do
    case ${ARGS[i]} in
        --outdir)
            if [[ $((i+1)) -lt ${#ARGS[@]} ]]; then
                OUTDIR="${ARGS[$((i+1))]}"
            fi
            ;;
        --config)
            if [[ $((i+1)) -lt ${#ARGS[@]} ]]; then
                CONFIG_FILE="${ARGS[$((i+1))]}"
            fi
            ;;
    esac
done

# If config file was used, extract outdir from it
if [[ -n "$CONFIG_FILE" && -z "$OUTDIR" ]]; then
    if command -v jq >/dev/null 2>&1; then
        OUTDIR=$(jq -r '.outdir // empty' "$CONFIG_FILE" 2>/dev/null || echo "")
    elif command -v python3 >/dev/null 2>&1; then
        OUTDIR=$(python3 -c "import json; print(json.load(open('$CONFIG_FILE')).get('outdir', ''))" 2>/dev/null || echo "")
    fi
fi

# Setup R environment and get Rscript command
RSCRIPT_CMD="$(setup_r_environment "$SCRIPT_DIR")"

echo "Executing R script..." >&2
$RSCRIPT_CMD "$R_SCRIPT" "$@"

# Check if analysis completed successfully and run summarization
if [[ $? -eq 0 ]]; then
    # Run summarization if we have an output directory
    if [[ -n "$OUTDIR" && -d "$OUTDIR" ]]; then
        echo "Generating summary report..." >&2
        SUMMARIZE_SCRIPT="scripts/fslmer_summarize.R"
        if [[ -f "$SUMMARIZE_SCRIPT" ]]; then
            # Use sexM as default effect for single ROI analyses
            $RSCRIPT_CMD "$SUMMARIZE_SCRIPT" --results-dir "$OUTDIR" --effect "sexM" --verbose
        else
            echo "Warning: Summary script not found at $SUMMARIZE_SCRIPT" >&2
        fi
    fi
fi
