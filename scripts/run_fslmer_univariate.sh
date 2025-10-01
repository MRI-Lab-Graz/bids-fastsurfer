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

# Show error function
show_error() {
    echo "ERROR: $1" >&2
    echo "Use --help for usage examples and argument details." >&2
    exit 1
}

# Change to project root directory
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR" || show_error "Cannot change to project directory: $ROOT_DIR"

# Validate that R script exists
R_SCRIPT="scripts/fslmer_univariate.R"
if [[ ! -f "$R_SCRIPT" ]]; then
    show_error "R script not found: $R_SCRIPT"
fi

# Environment setup - try micromamba first
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_RECORD="$SCRIPT_DIR/.mamba_env"
ENV_NAME="fastsurfer-r"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/micromamba}"
MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"

# Load environment record if it exists
if [[ -f "$ENV_RECORD" ]]; then
    # shellcheck disable=SC1090
    source "$ENV_RECORD"
    ENV_NAME="${ENV_NAME:-fastsurfer-r}"
    MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/.local/micromamba}"
    MAMBA_BIN="$MAMBA_ROOT_PREFIX/bin/micromamba"
fi

# Check for required R packages function
check_r_packages() {
    local rscript_cmd="$1"
    echo "Checking R package dependencies..." >&2
    
    if ! $rscript_cmd -e "
        required_pkgs <- c('optparse', 'jsonlite', 'stats')
        missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]
        if (length(missing) > 0) {
            cat('Missing required R packages:', paste(missing, collapse=', '), '\n', file=stderr())
            cat('Run scripts/install.sh to install required packages.\n', file=stderr())
            quit(status=1)
        }
        cat('All required R packages available.\n', file=stderr())
    " 2>/dev/null; then
        echo "ERROR: Missing required R packages. Run scripts/install.sh to install dependencies." >&2
        return 1
    fi
}

# Try micromamba environment first
if [[ -x "$MAMBA_BIN" && -d "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME" ]]; then
    echo "Using micromamba environment: $ENV_NAME" >&2
    RSCRIPT_CMD="$MAMBA_BIN run -p $MAMBA_ROOT_PREFIX/envs/$ENV_NAME Rscript"
    
    # Verify environment works
    if check_r_packages "$RSCRIPT_CMD"; then
        echo "Executing R script with micromamba..." >&2
        exec $RSCRIPT_CMD "$R_SCRIPT" "$@"
    else
        echo "WARNING: micromamba environment has missing packages, falling back to system R" >&2
    fi
fi

# Fallback: system Rscript
if ! command -v Rscript >/dev/null 2>&1; then
    show_error "Rscript not found in PATH and micromamba environment unavailable. Install R or run scripts/install.sh to set up environment."
fi

echo "Using system R installation" >&2
RSCRIPT_CMD="Rscript"

# Check system R packages
if ! check_r_packages "$RSCRIPT_CMD"; then
    show_error "System R missing required packages. Run scripts/install.sh or install packages manually: install.packages(c('optparse', 'jsonlite'))"
fi

echo "Executing R script with system R..." >&2
exec Rscript "$R_SCRIPT" "$@"
