#!/bin/bash


# Usage: ./bids_fastsurfer.sh /bids-folder /outputfolder -c fastsurfer_options.json [--dry_run]
#
# Arguments:
#   /bids-folder           Path to BIDS input directory
#   /outputfolder          Path to output directory
#   -c <config.json>       Path to JSON config file with FastSurfer options
#   --dry_run              (Optional) Print the Singularity command instead of running it
#
# Example:
#   ./bids_fastsurfer.sh ./data ./output -c fastsurfer_options.json --dry_run

# Show help if no arguments are provided
if [[ $# -eq 0 ]]; then
    echo "\nUsage: $0 <bids_data_dir> <output_dir> -c <config.json> [--dry_run]"
    echo "\nArguments:"
    echo "  <bids_data_dir>   Path to BIDS input directory"
    echo "  <output_dir>      Path to output directory"
    echo "  -c <config.json>  Path to JSON config file with FastSurfer options"
    echo "  --dry_run         (Optional) Print the Singularity command instead of running it"
    echo "\nExample:"
    echo "  $0 ./data ./output -c fastsurfer_options.json --dry_run"
    exit 0
fi

set -e

if ! command -v jq &> /dev/null; then
    echo "Error: jq is required but not installed. Please install jq."
    exit 1
fi

BIDS_DATA="$1"
OUTPUT_DIR="$2"
shift 2


CONFIG=""
DRY_RUN=0


while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        --dry_run)
            DRY_RUN=1
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Resolve CONFIG to absolute path
if [[ -n "$CONFIG" && ! "$CONFIG" =~ ^/ ]]; then
    CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
fi


# Check required arguments
if [[ -z "$BIDS_DATA" || -z "$OUTPUT_DIR" || -z "$CONFIG" ]]; then
    echo "Usage: $0 <bids_data_dir> <output_dir> -c <config.json> [--dry_run]"
    exit 1
fi

# Check if input/output folders and license file exist
if [[ ! -d "$BIDS_DATA" ]]; then
    echo "Error: BIDS input directory '$BIDS_DATA' does not exist."
    exit 1
fi
if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory '$OUTPUT_DIR' does not exist."
    exit 1
fi

LICENSE_PATH=$(jq -r .fs_license "$CONFIG")
if [[ ! -f "$LICENSE_PATH" ]]; then
    echo "Error: FreeSurfer license file '$LICENSE_PATH' does not exist."
    exit 1
fi

# Parse options from JSON config (except sid, t1, py)
parse_json_options() {
    jq -r 'to_entries[] | select(.value != null and .value != false and .value != "") | "--" + .key + (if (.value|type) == "boolean" then "" else " \"" + (.value|tostring) + "\"" end)' "$1"
}

# Find all T1w images (including preproc)
find "$BIDS_DATA" -type f \( -name "*_T1w.nii" -o -name "*_T1w.nii.gz" -o -name "*_desc-preproc_T1w.nii.gz" \) | while read -r t1w_img; do
    fname=$(basename "$t1w_img")
    subj=$(echo "$fname" | grep -o 'sub-[^_]*')
    sess=$(echo "$fname" | grep -o 'ses-[^_]*')
    if [[ -n "$sess" ]]; then
        sid="${subj}_${sess}"
    else
        sid="$subj"
    fi

    # Try to find T2w image for this subject/session
    t2w_img=""
    t2_pattern="${subj}"
    if [[ -n "$sess" ]]; then
        t2_pattern="${subj}_${sess}"
    fi
    t2w_img=$(find "$BIDS_DATA" -type f -name "${t2_pattern}_*T2w.nii.gz" | head -n1)

    # Build options from JSON config
    extra_opts=$(parse_json_options "$CONFIG")

    # If t2w_img exists, add --t2
    if [[ -n "$t2w_img" ]]; then
        extra_opts="$extra_opts --t2 /data/$(basename "$t2w_img")"
    fi

    # Always set --sid, --t1, --sd, --fs_license, --3T
    SIF_FILE=$(jq -r .sif_file "$CONFIG")
    if [[ ! -f "$SIF_FILE" ]]; then
        echo "Error: Singularity image file '$SIF_FILE' does not exist."
        exit 1
    fi
    cmd="singularity exec --nv \
        --no-home \
        -B \"$BIDS_DATA\":/data \
        -B \"$OUTPUT_DIR\":/output \
        -B \"$LICENSE_PATH\":/fs_license \
        $SIF_FILE \
        /fastsurfer/run_fastsurfer.sh \
        --t1 /data/$(basename \"$t1w_img\") \
        --sid \"$sid\" \
        --sd /output \
        --fs_license /fs_license/license.txt \
        --3T $extra_opts"

    echo "Processing $t1w_img with --sid $sid"
    if [[ $DRY_RUN -eq 1 ]]; then
        echo "$cmd"
    else
        eval $cmd
    fi
done