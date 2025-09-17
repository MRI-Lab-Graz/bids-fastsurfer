#!/bin/bash

# Usage: ./bids_fastsurfer.sh /bids-folder /outputfolder -c fastsurfer_options.json [--dry_run]

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

if [[ -z "$BIDS_DATA" || -z "$OUTPUT_DIR" || -z "$CONFIG" ]]; then
    echo "Usage: $0 <bids_data_dir> <output_dir> -c <config.json> [--dry_run]"
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
    cmd="singularity exec --nv \
        --no-home \
        -B \"$BIDS_DATA\":/data \
        -B \"$OUTPUT_DIR\":/output \
        -B \"$(jq -r .fs_license "$CONFIG")\":/fs_license \
        ./fastsurfer-gpu.sif \
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