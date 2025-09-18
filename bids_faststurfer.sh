#!/bin/bash



# Usage: ./bids_fastsurfer.sh /bids-folder /outputfolder -c fastsurfer_options.json [--dry_run] [--pilot] [--debug]
#
# Arguments:
#   /bids-folder           Path to BIDS input directory
#   /outputfolder          Path to output directory
#   -c <config.json>       Path to JSON config file with FastSurfer options
#   --dry_run              (Optional) Print the Singularity command instead of running it
#   --pilot                (Optional) Randomly pick only one subject for testing
#   --sub <subject>        (Optional) Process only the specified subject (e.g. sub-001)
#   --ses <session>        (Optional) Process only the specified session for the subject (e.g. ses-1)
#   --debug                (Optional) Print debug information about parsed options and paths
#
# Example:
#   ./bids_fastsurfer.sh ./data ./output -c fastsurfer_options.json --pilot --dry_run --debug
    echo "\nExample:"
    echo "  $0 ./data ./output -c fastsurfer_options.json --pilot --dry_run --debug"
    echo "  $0 ./data ./output -c fastsurfer_options.json --sub sub-001 --ses ses-1"

# Show help if no arguments are provided
SUBJECT=""
SESSION=""
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
        --pilot)
            PILOT=1
            shift
            ;;
        --debug)
            DEBUG=1
            shift
            ;;
        --sub)
            SUBJECT="$2"
            shift 2
            ;;
        --ses)
            SESSION="$2"
            shift 2
            ;;
        *)
            if [[ -z "$BIDS_DATA" ]]; then
                BIDS_DATA="$1"
            elif [[ -z "$OUTPUT_DIR" ]]; then
                OUTPUT_DIR="$1"
            else
                echo "Unknown argument: $1"
                exit 1
            fi
            shift
            ;;
    esac
done
DEBUG=0


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

        --pilot)
            PILOT=1
            shift
            ;;
        --debug)
            DEBUG=1
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
LICENSE_DIR=$(dirname "$LICENSE_PATH")

# Parse options from JSON config (except sid, t1, py, fs_license, and sif_file)
parse_json_options() {
    jq -r 'to_entries[] | select(.key != "fs_license" and .key != "sif_file" and .value != null and .value != false and .value != "") | "--" + .key + (if (.value|type) == "boolean" then "" else " " + (.value|tostring) end)' "$1"
}


# Build T1w list based on --sub/--ses or all
if [[ -n "$SUBJECT" ]]; then
    # Check subject exists
    if [[ ! -d "$BIDS_DATA/$SUBJECT" ]]; then
        echo "Error: Subject $SUBJECT not found in $BIDS_DATA."
        exit 1
    fi
    if [[ -n "$SESSION" ]]; then
        # Check session exists
        if [[ ! -d "$BIDS_DATA/$SUBJECT/$SESSION" ]]; then
            echo "Error: Session $SESSION not found for $SUBJECT in $BIDS_DATA."
            exit 1
        fi
        T1W_LIST=( $(find "$BIDS_DATA/$SUBJECT/$SESSION" -type f \( -name "*_T1w.nii" -o -name "*_T1w.nii.gz" -o -name "*_desc-preproc_T1w.nii.gz" \)) )
    else
        T1W_LIST=( $(find "$BIDS_DATA/$SUBJECT" -type f \( -name "*_T1w.nii" -o -name "*_T1w.nii.gz" -o -name "*_desc-preproc_T1w.nii.gz" \)) )
    fi
elif [[ -n "$SESSION" ]]; then
    # No subject, but session specified: process all subjects with this session
    mapfile -t SUBJECTS < <(find "$BIDS_DATA" -maxdepth 1 -type d -name 'sub-*' -exec basename {} \;)
    T1W_LIST=()
    for subj in "${SUBJECTS[@]}"; do
        if [[ -d "$BIDS_DATA/$subj/$SESSION" ]]; then
            t1s=( $(find "$BIDS_DATA/$subj/$SESSION" -type f \( -name "*_T1w.nii" -o -name "*_T1w.nii.gz" -o -name "*_desc-preproc_T1w.nii.gz" \)) )
            T1W_LIST+=("${t1s[@]}")
        fi
    done
    if [[ ${#T1W_LIST[@]} -eq 0 ]]; then
        echo "Error: No T1w images found for session $SESSION in $BIDS_DATA."
        exit 1
    fi
else
    T1W_LIST=( $(find "$BIDS_DATA" -type f \( -name "*_T1w.nii" -o -name "*_T1w.nii.gz" -o -name "*_desc-preproc_T1w.nii.gz" \)) )
fi

# If --pilot, randomly pick one T1w image
if [[ $PILOT -eq 1 ]]; then
    if [[ ${#T1W_LIST[@]} -eq 0 ]]; then
        echo "No T1w images found."
        exit 1
    fi
    RANDOM_IDX=$(( RANDOM % ${#T1W_LIST[@]} ))
    T1W_LIST=( "${T1W_LIST[$RANDOM_IDX]}" )
    echo "[PILOT MODE] Only processing: ${T1W_LIST[0]}"
fi

for t1w_img in "${T1W_LIST[@]}"; do
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
    if [[ $DEBUG -eq 1 ]]; then
        echo "DEBUG: Extra options from JSON: $extra_opts"
    fi

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
    # Ensure BIDS_DATA ends with a single slash
    BIDS_DATA_SLASH="$BIDS_DATA"
    [[ "${BIDS_DATA_SLASH}" != */ ]] && BIDS_DATA_SLASH="${BIDS_DATA_SLASH}/"
    # Compute relative path from BIDS_DATA to t1w_img
    t1w_relpath="${t1w_img#${BIDS_DATA_SLASH}}"
    if [[ $DEBUG -eq 1 ]]; then
        echo "DEBUG: Using --t1 /data/$t1w_relpath"
    fi
    cmd=(singularity exec --nv \
        --no-home \
        -B "$BIDS_DATA":/data \
        -B "$OUTPUT_DIR":/output \
        -B "$LICENSE_DIR":/fs_license \
        "$SIF_FILE" \
        /fastsurfer/run_fastsurfer.sh \
        --t1 "/data/$t1w_relpath" \
        --sid "$sid" \
        --sd /output \
        --fs_license /fs_license/license.txt \
        --3T)
    # Add extra options, splitting on spaces
    if [[ -n "$extra_opts" ]]; then
        # shellcheck disable=SC2206
        extra_opts_arr=($extra_opts)
        cmd+=("${extra_opts_arr[@]}")
    fi

    echo "Processing $t1w_img with --sid $sid"
    if [[ $DRY_RUN -eq 1 ]]; then
        printf '%q ' "${cmd[@]}"
        echo
    else
        "${cmd[@]}"
    fi
done