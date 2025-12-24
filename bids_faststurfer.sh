#!/bin/bash

# Source shared validation functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/scripts/validate_fastsurfer.sh"

############################################
# Main Script
############################################



# Usage: ./bids_fastsurfer.sh /bids-folder /outputfolder -c fastsurfer_options.json [--dry_run] [--pilot] [--debug] [--nohup] [--batch_size N]
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
#   --nohup                (Optional) Run each subject in the background. Without --batch_size, defaults to sequential (batch_size=1)
#   --batch_size N         (Optional) Process subjects in batches of N (recommended for large re-runs)
#
# Example:
#   ./bids_fastsurfer.sh ./data ./output -c fastsurfer_options.json --pilot --dry_run --debug
#   ./bids_fastsurfer.sh ./data ./output -c fastsurfer_options.json --sub sub-001 --ses ses-1
#   ./bids_fastsurfer.sh ./data ./output -c fastsurfer_options.json --nohup
#   ./bids_fastsurfer.sh ./data ./output -c fastsurfer_options.json --nohup --batch_size 4

# Show help if no arguments are provided
SUBJECT=""
SESSION=""
RERUN_FILE=""
NOHUP=0
BATCH_SIZE=""
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
        --nohup)
            NOHUP=1
            shift
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --re-run)
            RERUN_FILE="$2"
            shift 2
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

# Initialize variables that might not be set
DRY_RUN=${DRY_RUN:-0}
PILOT=${PILOT:-0}
DEBUG=${DEBUG:-0}

# Resolve CONFIG to absolute path
if [[ -n "$CONFIG" && ! "$CONFIG" =~ ^/ ]]; then
    CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
fi


# Check required arguments
if [[ -z "$BIDS_DATA" || -z "$OUTPUT_DIR" || -z "$CONFIG" ]]; then
    echo "Usage: $0 <bids_data_dir> <output_dir> -c <config.json> [--dry_run]"
    exit 1
fi

# Resolve CONFIG to absolute path
if [[ -n "$CONFIG" && ! "$CONFIG" =~ ^/ ]]; then
    CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
fi

# Quick check for jq before validation
if ! command -v jq >/dev/null 2>&1; then
    echo "Error: 'jq' is required but not found in PATH."
    exit 1
fi

# Extract config values for validation
SIF_FILE=$(jq -r .sif_file "$CONFIG" 2>/dev/null)
LICENSE_PATH=$(jq -r .fs_license "$CONFIG" 2>/dev/null)

if [[ -z "$SIF_FILE" || -z "$LICENSE_PATH" ]]; then
    echo "Error: Could not parse sif_file or fs_license from config: $CONFIG"
    exit 1
fi

# Run comprehensive validation
if ! validate_requirements "$BIDS_DATA" "$OUTPUT_DIR" "$CONFIG" "$SIF_FILE" "$LICENSE_PATH"; then
    echo ""
    echo "Validation failed. Please fix the errors above."
    exit 1
fi
echo ""

LICENSE_DIR=$(dirname "$LICENSE_PATH")


# Parse options from JSON config (except sid, t1, py, fs_license, and sif_file) from 'cross' section
parse_json_options_cross() {
    jq -r '.cross | to_entries[] | select(.key != "fs_license" and .key != "sif_file" and .value != null and .value != false and .value != "") | "--" + .key + (if (.value|type) == "boolean" then "" else " " + (.value|tostring) end)' "$1"
}


# Support --re-run JSON (list of subjects) to build T1W list from specific subjects
if [[ -n "$RERUN_FILE" ]]; then
    if [[ ! -f "$RERUN_FILE" ]]; then
        echo "Error: Re-run file '$RERUN_FILE' not found." >&2
        exit 1
    fi

    declare -a RERUN_SUBJECTS=()
    while IFS= read -r s; do
        # accept either full 'sub-XXX' or bare IDs
        RERUN_SUBJECTS+=("$s")
    done < <(jq -r '.subjects[]' "$RERUN_FILE")

    if [[ ${#RERUN_SUBJECTS[@]} -eq 0 ]]; then
        echo "Error: No subjects found in $RERUN_FILE" >&2
        exit 1
    fi

    T1W_LIST=()
    for subj in "${RERUN_SUBJECTS[@]}"; do
        subj_dir="$BIDS_DATA/$subj"
        if [[ ! -d "$subj_dir" ]]; then
            echo "[WARN] Subject directory not found for $subj at $subj_dir - skipping"
            continue
        fi
        t1s=( $(find "$subj_dir" -type f \( -name "*_T1w.nii" -o -name "*_T1w.nii.gz" -o -name "*_desc-preproc_T1w.nii.gz" \)) )
        if [[ ${#t1s[@]} -eq 0 ]]; then
            echo "[WARN] No T1w images found for $subj - skipping"
            continue
        fi
        T1W_LIST+=("${t1s[@]}")
    done

    if [[ ${#T1W_LIST[@]} -eq 0 ]]; then
        echo "Error: No T1w images found for re-run subjects in $RERUN_FILE" >&2
        exit 1
    fi
else
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

# Pre-flight: detect optional T2s for discovered T1s and validate they exist/readable
declare -A T2_INDEX=()
declare -a T2_FOUND_LIST=()
declare -a T2_ERRORS=()

# Only attempt T2 autodetection when enabled in config (cross.t2)
CROSS_T2=$(jq -r '.cross.t2 // false' "$CONFIG")
if [[ "$CROSS_T2" == "true" ]]; then
    idx=0
    total=${#T1W_LIST[@]}
    for t1w_img in "${T1W_LIST[@]}"; do
        idx=$((idx+1))
        fname=$(basename "$t1w_img")
        subj=$(echo "$fname" | grep -o 'sub-[^_]*')
        sess=$(echo "$fname" | grep -o 'ses-[^_]*')
        key="${subj}"
        if [[ -n "$sess" ]]; then
            key="${subj}_${sess}"
        fi

        # subject-level dir to limit search (much faster than searching entire BIDS root)
        if [[ -n "$sess" ]]; then
            search_root="$BIDS_DATA/${subj}/${sess}"
        else
            search_root="$BIDS_DATA/${subj}"
        fi
        # If the subject/session dir doesn't exist, skip
        if [[ ! -d "$search_root" ]]; then
            [[ $DEBUG -eq 1 ]] && echo "[DEBUG] Missing dir for $key: $search_root (skipping)"
            continue
        fi

        [[ $DEBUG -eq 1 ]] && echo "[DEBUG] T2 autodetect ($idx/$total): searching $search_root for ${key}_*T2w.nii.gz"

        # find all matching T2s for this subject/session under the subject dir
        mapfile -t found_t2s < <(find "$search_root" -maxdepth 3 -type f -name "${key}_*T2w.nii.gz" 2>/dev/null | sort)
        if [[ ${#found_t2s[@]} -gt 0 ]]; then
            # convert to relative paths under BIDS_DATA
            BIDS_DATA_SLASH="$BIDS_DATA"
            [[ "${BIDS_DATA_SLASH}" != */ ]] && BIDS_DATA_SLASH="${BIDS_DATA_SLASH}/"
            rels=()
            for f in "${found_t2s[@]}"; do
                if [[ -r "$f" ]]; then
                    rel="${f#${BIDS_DATA_SLASH}}"
                    rels+=("$rel")
                    T2_FOUND_LIST+=("$rel")
                else
                    T2_ERRORS+=("$f")
                fi
            done
            # store space-separated relative paths for this key
            T2_INDEX["$key"]="${rels[*]}"
        fi
    done

    if [[ ${#T2_ERRORS[@]} -gt 0 ]]; then
        echo "[VALIDATION] ERROR: The following T2 files were detected but are not readable:"
        for e in "${T2_ERRORS[@]}"; do
            echo "  - $e"
        done
        echo "Please fix permissions or remove the files from the dataset."
        exit 1
    fi

    if [[ ${#T2_FOUND_LIST[@]} -gt 0 ]]; then
        echo "[VALIDATION] Detected ${#T2_FOUND_LIST[@]} T2w image(s) across subjects/sessions."
    else
        echo "[VALIDATION] No T2w images detected for these T1w inputs; no --t2 will be passed to the container."
    fi
else
    echo "[VALIDATION] T2 autodetection disabled in config (cross.t2=false); no --t2 will be used."
fi


# If nohup is set without batch_size, default to sequential processing (batch_size=1)
if [[ $NOHUP -eq 1 && -z "$BATCH_SIZE" ]]; then
    echo "[INFO] --nohup specified without --batch_size, defaulting to sequential processing (batch_size=1)"
    BATCH_SIZE=1
fi

# Batch mode
if [[ -n "$BATCH_SIZE" ]]; then
    echo "[BATCH] Processing subjects in batches of $BATCH_SIZE"
    
    # Create unique log file with timestamp and type
    TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
    LOG_FILE="$OUTPUT_DIR/batch_cross_${TIMESTAMP}.log"
    
    # Process in batches
    total_subjects=${#T1W_LIST[@]}
    for ((i=0; i<total_subjects; i+=BATCH_SIZE)); do
        batch_start=$i
        batch_end=$((i+BATCH_SIZE-1))
        if [[ $batch_end -ge $total_subjects ]]; then
            batch_end=$((total_subjects-1))
        fi
        
        echo "Processing batch $((i/BATCH_SIZE + 1)): subjects $batch_start to $batch_end" >> "$LOG_FILE"
        
        # Process this batch
        batch_pids=()
        for ((j=batch_start; j<=batch_end; j++)); do
            t1w_img="${T1W_LIST[$j]}"
            fname=$(basename "$t1w_img")
            subj=$(echo "$fname" | grep -o 'sub-[^_]*')
            sess=$(echo "$fname" | grep -o 'ses-[^_]*')
            if [[ -n "$sess" ]]; then
                sid="${subj}_${sess}"
            else
                sid="$subj"
            fi
            
            # Build command (use prevalidated T2_INDEX if available)
            extra_opts=$(parse_json_options_cross "$CONFIG")
            if [[ "$CROSS_T2" == "true" ]]; then
                key="${subj}"
                if [[ -n "$sess" ]]; then
                    key="${subj}_${sess}"
                fi
                t2_list="${T2_INDEX[$key]:-}"
                if [[ -n "$t2_list" ]]; then
                    # take first validated T2 for now
                    read -r first_t2 _ <<< "$t2_list"
                    extra_opts="$extra_opts --t2 /data/${first_t2}"
                    [[ $DEBUG -eq 1 ]] && echo "[DEBUG] Using prevalidated T2 for $key: $first_t2" >> "$LOG_FILE"
                else
                    [[ $DEBUG -eq 1 ]] && echo "[DEBUG] No prevalidated T2 for $key" >> "$LOG_FILE"
                fi
            fi
            
            BIDS_DATA_SLASH="$BIDS_DATA"
            [[ "${BIDS_DATA_SLASH}" != */ ]] && BIDS_DATA_SLASH="${BIDS_DATA_SLASH}/"
            t1w_relpath="${t1w_img#${BIDS_DATA_SLASH}}"
            
            cmd=(singularity exec --nv --no-home -B "$BIDS_DATA":/data -B "$OUTPUT_DIR":/output -B "$LICENSE_DIR":/fs_license "$SIF_FILE" /fastsurfer/run_fastsurfer.sh --t1 "/data/$t1w_relpath" --sid "$sid" --sd /output --fs_license /fs_license/license.txt --3T)
            if [[ -n "$extra_opts" ]]; then
                extra_opts_arr=($extra_opts)
                cmd+=("${extra_opts_arr[@]}")
            fi
            
            echo "Starting $sid in batch..." >> "$LOG_FILE"
            nohup "${cmd[@]}" >> "$LOG_FILE" 2>&1 &
            batch_pids+=($!)
        done
        
        # Wait for this batch to complete
        echo "Waiting for batch to complete..." >> "$LOG_FILE"
        wait_time=0
        max_wait=3600
        
        while ps aux | grep -q "singularity.*fastsurfer" | grep -v grep; do
            sleep 30
            wait_time=$((wait_time + 30))
            running=$(ps aux | grep "singularity.*fastsurfer" | grep -v grep | wc -l)
            echo "Still running: $running FastSurfer processes (waited ${wait_time}s)" >> "$LOG_FILE"
            
            if [[ $running -eq 0 ]]; then
                break
            fi
            
            if [[ $wait_time -ge $max_wait ]]; then
                echo "WARNING: Timeout reached (${max_wait}s), proceeding to next batch..." >> "$LOG_FILE"
                break
            fi
        done
        
        echo "Batch completed after ${wait_time} seconds" >> "$LOG_FILE"
        
        # Wait between batches
        if [[ $((i+BATCH_SIZE)) -lt $total_subjects ]]; then
            echo "Waiting 5 minutes before next batch..." >> "$LOG_FILE"
            sleep 300
        fi
    done
    
    echo "All batches completed! Monitor with: tail -f $LOG_FILE"
    exit 0
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

    # Build options from JSON config and consult prevalidated T2_INDEX
    extra_opts=$(parse_json_options_cross "$CONFIG")
    if [[ $DEBUG -eq 1 ]]; then
        echo "DEBUG: Extra options from JSON: $extra_opts"
    fi
    if [[ "$CROSS_T2" == "true" ]]; then
        key="${subj}"
        if [[ -n "$sess" ]]; then
            key="${subj}_${sess}"
        fi
        t2_list="${T2_INDEX[$key]:-}"
        if [[ -n "$t2_list" ]]; then
            read -r first_t2 _ <<< "$t2_list"
            extra_opts="$extra_opts --t2 /data/${first_t2}"
            [[ $DEBUG -eq 1 ]] && echo "[DEBUG] Using prevalidated T2 for $key: $first_t2"
        else
            [[ $DEBUG -eq 1 ]] && echo "[DEBUG] No prevalidated T2 for $key"
        fi
    fi

    # Always set --sid, --t1, --sd, --fs_license, --3T
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
    elif [[ $NOHUP -eq 1 ]]; then
        log_file="$OUTPUT_DIR/fastsurfer_${sid}.log"
        echo "Running with nohup, output redirected to: $log_file"
        nohup "${cmd[@]}" > "$log_file" 2>&1 &
        pid=$!
        echo "Started process PID: $pid"

        # Sequential by default (batch_size=1 or empty)
        if [[ -z "$BATCH_SIZE" || "$BATCH_SIZE" -le 1 ]]; then
            wait $pid
        else
            # Batching: wait if we have reached the batch size limit
            while [[ $(jobs -rp | wc -l) -ge "$BATCH_SIZE" ]]; do
                sleep 10
            done
        fi
    else
        "${cmd[@]}"
    fi
done