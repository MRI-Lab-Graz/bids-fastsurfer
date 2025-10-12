#!/bin/bash
set -euo pipefail

# Batch processing script for FastSurfer re-runs
# Usage: ./batch_fastsurfer.sh <batch_size> <json_file>

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <batch_size> <json_file>"
    echo "Example: $0 5 missing_subjects.json"
    exit 1
fi

BATCH_SIZE=$1
JSON_FILE=$2

# Extract subjects from JSON
SUBJECTS=($(jq -r '.subjects[]' "$JSON_FILE"))
TOTAL_SUBJECTS=${#SUBJECTS[@]}

echo "Found $TOTAL_SUBJECTS subjects to process in batches of $BATCH_SIZE"

# Process in batches
for ((i=0; i<TOTAL_SUBJECTS; i+=BATCH_SIZE)); do
    BATCH_START=$i
    BATCH_END=$((i+BATCH_SIZE-1))
    if [[ $BATCH_END -ge $TOTAL_SUBJECTS ]]; then
        BATCH_END=$((TOTAL_SUBJECTS-1))
    fi

    echo "Processing batch $((i/BATCH_SIZE + 1)): subjects $BATCH_START to $BATCH_END"

    # Create batch JSON
    BATCH_FILE="batch_${BATCH_START}_${BATCH_END}.json"
    BATCH_SUBJECTS=("${SUBJECTS[@]:BATCH_START:BATCH_SIZE}")
    jq -n --argjson subjects "$(printf '%s\n' "${BATCH_SUBJECTS[@]}" | jq -R . | jq -s .)" \
       '{"subjects": $subjects}' > "$BATCH_FILE"

    echo "Created $BATCH_FILE with ${#BATCH_SUBJECTS[@]} subjects"

    # Run the batch
    echo "Starting batch processing..."
    bash bids_long_fastsurfer.sh /data/mrivault/_0_STAGING/129_PK01/rawdata/ /data/local/129_PK01/derivatives/fastsurfer/ -c fastsurfer_options.json --re-run "$BATCH_FILE" --nohup

    # Wait for batch to complete
    echo "Waiting for batch to complete..."
    while ps aux | grep -q "long_fastsurfer.sh"; do
        sleep 60
        RUNNING=$(ps aux | grep "long_fastsurfer.sh" | grep -v grep | wc -l)
        echo "Still running: $RUNNING processes"
    done

    echo "Batch completed. Cleaning up $BATCH_FILE"
    rm "$BATCH_FILE"

    # Optional: wait between batches
    if [[ $((i+BATCH_SIZE)) -lt $TOTAL_SUBJECTS ]]; then
        echo "Waiting 5 minutes before next batch..."
        sleep 300
    fi
done

echo "All batches completed!"