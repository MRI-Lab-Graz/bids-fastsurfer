#!/bin/bash
# Usage: fastsurfer_check_results.sh <BIDS_FOLDER> <RESULTS_FOLDER>
# Compares subjects in BIDS_FOLDER to completed outputs in RESULTS_FOLDER and lists missing or failed subjects.
# Also checks session counts and suggests processing methods.

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <BIDS_FOLDER> <RESULTS_FOLDER>"
    exit 1
fi

BIDS_FOLDER="$1"
RESULTS_FOLDER="$2"

# Find all subjects in BIDS folder
mapfile -t input_subjects < <(find "$BIDS_FOLDER" -maxdepth 1 -type d -name 'sub-*' -exec basename {} \; | sort)

echo "Checking results for ${#input_subjects[@]} subjects in $BIDS_FOLDER"
echo

missing_long=()
failed_long=()
completed_long=()
single_session=()
no_session=()

for subj in "${input_subjects[@]}"; do
    subj_dir="$BIDS_FOLDER/$subj"
    # Count sessions
    ses_count=$(find "$subj_dir" -maxdepth 1 -type d -name 'ses-*' | wc -l)
    
    if [[ $ses_count -eq 0 ]]; then
        no_session+=("$subj")
        continue
    elif [[ $ses_count -eq 1 ]]; then
        single_session+=("$subj")
        continue
    fi
    
    # For >=2 sessions, check longitudinal results
    # Check if aseg.stats exists in longitudinal output
    aseg_file=$(find "$RESULTS_FOLDER" -maxdepth 4 -type f -name 'aseg.stats' | grep "$subj" | head -n1)
    if [[ -z "$aseg_file" ]]; then
        missing_long+=("$subj")
        continue
    fi

    # Check recon-all.log for successful completion
    log_file="$RESULTS_FOLDER/$subj/scripts/recon-all.log"
    if [[ ! -f "$log_file" ]]; then
        failed_long+=("$subj (log file missing)")
        continue
    fi

    # Get last line of log
    last_line=$(tail -n1 "$log_file")
    if [[ "$last_line" =~ finished\ without\ error ]]; then
        completed_long+=("$subj")
    else
        failed_long+=("$subj (log does not indicate success: $last_line)")
    fi
done

if [[ ${#no_session[@]} -gt 0 ]]; then
    echo "Subjects with no sessions (check data integrity):"
    printf '  %s\n' "${no_session[@]}"
    echo
fi

if [[ ${#single_session[@]} -gt 0 ]]; then
    echo "Subjects with 1 session (use cross-sectional script: bids_faststurfer.sh):"
    printf '  %s\n' "${single_session[@]}"
    echo "  Example: bash bids_faststurfer.sh $BIDS_FOLDER $RESULTS_FOLDER -c fastsurfer_options.json --sub ${single_session[0]}"
    echo
fi

if [[ ${#missing_long[@]} -gt 0 ]]; then
    echo "Missing longitudinal subjects (>=2 sessions, no aseg.stats):"
    printf '  %s\n' "${missing_long[@]}"
    echo "  Example: bash bids_long_fastsurfer.sh $BIDS_FOLDER $RESULTS_FOLDER -c fastsurfer_options.json --tid ${missing_long[0]%%_*} --tpids $(find "$BIDS_FOLDER/${missing_long[0]}" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sed 's/ses-//' | xargs -I {} echo "${missing_long[0]}_ses-{}" | tr '\n' ' ')"
    echo
fi

if [[ ${#failed_long[@]} -gt 0 ]]; then
    echo "Failed longitudinal subjects (aseg.stats exists but log indicates failure):"
    printf '  %s\n' "${failed_long[@]}"
    echo
fi

if [[ ${#completed_long[@]} -gt 0 ]]; then
    echo "Completed longitudinal subjects (${#completed_long[@]}):"
    printf '  %s\n' "${completed_long[@]}"
else
    echo "No longitudinal subjects completed successfully."
fi

# Write JSON files for missing subjects
if [[ ${#single_session[@]} -gt 0 ]]; then
    printf '%s\n' "${single_session[@]}" | jq -R -s 'split("\n") | map(select(. != "")) | {"subjects": .}' > missing_cross_subjects.json
    echo "Wrote missing cross-sectional subjects to missing_cross_subjects.json"
fi

if [[ ${#missing_long[@]} -gt 0 ]]; then
    printf '%s\n' "${missing_long[@]}" | jq -R -s 'split("\n") | map(select(. != "")) | {"subjects": .}' > missing_long_subjects.json
    echo "Wrote missing longitudinal subjects to missing_long_subjects.json"
fi