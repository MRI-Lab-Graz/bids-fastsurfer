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
    if [[ $ses_count -ge 2 ]]; then
        # Build list of expected timepoint IDs from BIDS (sub-XXX_ses-YYY)
        mapfile -t ses_list < <(find "$BIDS_FOLDER/$subj" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sort)
        expected_tpids=()
        for ses in "${ses_list[@]}"; do
            expected_tpids+=("${subj}_${ses}")
        done

        # Find .long.<template> symlinks in top-level results directory
        mapfile -t long_links < <(find "$RESULTS_FOLDER" -maxdepth 1 -type l -name "*.long.${subj}" -exec basename {} \; 2>/dev/null || true)

        # If number of long symlinks equals number of sessions, assume longitudinal processing completed
        if [[ ${#long_links[@]} -eq ${#expected_tpids[@]} && ${#long_links[@]} -gt 0 ]]; then
            # Further sanity-check: ensure each linked directory has a stats/aseg.stats
            all_have_aseg=1
            for link in "${long_links[@]}"; do
                if [[ ! -f "${RESULTS_FOLDER%/}/${link}/stats/aseg.stats" ]]; then
                    all_have_aseg=0
                    break
                fi
            done
            if [[ $all_have_aseg -eq 1 ]]; then
                completed_long+=("$subj")
            else
                missing_long+=("$subj")
            fi
            continue
        fi

        # If we reached here, longitudinal outputs are incomplete or missing
        missing_long+=("$subj")
        continue
    fi

    # For single-session subjects, check cross-sectional outputs
    if [[ $ses_count -eq 1 ]]; then
        # Expect a timepoint output directory like <sub>_ses-*/
        timepoint_dir=$(find "$RESULTS_FOLDER" -maxdepth 1 -type d -name "${subj}_ses-*" | head -n1 || true)
        if [[ -z "$timepoint_dir" ]]; then
            # No output at all for this session
            missing_cross+=("$subj")
            continue
        fi

        # Check for aseg.stats in that timepoint's stats folder
        if [[ -f "${timepoint_dir}/stats/aseg.stats" ]]; then
            completed_long+=("$subj")
        else
            missing_cross+=("$subj")
        fi
        continue
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

if [[ ${#missing_cross[@]} -gt 0 ]]; then
    printf '%s\n' "${missing_cross[@]}" | jq -R -s 'split("\n") | map(select(. != "")) | {"subjects": .}' > missing_cross_subjects.json
    echo "Wrote missing cross-sectional subjects to missing_cross_subjects.json"
fi

if [[ ${#missing_long[@]} -gt 0 ]]; then
    printf '%s\n' "${missing_long[@]}" | jq -R -s 'split("\n") | map(select(. != "")) | {"subjects": .}' > missing_long_subjects.json
    echo "Wrote missing longitudinal subjects to missing_long_subjects.json"
fi

# Summary
total_subjects=${#input_subjects[@]}
cross_count=${#single_session[@]}
no_session_count=${#no_session[@]}
missing_cross_count=${#missing_cross[@]}
missing_long_count=${#missing_long[@]}
long_count=$(( total_subjects - cross_count - no_session_count ))
echo
echo "total: $total_subjects"
echo "cross: $cross_count"
echo "missing-cross: $missing_cross_count"
echo "long: $long_count"
echo "missing-long: $missing_long_count"