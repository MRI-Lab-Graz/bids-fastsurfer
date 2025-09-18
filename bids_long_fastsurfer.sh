#!/bin/bash

# Usage: ./bids_long_fastsurfer.sh /bids-folder /outputfolder -c long_fastsurfer_options.json --tid <templateID> --tpids <tpID1> <tpID2> ... --t1s <T1_1> <T1_2> ... [OPTIONS]
#
# Arguments:
#   /bids-folder           Path to BIDS input directory
#   /outputfolder          Path to output directory (subjects_dir)
#   -c <config.json>       Path to JSON config file with long FastSurfer options
#   --tid <templateID>     ID for person-specific template directory
#   --tpids <tpID1> ...    IDs for time points
#   --t1s <T1_1> ...       Absolute paths to T1w images for each time point
#   --py <python_cmd>      Python command (optional)
#   --parallel <n>|max     Parallel processing pool size (optional)
#   --parallel_seg <n>|max Parallel segmentation pool size (optional)
#   --parallel_surf <n>|max Parallel surface pool size (optional)
#   --dry_run              Print the Singularity command instead of running it
#   --debug                Print debug information
#
# Example:
#   ./bids_long_fastsurfer.sh ./data ./output -c long_fastsurfer_options.json --tid sub-01 --tpids tp1 tp2 --t1s /abs/path/tp1.nii.gz /abs/path/tp2.nii.gz --dry_run --debug

set -e

show_help() {
	echo "\nUsage: $0 /bids-folder /outputfolder -c long_fastsurfer_options.json --tid <subject> --tpids <tpID1> <tpID2> ... --t1s <sub-xx ses-yy> <sub-xx ses-yy> ... [OPTIONS]"
	echo "\nArguments:"
	echo "  /bids-folder           Path to BIDS input directory"
	echo "  /outputfolder          Path to output directory (subjects_dir)"
	echo "  -c <config.json>       Path to JSON config file with long FastSurfer options"
	echo "  --tid <subject>        Subject ID for template directory (e.g. sub-01)"
	echo "  --tpids <tpID1> ...    Time point IDs (e.g. sub-01_ses-01 sub-01_ses-02)"
	echo "  --t1s <sub ses> ...    Subject/session pairs for each time point (e.g. sub-01 ses-01 sub-01 ses-02)"
	echo "  --py <python_cmd>      Python command (optional)"
	echo "  --parallel <n>|max     Parallel processing pool size (optional)"
	echo "  --parallel_seg <n>|max Parallel segmentation pool size (optional)"
	echo "  --parallel_surf <n>|max Parallel surface pool size (optional)"
	echo "  --dry_run              Print the Singularity command instead of running it"
	echo "  --debug                Print debug information"
	echo "\nExample:"
	echo "  $0 ./data ./output -c long_fastsurfer_options.json --tid sub-01 --tpids sub-01_ses-01 sub-01_ses-02 --t1s sub-01 ses-01 sub-01 ses-02 --dry_run --debug"
	exit 1
}

# Defaults
DRY_RUN=0
DEBUG=0
PYTHON="python3"
CONFIG=""
BIDS_DATA=""
OUTPUT_DIR=""
TID=""
TPIDS=()
T1S_SUBJ=()
T1S_SES=()
PARALLEL=""
PARALLEL_SEG=""
PARALLEL_SURF=""

# Parse positional and named arguments
while [[ $# -gt 0 ]]; do
	case "$1" in
		-h|--help)
			show_help
			;;
		-c|--config)
			CONFIG="$2"; shift 2;;
		--tid)
			TID="$2"; shift 2;;
		--tpids)
			shift
			while [[ $# -gt 0 && ! $1 =~ ^-- ]]; do
				TPIDS+=("$1")
				shift
			done
			;;
		--t1s)
			shift
			while [[ $# -gt 1 && ! $1 =~ ^-- ]]; do
				T1S_SUBJ+=("$1")
				T1S_SES+=("$2")
				shift 2
			done
			;;
		--py)
			PYTHON="$2"; shift 2;;
		--parallel)
			PARALLEL="$2"; shift 2;;
		--parallel_seg)
			PARALLEL_SEG="$2"; shift 2;;
		--parallel_surf)
			PARALLEL_SURF="$2"; shift 2;;
		--dry_run)
			DRY_RUN=1; shift;;
		--debug)
			DEBUG=1; shift;;
		*)
			if [[ -z "$BIDS_DATA" ]]; then
				BIDS_DATA="$1"
			elif [[ -z "$OUTPUT_DIR" ]]; then
				OUTPUT_DIR="$1"
			else
				echo "Unknown argument: $1"; show_help
			fi
			shift;;
	esac
done

# Validate required arguments
if [[ -z "$BIDS_DATA" || -z "$OUTPUT_DIR" || -z "$CONFIG" || -z "$TID" || ${#TPIDS[@]} -eq 0 ]]; then
	echo "Missing required arguments."; show_help
fi

# For each --tpids entry, find the corresponding T1w image in BIDS
T1S_PATHS=()
for tpid in "${TPIDS[@]}"; do
	subj=$(echo "$tpid" | grep -o 'sub-[^_]*')
	ses=$(echo "$tpid" | grep -o 'ses-[^_]*')
	if [[ -n "$ses" ]]; then
		t1w=$(find "$BIDS_DATA/$subj/$ses/anat" -type f -name "${subj}_${ses}_desc-preproc_T1w.nii.gz" 2>/dev/null | head -n1)
		if [[ -z "$t1w" ]]; then
			t1w=$(find "$BIDS_DATA/$subj/$ses/anat" -type f -name "${subj}_${ses}_T1w.nii.gz" 2>/dev/null | head -n1)
		fi
	else
		t1w=$(find "$BIDS_DATA/$subj/anat" -type f -name "${subj}_desc-preproc_T1w.nii.gz" 2>/dev/null | head -n1)
		if [[ -z "$t1w" ]]; then
			t1w=$(find "$BIDS_DATA/$subj/anat" -type f -name "${subj}_T1w.nii.gz" 2>/dev/null | head -n1)
		fi
	fi
	if [[ -z "$t1w" ]]; then
		echo "Error: No T1w image found for $tpid in $BIDS_DATA."; exit 1
	fi
	T1S_PATHS+=("$t1w")
done

if [[ ${#TPIDS[@]} -ne ${#T1S_PATHS[@]} ]]; then
	echo "Number of --tpids and found T1w images does not match."; exit 1
fi

if [[ $DEBUG -eq 1 ]]; then
	echo "DEBUG: BIDS_DATA=$BIDS_DATA"
	echo "DEBUG: OUTPUT_DIR=$OUTPUT_DIR"
	echo "DEBUG: CONFIG=$CONFIG"
	echo "DEBUG: TID=$TID"
	echo "DEBUG: TPIDS=${TPIDS[*]}"
	echo "DEBUG: T1S_PATHS=${T1S_PATHS[*]}"
	echo "DEBUG: PYTHON=$PYTHON"
	echo "DEBUG: PARALLEL=$PARALLEL"
	echo "DEBUG: PARALLEL_SEG=$PARALLEL_SEG"
	echo "DEBUG: PARALLEL_SURF=$PARALLEL_SURF"
fi
