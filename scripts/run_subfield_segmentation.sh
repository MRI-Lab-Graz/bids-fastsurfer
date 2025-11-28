#!/bin/bash
set -euo pipefail

# Script to run FreeSurfer subfield segmentations using the Python-based tool (FS 7.3.2+)
# This avoids the need for Matlab Runtime (MCR).
# Usage: ./run_subfield_segmentation.sh <SUBJECTS_DIR> [options] [subject_id]

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <SUBJECTS_DIR> [subject_id] [--hippo-amygdala] [--brainstem] [--thalamus]"
    echo "If subject_id is not provided, runs on all subjects in SUBJECTS_DIR."
    exit 1
fi

SUBJECTS_DIR="$1"
shift
export SUBJECTS_DIR

# Parse arguments
SUBJECT_ID=""
RUN_HIPPO=false
RUN_BRAINSTEM=false
RUN_THALAMUS=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --hippo-amygdala)
      RUN_HIPPO=true
      shift
      ;;
    --brainstem)
      RUN_BRAINSTEM=true
      shift
      ;;
    --thalamus)
      RUN_THALAMUS=true
      shift
      ;;
    -*)
      echo "Unknown option: $1"
      exit 1
      ;;
    *)
      if [[ -z "$SUBJECT_ID" ]]; then
          SUBJECT_ID="$1"
          shift
      else
          echo "Error: Multiple subject IDs provided or unknown argument: $1"
          exit 1
      fi
      ;;
  esac
done

# If no specific structures requested, run all
if [[ "$RUN_HIPPO" == "false" && "$RUN_BRAINSTEM" == "false" && "$RUN_THALAMUS" == "false" ]]; then
    RUN_HIPPO=true
    RUN_BRAINSTEM=true
    RUN_THALAMUS=true
fi

# Check for segment_subregions tool
if ! command -v segment_subregions &> /dev/null; then
    echo "Error: 'segment_subregions' command not found."
    echo "Please ensure you have sourced FreeSurfer 7.3.2 or later."
    exit 1
fi

if [[ -n "$SUBJECT_ID" ]]; then
    SUBJECTS=("$SUBJECT_ID")
else
    # Find all subjects (directories starting with sub-)
    SUBJECTS=($(find "$SUBJECTS_DIR" -maxdepth 1 -type d -name "sub-*" -exec basename {} \; | sort))
fi

if [[ ${#SUBJECTS[@]} -eq 0 ]]; then
    echo "No subjects found in $SUBJECTS_DIR"
    exit 1
fi

echo "Found ${#SUBJECTS[@]} directories to check."

# Function to check and create wmparc.mgz if missing
check_and_create_wmparc() {
    local sub=$1
    local sd=$2
    local mri_dir="$sd/$sub/mri"
    
    if [[ ! -f "$mri_dir/wmparc.mgz" ]]; then
        echo "  wmparc.mgz missing for $sub. Attempting to generate it..."
        
        # 1. Ensure aseg.mgz exists (required by mris_volmask)
        if [[ ! -f "$mri_dir/aseg.mgz" ]]; then
            if [[ -f "$mri_dir/aseg.auto.mgz" ]]; then
                echo "    Linking aseg.auto.mgz -> aseg.mgz"
                ln -s aseg.auto.mgz "$mri_dir/aseg.mgz"
            else
                echo "    Error: aseg.auto.mgz not found. Cannot create aseg.mgz."
                return 1
            fi
        fi
        
        # 2. Ensure ribbon.mgz exists (required by mri_aparc2aseg)
        if [[ ! -f "$mri_dir/ribbon.mgz" ]]; then
            echo "    Generating ribbon.mgz (mris_volmask)..."
            (export SUBJECTS_DIR="$sd"; mris_volmask --save_ribbon "$sub" > "$mri_dir/mris_volmask.log" 2>&1)
            if [[ $? -ne 0 ]]; then
                echo "    Error: mris_volmask failed. Check $mri_dir/mris_volmask.log"
                return 1
            fi
        fi
        
        # 3. Generate wmparc.mgz
        echo "    Generating wmparc.mgz (mri_aparc2aseg)..."
        local annot_flag=""
        if [[ -f "$sd/$sub/label/lh.aparc.DKTatlas.mapped.annot" ]]; then
            annot_flag="--annot aparc.DKTatlas.mapped"
        elif [[ -f "$sd/$sub/label/lh.aparc.mapped.annot" ]]; then
            annot_flag="--annot aparc.mapped"
        fi
        
        (export SUBJECTS_DIR="$sd"; mri_aparc2aseg --s "$sub" --labelwm --hypo-as-wm --rip-unknown --ctxseg "$mri_dir/aparc.DKTatlas+aseg.deep.mgz" --o "$mri_dir/wmparc.mgz" --aseg "$mri_dir/aseg.mgz" $annot_flag > "$mri_dir/mri_aparc2aseg.log" 2>&1)
        if [[ $? -ne 0 ]]; then
            echo "    Error: mri_aparc2aseg failed. Check $mri_dir/mri_aparc2aseg.log"
            return 1
        fi
        echo "    Successfully created wmparc.mgz"
    fi
}

for sub in "${SUBJECTS[@]}"; do
    SUB_DIR="$SUBJECTS_DIR/$sub"
    
    # Check if it is a longitudinal base template (FastSurfer creates base-tps.fastsurfer)
    if [[ -f "$SUB_DIR/base-tps.fastsurfer" ]]; then
        echo "Processing Longitudinal Base: $sub"
        
        # Fix: segment_subregions expects 'base-tps', but FastSurfer might name it 'base-tps.fastsurfer'
        if [[ ! -f "$SUB_DIR/base-tps" ]]; then
            echo "  Creating symlink: base-tps -> base-tps.fastsurfer"
            ln -s base-tps.fastsurfer "$SUB_DIR/base-tps"
        fi

        # Ensure wmparc.mgz exists (needed for hippo-amygdala)
        if [[ "$RUN_HIPPO" == "true" ]]; then
            check_and_create_wmparc "$sub" "$SUBJECTS_DIR"
        fi

        # 1. Hippocampal Subfields and Amygdala
        if [[ "$RUN_HIPPO" == "true" ]]; then
            echo "  [Longitudinal] Running hippo-amygdala..."
            segment_subregions hippo-amygdala --long-base "$sub" --sd "$SUBJECTS_DIR"
        fi

        # 2. Brainstem
        if [[ "$RUN_BRAINSTEM" == "true" ]]; then
            echo "  [Longitudinal] Running brainstem..."
            segment_subregions brainstem --long-base "$sub" --sd "$SUBJECTS_DIR"
        fi

        # 3. Thalamus
        if [[ "$RUN_THALAMUS" == "true" ]]; then
            echo "  [Longitudinal] Running thalamus..."
            segment_subregions thalamus --long-base "$sub" --sd "$SUBJECTS_DIR"
        fi

    elif [[ "$sub" == *".long."* ]]; then
        # This is a longitudinal timepoint directory (e.g., sub-01_ses-01.long.sub-01)
        # These are processed automatically when running on the base, so we skip them here.
        echo "Skipping longitudinal timepoint directory: $sub (processed via base)"
        
    elif [[ "$sub" == *"_ses-"* ]]; then
        # This is a cross-sectional timepoint directory (e.g., sub-01_ses-01)
        # If you want to process these cross-sectionally, uncomment the block below.
        echo "Skipping cross-sectional timepoint: $sub (focusing on longitudinal stream)"
        
    else
        # Assume this is a pure cross-sectional subject (no sessions, no .long.)
        echo "Processing Cross-Sectional Subject: $sub"
        
        # Ensure wmparc.mgz exists (needed for hippo-amygdala)
        if [[ "$RUN_HIPPO" == "true" ]]; then
            check_and_create_wmparc "$sub" "$SUBJECTS_DIR"
        fi
        
        if [[ "$RUN_HIPPO" == "true" ]]; then
            echo "  [Cross-Sectional] Running hippo-amygdala..."
            segment_subregions hippo-amygdala --cross "$sub" --sd "$SUBJECTS_DIR"
        fi

        if [[ "$RUN_BRAINSTEM" == "true" ]]; then
            echo "  [Cross-Sectional] Running brainstem..."
            segment_subregions brainstem --cross "$sub" --sd "$SUBJECTS_DIR"
        fi

        if [[ "$RUN_THALAMUS" == "true" ]]; then
            echo "  [Cross-Sectional] Running thalamus..."
            segment_subregions thalamus --cross "$sub" --sd "$SUBJECTS_DIR"
        fi
    fi
done

echo "All done."
