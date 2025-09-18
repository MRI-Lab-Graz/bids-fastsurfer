#!/usr/bin/env bash
set -euo pipefail

############################################
# Usage / Help
############################################
usage() {
  cat <<'EOF'
Usage:
  bash bids_long_fastsurfer.sh <BIDS_ROOT> <OUTPUT_DIR> -c <config.json> [OPTIONS]
    (Default: auto-detect all longitudinal subjects with >=2 sessions)
  bash bids_long_fastsurfer.sh <BIDS_ROOT> <OUTPUT_DIR> -c <config.json> --tid <subject> --tpids <sub-XX_ses-YY> [<sub-XX_ses-ZZ> ...] [OPTIONS]

Required:
  BIDS_ROOT            Path to BIDS dataset root (must contain sub-*/).
  OUTPUT_DIR           Output directory (will be bind-mounted at /output).
  -c, --config FILE    JSON config with keys: fs_license, sif_file, and "long" section.
  --tid SUBJECT        Template subject ID (no session). Provide with or without 'sub-' prefix.
  --tpids LIST         One or more longitudinal timepoint IDs of form sub-XXX_ses-YYY (switches to manual mode).
  --tid SUBJECT        Template subject (switches to manual mode when used with --tpids).

Optional:
  --dry_run            Print the Singularity command only.
  --debug              Verbose internal debug output.

Behavior:
  - For each --tpids entry sub-XXX_ses-YYY the script locates:
      BIDS_ROOT/sub-XXX/ses-YYY/anat/*_T1w.nii.gz (or .nii)  (first match)
  - Builds: long_fastsurfer.sh --tid <tid> --t1s <T1 list> --tpids <TPID list>
  - Adds additional FastSurfer longitudinal options from the "long" JSON section
    (boolean true => --flag, value => --flag value).
  - Always binds:
       BIDS_ROOT -> /data
       OUTPUT_DIR -> /output
       license parent dir -> /fs_license
  - License passed as: --fs_license /fs_license/license.txt

JSON Structure Example:
{
  "fs_license": "/path/to/license.txt",
  "sif_file": "/path/to/fastsurfer-gpu.sif",
  "long": {
    "parallel": null,
    "parallel_seg": null,
    "parallel_surf": null,
    "reg_mode": "coreg",
    "qc_snap": false,
    "surf_only": false,
    "3T": true
  }
}

Notes:
  - Unknown keys inside "long" are ignored with a warning (debug mode only).
  - If no T1w is found for a timepoint, the script exits with error.
  - --tid is normalized to include sub- prefix internally if missing.

Examples:
  # Manual specification
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --tid sub-001 \
    --tpids sub-001_ses-01 sub-001_ses-02 --dry_run --debug

  # Automatic detection across all subjects with >=2 sessions (default)
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --dry_run

EOF
}

############################################
# Defaults / Vars
############################################
BIDS_ROOT=""
OUTPUT_DIR=""
CONFIG=""
TEMPLATE_SUBJECT=""
declare -a TPIDS=()
DRY_RUN=0
DEBUG=0
# AUTO default: enabled unless user provides --tid/--tpids
AUTO=1

############################################
# Argument Parsing
############################################
if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

POS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -c|--config)
      CONFIG="${2:-}"; shift 2 ;;
    --tid)
      TEMPLATE_SUBJECT="${2:-}"; AUTO=0; shift 2 ;;
    --tpids)
      shift
      # collect until next --* or end
      while [[ $# -gt 0 && "$1" != --* ]]; do
        TPIDS+=("$1")
        shift
      done
      AUTO=0
      ;;
    --auto) # retained for backward compatibility; no-op now (auto is default)
      AUTO=1; shift ;;
    --dry_run)
      DRY_RUN=1; shift ;;
    --debug)
      DEBUG=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    --*)
      echo "Error: Unknown option: $1" >&2
      usage; exit 1 ;;
    *)
      POS+=("$1"); shift ;;
  esac
done

# Positional expectation: BIDS_ROOT OUTPUT_DIR
if [[ ${#POS[@]} -lt 2 ]]; then
  echo "Error: Missing required positional arguments." >&2
  usage
  exit 1
fi
BIDS_ROOT="${POS[0]}"
OUTPUT_DIR="${POS[1]}"

############################################
# Basic Validations
############################################
if [[ -z "${CONFIG}" ]]; then
  echo "Error: Config file not specified. Use -c <config.json>" >&2
  exit 1
fi
if [[ ! -f "${CONFIG}" ]]; then
  echo "Error: Config file '${CONFIG}' not found." >&2
  exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
  echo "Error: 'jq' is required but not found in PATH." >&2
  exit 1
fi
if ! command -v singularity >/dev/null 2>&1; then
  echo "Error: 'singularity' not found in PATH." >&2
  exit 1
fi

if [[ ! -d "${BIDS_ROOT}" ]]; then
  echo "Error: BIDS root '${BIDS_ROOT}' does not exist or is not a directory." >&2
  exit 1
fi
if [[ ! -d "${OUTPUT_DIR}" ]]; then
  echo "Error: Output directory '${OUTPUT_DIR}' does not exist." >&2
  exit 1
fi

if [[ $AUTO -eq 0 ]]; then
  if [[ -z "${TEMPLATE_SUBJECT}" ]]; then
    echo "Error: Manual mode detected ( --tpids used ); --tid <subject> is required." >&2
    exit 1
  fi
  if [[ "${TEMPLATE_SUBJECT}" != sub-* ]]; then
    TEMPLATE_SUBJECT="sub-${TEMPLATE_SUBJECT}"
  fi
  if [[ ${#TPIDS[@]} -eq 0 ]]; then
    echo "Error: Manual mode requires at least one --tpids entry (sub-XXX_ses-YYY)." >&2
    exit 1
  fi
fi

############################################
# Extract Top-level Config Values
############################################
SIF_FILE=$(jq -r '.sif_file // empty' "${CONFIG}")
FS_LICENSE=$(jq -r '.fs_license // empty' "${CONFIG}")

if [[ -z "${SIF_FILE}" || ! -f "${SIF_FILE}" ]]; then
  echo "Error: Singularity image file '${SIF_FILE}' does not exist (sif_file in config)." >&2
  exit 1
fi
if [[ -z "${FS_LICENSE}" || ! -f "${FS_LICENSE}" ]]; then
  echo "Error: FreeSurfer license file '${FS_LICENSE}' does not exist (fs_license in config)." >&2
  exit 1
fi

LICENSE_DIR=$(dirname "${FS_LICENSE}")

############################################
if [[ $AUTO -eq 0 ]]; then
  # Manual mode processing (single subject)
  declare -a T1_PATHS=()
  for tpid in "${TPIDS[@]}"; do
    if [[ ! "${tpid}" =~ ^sub-[^_]+_ses-[^_]+$ ]]; then
      echo "Error: TPID '${tpid}' does not match pattern sub-XXX_ses-YYY." >&2
      exit 1
    fi
    subj_part="${tpid%%_ses-*}"         # sub-XXX
    ses_part="ses-${tpid##*_ses-}"      # ses-YYY
    subj_dir="${BIDS_ROOT%/}/${subj_part}"
    ses_dir="${subj_dir}/${ses_part}"
    anat_dir="${ses_dir}/anat"
    if [[ ! -d "${anat_dir}" ]]; then
      echo "Error: anat directory missing for timepoint '${tpid}' at '${anat_dir}'." >&2
      exit 1
    fi
    t1_candidate=$(ls -1 "${anat_dir}"/*_T1w.nii.gz 2>/dev/null | head -n1 || true)
    if [[ -z "${t1_candidate}" ]]; then
      t1_candidate=$(ls -1 "${anat_dir}"/*_T1w.nii 2>/dev/null | head -n1 || true)
    fi
    if [[ -z "${t1_candidate}" ]]; then
      echo "Error: No T1w image found for '${tpid}' in '${anat_dir}'." >&2
      exit 1
    fi
    rel="${t1_candidate#${BIDS_ROOT%/}/}"
    container_t1="/data/${rel}"
    T1_PATHS+=( "${container_t1}" )
    [[ $DEBUG -eq 1 ]] && {
      echo "[DEBUG] TPID=${tpid}"; echo "[DEBUG] Host T1=${t1_candidate}"; echo "[DEBUG] Container T1=${container_t1}"; }
  done
  cmd=( singularity exec --nv --no-home -B "${BIDS_ROOT%/}":/data -B "${OUTPUT_DIR%/}":/output -B "${LICENSE_DIR}":/fs_license "${SIF_FILE}" /fastsurfer/long_fastsurfer.sh --tid "${TEMPLATE_SUBJECT}" --t1s "${T1_PATHS[@]}" --tpids "${TPIDS[@]}" --sd /output --fs_license /fs_license/license.txt )
  if [[ ${#LONG_OPTS[@]} -gt 0 ]]; then cmd+=( "${LONG_OPTS[@]}" ); fi
  if [[ $DEBUG -eq 1 ]]; then
    echo "[DEBUG] TEMPLATE_SUBJECT: ${TEMPLATE_SUBJECT}"; echo "[DEBUG] TPIDS: ${TPIDS[*]}"; echo "[DEBUG] T1_PATHS (container): ${T1_PATHS[*]}"; echo "[DEBUG] LONG_OPTS: ${LONG_OPTS[*]}"; echo "[DEBUG] Singularity image: ${SIF_FILE}"; echo "[DEBUG] License: ${FS_LICENSE}"; fi
  echo "Running longitudinal FastSurfer for template '${TEMPLATE_SUBJECT}' with ${#TPIDS[@]} timepoints."; echo "Command:"; printf ' %q' "${cmd[@]}"; echo
  if [[ $DRY_RUN -eq 1 ]]; then echo "[DRY RUN] Not executing."; exit 0; fi
  "${cmd[@]}"
else
  # Auto mode: iterate subjects with >=2 sessions
  echo "[AUTO] Default auto-detection: scanning subjects with >=2 sessions..."
  shopt -s nullglob
  subj_dirs=("${BIDS_ROOT%/}"/sub-*)
  if [[ ${#subj_dirs[@]} -eq 0 ]]; then
    echo "[AUTO] No subjects (sub-*) found in ${BIDS_ROOT}"; exit 1
  fi
  total_processed=0
  for subj_path in "${subj_dirs[@]}"; do
    [[ -d "$subj_path" ]] || continue
    subj=$(basename "$subj_path")
    # gather sessions
    mapfile -t sessions < <(find "$subj_path" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sort)
    if [[ ${#sessions[@]} -lt 2 ]]; then
      [[ $DEBUG -eq 1 ]] && echo "[AUTO][SKIP] $subj has <2 sessions"
      continue
    fi
    # Build TPIDs and T1 paths
    TPIDS_LOCAL=()
    T1_PATHS_LOCAL=()
    missing_any=0
    for ses in "${sessions[@]}"; do
      anat_dir="$subj_path/$ses/anat"
      if [[ ! -d "$anat_dir" ]]; then
        echo "[AUTO][WARN] Missing anat dir for $subj $ses; skipping subject."; missing_any=1; break
      fi
      t1_candidate=$(ls -1 "$anat_dir"/*_T1w.nii.gz 2>/dev/null | head -n1 || true)
      if [[ -z "$t1_candidate" ]]; then
        t1_candidate=$(ls -1 "$anat_dir"/*_T1w.nii 2>/dev/null | head -n1 || true)
      fi
      if [[ -z "$t1_candidate" ]]; then
        echo "[AUTO][WARN] No T1w for $subj $ses; skipping subject."; missing_any=1; break
      fi
      rel="${t1_candidate#${BIDS_ROOT%/}/}"
      container_t1="/data/${rel}"
      TPIDS_LOCAL+=("${subj}_${ses}")
      T1_PATHS_LOCAL+=("${container_t1}")
    done
    if [[ $missing_any -eq 1 ]]; then continue; fi
    if [[ ${#TPIDS_LOCAL[@]} -lt 2 ]]; then
      [[ $DEBUG -eq 1 ]] && echo "[AUTO][SKIP] $subj fewer than 2 valid T1 timepoints"
      continue
    fi
    # Build command per subject
    cmd=( singularity exec --nv --no-home -B "${BIDS_ROOT%/}":/data -B "${OUTPUT_DIR%/}":/output -B "${LICENSE_DIR}":/fs_license "${SIF_FILE}" /fastsurfer/long_fastsurfer.sh --tid "$subj" --t1s "${T1_PATHS_LOCAL[@]}" --tpids "${TPIDS_LOCAL[@]}" --sd /output --fs_license /fs_license/license.txt )
    if [[ ${#LONG_OPTS[@]} -gt 0 ]]; then cmd+=( "${LONG_OPTS[@]}" ); fi
    echo "[AUTO] Subject $subj with ${#TPIDS_LOCAL[@]} sessions -> TPIDs: ${TPIDS_LOCAL[*]}"
    printf '  CMD:'; printf ' %q' "${cmd[@]}"; echo
    if [[ $DRY_RUN -eq 1 ]]; then
      continue
    else
      "${cmd[@]}"
      total_processed=$((total_processed+1))
    fi
  done
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[AUTO][DRY RUN] Completed listing commands."
  else
    echo "[AUTO] Processed $total_processed subjects with longitudinal data."
  fi
fi

############################################
# Build Singularity Command
############################################
cmd=( singularity exec --nv --no-home
      -B "${BIDS_ROOT%/}":/data
      -B "${OUTPUT_DIR%/}":/output
      -B "${LICENSE_DIR}":/fs_license
      "${SIF_FILE}"
      /fastsurfer/long_fastsurfer.sh
      --tid "${TEMPLATE_SUBJECT}"
      --t1s "${T1_PATHS[@]}"
      --tpids "${TPIDS[@]}"
      --sd /output
      --fs_license /fs_license/license.txt
)

# Append long opts
if [[ ${#LONG_OPTS[@]} -gt 0 ]]; then
  cmd+=( "${LONG_OPTS[@]}" )
fi

############################################
# Debug / Dry Run / Execute
############################################
if [[ $DEBUG -eq 1 ]]; then
  echo "[DEBUG] TEMPLATE_SUBJECT: ${TEMPLATE_SUBJECT}"
  echo "[DEBUG] TPIDS: ${TPIDS[*]}"
  echo "[DEBUG] T1_PATHS (container): ${T1_PATHS[*]}"
  echo "[DEBUG] LONG_OPTS: ${LONG_OPTS[*]}"
  echo "[DEBUG] Singularity image: ${SIF_FILE}"
  echo "[DEBUG] License: ${FS_LICENSE}"
fi

echo "Running longitudinal FastSurfer for template '${TEMPLATE_SUBJECT}' with ${#TPIDS[@]} timepoints."
echo "Command:"
printf ' %q' "${cmd[@]}"; echo

if [[ $DRY_RUN -eq 1 ]]; then
  echo "[DRY RUN] Not executing."
  exit 0
fi

# Execute
"${cmd[@]}"