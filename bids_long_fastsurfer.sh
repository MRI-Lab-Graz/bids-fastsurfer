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
  bash bids_long_fastsurfer.sh <BIDS_ROOT> <OUTPUT_DIR> -c <config.json> --re-run <subjects.json> [OPTIONS]

Required:
  BIDS_ROOT            Path to BIDS dataset root (must contain sub-*/).
  OUTPUT_DIR           Output directory (will be bind-mounted at /output).
  -c, --config FILE    JSON config with keys: fs_license, sif_file, and "long" section.
  --tid SUBJECT        Template subject ID (no session). Provide with or without 'sub-' prefix.
  --tpids LIST         One or more longitudinal timepoint IDs of form sub-XXX_ses-YYY (switches to manual mode).
  --tid SUBJECT        Template subject (switches to manual mode when used with --tpids).

Optional:
  --pilot              (Auto mode) Randomly select one eligible longitudinal subject (>=2 sessions) and process only that subject.
  --re-run FILE        JSON file with subjects to re-run. Format: {"subjects": ["sub-001", "sub-002", ...]}
  --batch_size N       Process N subjects in parallel (requires --re-run). If --nohup is used without --batch_size, defaults to 1 (sequential).
  --nohup              Run commands with nohup for long-running jobs (redirects output to log files). 
                       Without --batch_size, defaults to sequential processing (batch_size=1).
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
  - Skips subjects that are already fully processed (all timepoint directories exist in OUTPUT_DIR AND have longitudinal aseg.stats)
  - With --nohup: runs commands in background with output redirected to log files

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

Re-run JSON Structure Example:
{
  "subjects": ["sub-1291056", "sub-1292036", "sub-1292037"]
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

  # Pilot (one random longitudinal subject)
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --pilot --dry_run

  # Re-run specific subjects from JSON file
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --re-run missing_subjects.json --dry_run

  # Re-run with nohup (defaults to sequential processing, batch_size=1)
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --re-run missing_subjects.json --nohup

  # Re-run with parallel processing (6 subjects at a time)
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --re-run missing_subjects.json --nohup --batch_size 6

EOF
}

############################################
# Helper Functions
############################################

# Create .long symlinks for a subject's timepoints
create_long_symlinks() {
  local output_dir="$1"
  local template_subject="$2"
  shift 2
  local tpids=("$@")
  
  if [[ ${#tpids[@]} -eq 0 ]]; then
    return 0
  fi
  
  echo "[INFO] Creating .long symlinks for ${template_subject}..."
  local created=0
  for tpid in "${tpids[@]}"; do
    local tp_dir="${output_dir%/}/${tpid}"
    local long_link="${output_dir%/}/${tpid}.long.${template_subject}"
    
    if [[ ! -d "$tp_dir" ]]; then
      echo "  [WARN] Timepoint directory not found: $tp_dir"
      continue
    fi
    
    if [[ -e "$long_link" ]]; then
      echo "  [SKIP] ${tpid}.long.${template_subject} (already exists)"
    else
      if ln -s "$tp_dir" "$long_link" 2>/dev/null; then
        echo "  [OK] Created: ${tpid}.long.${template_subject} -> ${tpid}"
        created=$((created + 1))
      else
        echo "  [ERROR] Failed to create symlink: $long_link"
      fi
    fi
  done
  
  if [[ $created -gt 0 ]]; then
    echo "[INFO] Created $created .long symlink(s) for ${template_subject}"
  fi
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
PILOT=0
RERUN_FILE=""
NOHUP=0
PYTHON_CMD="python3"
# AUTO default: enabled unless user provides --tid/--tpids
AUTO=1

BATCH_SIZE=""

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
      # collect until next -* or end
      while [[ $# -gt 0 && "$1" != -* ]]; do
        TPIDS+=("$1")
        shift
      done
      AUTO=0
      ;;
    --auto) # retained for backward compatibility; no-op now (auto is default)
      AUTO=1; shift ;;
    --pilot)
      PILOT=1; shift ;;
    --re-run)
      RERUN_FILE="${2:-}"; shift 2 ;;
      --batch_size)
        BATCH_SIZE="${2:-}"; shift 2 ;;
    --nohup)
      NOHUP=1; shift ;;
    --dry_run)
      DRY_RUN=1; shift ;;
    --debug)
      DEBUG=1; shift ;;
    --py)
      PYTHON_CMD="${2:-python3}"; shift 2 ;;
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
  echo "Error: Missing required positional arguments. Expected <BIDS_ROOT> and <OUTPUT_DIR>." >&2
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
  echo "[INFO] Creating output directory: ${OUTPUT_DIR}"
  mkdir -p "${OUTPUT_DIR}"
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

# --pilot only valid in auto mode
if [[ $PILOT -eq 1 && $AUTO -eq 0 ]]; then
  echo "Error: --pilot can only be used in auto mode (omit --tid/--tpids)." >&2
  exit 1
fi

# --re-run validation
if [[ -n "${RERUN_FILE}" ]]; then
  if [[ ! -f "${RERUN_FILE}" ]]; then
    echo "Error: Re-run file '${RERUN_FILE}' does not exist." >&2
    exit 1
  fi
  if [[ $AUTO -eq 0 ]]; then
    echo "Error: --re-run can only be used in auto mode (omit --tid/--tpids)." >&2
    exit 1
  fi
  AUTO=2  # Special mode for re-run

  # If nohup is set without batch_size, default to sequential processing (batch_size=1)
  if [[ $NOHUP -eq 1 && -z "$BATCH_SIZE" ]]; then
    echo "[INFO] --nohup specified without --batch_size, defaulting to sequential processing (batch_size=1)"
    BATCH_SIZE=1
  fi
  
  # If batch_size is set, trigger batch_fastsurfer.sh and exit
  if [[ -n "$BATCH_SIZE" ]]; then
    echo "[BATCH] Triggering batch_fastsurfer.sh with batch size $BATCH_SIZE"
    script_dir="$(dirname "$0")"
    nohup "$script_dir/batch_fastsurfer.sh" "$BATCH_SIZE" "$RERUN_FILE" > "$OUTPUT_DIR/batch_processing.log" 2>&1 &
    echo "Batch processing started in background. Monitor with: tail -f $OUTPUT_DIR/batch_processing.log"
    exit 0
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
# Collect Longitudinal Options (from JSON long section)
############################################
VALID_LONG_KEYS=(
  parallel parallel_seg parallel_surf
  reg_mode qc_snap surf_only 3T device viewagg_device
  threads threads_seg threads_surf batch ignore_fs_version
  fstess fsqsphere fsaparc no_fs_T1 no_surfreg allow_root base
)
declare -a LONG_OPTS=()
for key in "${VALID_LONG_KEYS[@]}"; do
  raw=$(jq -r --arg k "$key" '.long[$k]' "${CONFIG}" 2>/dev/null || echo "null")
  [[ "$raw" == "null" || -z "$raw" ]] && continue
  type=$(jq -r --arg k "$key" 'if (.long[$k]|type) then (.long[$k]|type) else "null" end' "${CONFIG}")
  if [[ "$type" == "boolean" ]]; then
    [[ "$raw" == "true" ]] && LONG_OPTS+=("--$key")
  else
    LONG_OPTS+=("--$key" "$raw")
  fi
done
if [[ $DEBUG -eq 1 ]]; then
  mapfile -t present_keys < <(jq -r '.long | keys[]' "${CONFIG}")
  for k in "${present_keys[@]}"; do
    if ! printf '%s\n' "${VALID_LONG_KEYS[@]}" | grep -q "^${k}$"; then
      echo "[DEBUG] Ignoring unknown long option key '$k'"
    fi
  done
fi

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
  
  # Check if subject is already fully processed (all timepoint directories exist AND have longitudinal stats)
  subject_fully_processed=1
  for tpid in "${TPIDS[@]}"; do
    expected_dir="${OUTPUT_DIR%/}/${tpid}"
    long_dir="${OUTPUT_DIR%/}/${tpid}.long.${TEMPLATE_SUBJECT}"
    aseg_stats="${long_dir}/stats/aseg.stats"
    if [[ ! -d "$expected_dir" || ! -f "$aseg_stats" ]]; then
      subject_fully_processed=0
      break
    fi
  done
  
  if [[ $subject_fully_processed -eq 1 ]]; then
    echo "[SKIP] ${TEMPLATE_SUBJECT} already fully processed (${#TPIDS[@]} sessions) - skipping"
    exit 0
  fi
  
  cmd=( singularity exec --nv --no-home -B "${BIDS_ROOT%/}":/data -B "${OUTPUT_DIR%/}":/output -B "${LICENSE_DIR}":/fs_license "${SIF_FILE}" /fastsurfer/long_fastsurfer.sh --tid "${TEMPLATE_SUBJECT}" --t1s "${T1_PATHS[@]}" --tpids "${TPIDS[@]}" --sd /output --fs_license /fs_license/license.txt --py "${PYTHON_CMD}" )
  if [[ ${#LONG_OPTS[@]} -gt 0 ]]; then cmd+=( "${LONG_OPTS[@]}" ); fi
  if [[ $DEBUG -eq 1 ]]; then
    echo "[DEBUG] TEMPLATE_SUBJECT: ${TEMPLATE_SUBJECT}"; echo "[DEBUG] TPIDS: ${TPIDS[*]}"; echo "[DEBUG] T1_PATHS (container): ${T1_PATHS[*]}"; echo "[DEBUG] LONG_OPTS: ${LONG_OPTS[*]}"; echo "[DEBUG] Singularity image: ${SIF_FILE}"; echo "[DEBUG] License: ${FS_LICENSE}"; fi
  echo "Running longitudinal FastSurfer for template '${TEMPLATE_SUBJECT}' with ${#TPIDS[@]} timepoints."; echo "Command:"; printf ' %q' "${cmd[@]}"; echo
  if [[ $DRY_RUN -eq 1 ]]; then 
    if [[ $NOHUP -eq 1 ]]; then
      echo "[DRY RUN] Would run with nohup, output to: ${OUTPUT_DIR%/}/long_fastsurfer_${TEMPLATE_SUBJECT}.log"
    else
      echo "[DRY RUN] Would run directly"
    fi
    echo "[DRY RUN] Not executing."; exit 0; 
  fi
  
  if [[ $NOHUP -eq 1 ]]; then
    log_file="${OUTPUT_DIR%/}/long_fastsurfer_${TEMPLATE_SUBJECT}.log"
    echo "Running with nohup, output redirected to: $log_file"
    
    # Use a subshell to run the command and then do the monitoring/cleanup
    # This avoids the "wait: pid is not a child of this shell" error
    (
      if nohup "${cmd[@]}" > "$log_file" 2>&1; then
        # Success: create symlinks (inherited function)
        create_long_symlinks "${OUTPUT_DIR}" "${TEMPLATE_SUBJECT}" "${TPIDS[@]}"
      else
        rc=$?
        ERROR_LOG="${OUTPUT_DIR%/}/fastsurfer_errors.log"
        echo "$(date) [ERROR] ${TEMPLATE_SUBJECT} exited with code $rc (nohup)" >> "$ERROR_LOG"
        echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
        echo "--- tail of ${log_file} ---" >> "$ERROR_LOG"
        tail -n 300 "$log_file" >> "$ERROR_LOG" 2>/dev/null
        echo "--- end tail ---" >> "$ERROR_LOG"
      fi
    ) &
    pid=$!
    echo "Started process PID: $pid"
    echo "Monitor progress with: tail -f $log_file"
  else
    if "${cmd[@]}"; then
      # Create .long symlinks after successful processing
      create_long_symlinks "${OUTPUT_DIR}" "${TEMPLATE_SUBJECT}" "${TPIDS[@]}"
    else
      rc=$?
      ERROR_LOG="${OUTPUT_DIR%/}/fastsurfer_errors.log"
      echo "$(date) [ERROR] ${TEMPLATE_SUBJECT} exited with code $rc (foreground)" >> "$ERROR_LOG"
      echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
      echo "Nohup/direct log not available; consider re-running with --nohup to capture logs." >> "$ERROR_LOG"
      exit $rc
    fi
  fi
else
  # Auto mode (with optional pilot or re-run)
  if [[ $AUTO -eq 2 ]]; then
    # Re-run mode: read subjects from JSON
    echo "[RE-RUN] Reading subjects from ${RERUN_FILE}"
    if ! command -v jq >/dev/null 2>&1; then
      echo "Error: 'jq' is required for --re-run but not found in PATH." >&2
      exit 1
    fi
    
    # Parse JSON and get subjects array
    declare -a RERUN_SUBJECTS=()
    while IFS= read -r subject; do
      RERUN_SUBJECTS+=("${BIDS_ROOT%/}/$subject")
    done < <(jq -r '.subjects[]' "${RERUN_FILE}")
    
    if [[ ${#RERUN_SUBJECTS[@]} -eq 0 ]]; then
      echo "Error: No subjects found in ${RERUN_FILE}" >&2
      exit 1
    fi
    
    echo "[RE-RUN] Found ${#RERUN_SUBJECTS[@]} subjects to re-run"
    
    # Apply --pilot to re-run mode if specified
    if [[ $PILOT -eq 1 ]]; then
      pick_idx=$(( RANDOM % ${#RERUN_SUBJECTS[@]} ))
      echo "[RE-RUN][PILOT] Selected $(basename "${RERUN_SUBJECTS[$pick_idx]}") from ${#RERUN_SUBJECTS[@]} subjects to re-run"
      subjects=("${RERUN_SUBJECTS[$pick_idx]}")
    else
      subjects=("${RERUN_SUBJECTS[@]}")
    fi
  else
    # Standard auto mode
    shopt -s nullglob
    all_subj=("${BIDS_ROOT%/}"/sub-*)
    if [[ ${#all_subj[@]} -eq 0 ]]; then echo "[AUTO] No subjects found."; exit 1; fi
    eligible=()
    for sp in "${all_subj[@]}"; do
      [[ -d "$sp" ]] || continue
      sbase=$(basename "$sp")
      mapfile -t ses_list < <(find "$sp" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sort)
      if [[ ${#ses_list[@]} -ge 2 ]]; then
        eligible+=("$sp")
      else
        [[ $DEBUG -eq 1 ]] && echo "[AUTO][SKIP] $sbase has <2 sessions"
      fi
    done
    if [[ ${#eligible[@]} -eq 0 ]]; then echo "[AUTO] No longitudinal subjects (>=2 sessions)"; exit 1; fi
    if [[ $PILOT -eq 1 ]]; then
      pick_idx=$(( RANDOM % ${#eligible[@]} ))
      echo "[AUTO][PILOT] Selected $(basename "${eligible[$pick_idx]}") from ${#eligible[@]} eligible subjects"
      subjects=("${eligible[$pick_idx]}")
    else
      echo "[AUTO] Processing ${#eligible[@]} eligible subjects"
      subjects=("${eligible[@]}")
    fi
  fi
  total_processed=0
  declare -a pids=()
  for sp in "${subjects[@]}"; do
    sbase=$(basename "$sp")
    mapfile -t ses_list < <(find "$sp" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sort)
    TPIDS_LOCAL=()
    T1_PATHS_LOCAL=()
    skip=0
    for ses in "${ses_list[@]}"; do
      anat_dir="$sp/$ses/anat"
      if [[ ! -d "$anat_dir" ]]; then echo "[AUTO][WARN] Missing anat for $sbase $ses -> skip subject"; skip=1; break; fi
      t1=$(ls -1 "$anat_dir"/*_T1w.nii.gz 2>/dev/null | head -n1 || true)
      [[ -z "$t1" ]] && t1=$(ls -1 "$anat_dir"/*_T1w.nii 2>/dev/null | head -n1 || true)
      if [[ -z "$t1" ]]; then echo "[AUTO][WARN] No T1w for $sbase $ses -> skip subject"; skip=1; break; fi
      rel="${t1#${BIDS_ROOT%/}/}"; container_t1="/data/${rel}"
      TPIDS_LOCAL+=("${sbase}_${ses}")
      T1_PATHS_LOCAL+=("${container_t1}")
    done
    [[ $skip -eq 1 ]] && continue
    if [[ ${#TPIDS_LOCAL[@]} -lt 2 ]]; then [[ $DEBUG -eq 1 ]] && echo "[AUTO][SKIP] $sbase insufficient valid sessions"; continue; fi
    
    # Check if subject is already fully processed (all timepoint directories exist AND have longitudinal stats)
    subject_fully_processed=1
    for tpid in "${TPIDS_LOCAL[@]}"; do
      expected_dir="${OUTPUT_DIR%/}/${tpid}"
      long_dir="${OUTPUT_DIR%/}/${tpid}.long.${sbase}"
      aseg_stats="${long_dir}/stats/aseg.stats"
      if [[ ! -d "$expected_dir" || ! -f "$aseg_stats" ]]; then
        subject_fully_processed=0
        break
      fi
    done
    
    if [[ $subject_fully_processed -eq 1 ]]; then
      echo "[AUTO][SKIP] $sbase already fully processed (${#TPIDS_LOCAL[@]} sessions) - skipping"
      continue
    fi
    
    cmd=( singularity exec --nv --no-home -B "${BIDS_ROOT%/}":/data -B "${OUTPUT_DIR%/}":/output -B "${LICENSE_DIR}":/fs_license "${SIF_FILE}" /fastsurfer/long_fastsurfer.sh --tid "$sbase" --t1s "${T1_PATHS_LOCAL[@]}" --tpids "${TPIDS_LOCAL[@]}" --sd /output --fs_license /fs_license/license.txt --py "${PYTHON_CMD}" )
    if [[ ${#LONG_OPTS[@]} -gt 0 ]]; then cmd+=("${LONG_OPTS[@]}"); fi
    echo "[AUTO] Subject $sbase (${#TPIDS_LOCAL[@]} sessions)"
    printf '  CMD:'; printf ' %q' "${cmd[@]}"; echo
    if [[ $DRY_RUN -eq 1 ]]; then 
      if [[ $NOHUP -eq 1 ]]; then
        echo "  [DRY RUN] Would run with nohup, output to: ${OUTPUT_DIR%/}/long_fastsurfer_${sbase}.log"
      else
        echo "  [DRY RUN] Would run directly"
      fi
      continue; 
    fi
    
    if [[ $NOHUP -eq 1 ]]; then
      log_file="${OUTPUT_DIR%/}/long_fastsurfer_${sbase}.log"
      echo "  Running with nohup, output redirected to: $log_file"
      
      # Use a subshell to run the command and then do the monitoring/cleanup
      # This avoids the "wait: pid is not a child of this shell" error
      (
        if nohup "${cmd[@]}" > "$log_file" 2>&1; then
          # Success: create symlinks (inherited function)
          create_long_symlinks "${OUTPUT_DIR}" "${sbase}" "${TPIDS_LOCAL[@]}"
        else
          rc=$?
          ERROR_LOG="${OUTPUT_DIR%/}/fastsurfer_errors.log"
          echo "$(date) [ERROR] ${sbase} exited with code $rc (nohup)" >> "$ERROR_LOG"
          echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
          echo "--- tail of ${log_file} ---" >> "$ERROR_LOG"
          tail -n 300 "$log_file" >> "$ERROR_LOG" 2>/dev/null
          echo "--- end tail ---" >> "$ERROR_LOG"
        fi
      ) &
      pid=$!
      echo "  Started process PID: $pid"
      pids+=($pid)
      total_processed=$((total_processed+1))
      
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
      if "${cmd[@]}"; then
        total_processed=$((total_processed+1))
        # Create .long symlinks after successful processing
        create_long_symlinks "${OUTPUT_DIR}" "${sbase}" "${TPIDS_LOCAL[@]}"
      else
        rc=$?
        ERROR_LOG="${OUTPUT_DIR%/}/fastsurfer_errors.log"
        echo "$(date) [ERROR] ${sbase} exited with code $rc (foreground)" >> "$ERROR_LOG"
        echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
        echo "Nohup/direct log not available; consider re-running with --nohup to capture logs." >> "$ERROR_LOG"
      fi
    fi
  done
  
  if [[ $NOHUP -eq 1 && $DRY_RUN -eq 0 ]]; then
    echo "[AUTO] Started ${#pids[@]} background processes"
    echo "Monitor individual logs with: tail -f ${OUTPUT_DIR%/}/long_fastsurfer_*.log"
    echo "Check running processes with: ps -p ${pids[*]}"
    echo "PIDs: ${pids[*]}"
    echo ""
    echo "[INFO] .long symlinks will NOT be created automatically for background jobs."
    echo "       After jobs complete, run the repair script:"
    echo "       bash scripts/create_missing_long_symlinks.sh ${OUTPUT_DIR}"
  else
    if [[ $DRY_RUN -eq 1 ]]; then
      echo "[AUTO][DRY RUN] Done."
    else
      echo "[AUTO] Processed $total_processed subject(s)."
    fi
  fi
fi

exit 0