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

Manual mode (single subject):
  --tid SUBJECT        Template subject ID (no session). Provide with or without 'sub-' prefix.
  --tpids LIST         One or more timepoint IDs of form sub-XXX_ses-YYY (enables manual mode).

Optional:
  --pilot              (Auto mode) Randomly select one eligible longitudinal subject (>=2 sessions).
  --re-run FILE        JSON file with subjects to re-run. Format: {"subjects": ["sub-001", ...]}
  --batch_size N       Process N subjects per batch; wait for each batch to finish before starting
                       the next. Requires --nohup. Without --batch_size, --nohup defaults to
                       sequential (N=1) to avoid GPU memory exhaustion.
  --nohup              Run with nohup for long jobs (output redirected to per-subject log files).
                       Always sequential by default (batch_size=1) to avoid GPU OOM.
  --force              Re-run even if subject is already fully processed. Deletes existing
                       .long.* directories before re-running. The template directory
                       (OUTPUT_DIR/<tid>) is preserved — only .long.* dirs are removed.
                       Safe with --dry_run: nothing is modified in dry-run mode.
  --dry_run            Print the Singularity command only; no files are run or deleted.
  --debug              Verbose internal debug output.

Behavior:
  - Skips subjects already fully processed (all .long.* dirs present with aseg.stats) unless --force.
  - --nohup without --batch_size defaults to sequential (batch_size=1) to avoid GPU OOM.
  - --force with --dry_run only prints what would be deleted; nothing is modified.
  - GPU is auto-detected; --nv is only passed to Singularity when a GPU is available.
  - With --nohup, .long symlinks are NOT created automatically; run after completion:
      bash scripts/create_missing_long_symlinks.sh <OUTPUT_DIR>

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

  # Re-run specific subjects sequentially (safe for single-GPU setups)
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --re-run missing_subjects.json --nohup

  # Re-run with 2 subjects in parallel (only if multiple GPUs are available)
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --re-run missing_subjects.json --nohup --batch_size 2

  # Force re-run — dry run first to verify what would be deleted
  bash bids_long_fastsurfer.sh /data/BIDS /data/derivatives/fastsurfer_long \
    -c fastsurfer_options.json --force --dry_run

EOF
}

############################################
# Helper: create .long symlinks for a subject's timepoints
############################################
create_long_symlinks() {
  local output_dir="$1"
  local template_subject="$2"
  shift 2
  local tpids=("$@")

  [[ ${#tpids[@]} -eq 0 ]] && return 0

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
    elif ln -s "$tp_dir" "$long_link" 2>/dev/null; then
      echo "  [OK] Created: ${tpid}.long.${template_subject} -> ${tpid}"
      created=$((created + 1))
    else
      echo "  [ERROR] Failed to create symlink: $long_link"
    fi
  done
  [[ $created -gt 0 ]] && echo "[INFO] Created $created .long symlink(s) for ${template_subject}"
}

############################################
# Helper: delete (or report) .long.* dirs for --force
# Never deletes in dry-run mode.
############################################
force_delete_long_outputs() {
  local template_subject="$1"
  shift
  local tpids=("$@")

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[DRY RUN][FORCE] Would delete longitudinal output for ${template_subject}:"
    for tpid in "${tpids[@]}"; do
      echo "  rm -rf ${OUTPUT_DIR%/}/${tpid}.long.${template_subject}"
    done
  else
    echo "[FORCE] Removing .long.* directories for ${template_subject}"
    for tpid in "${tpids[@]}"; do
      rm -rf "${OUTPUT_DIR%/}/${tpid}.long.${template_subject}"
    done
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
RERUN_MODE=0
NOHUP=0
FORCE=0
PYTHON_CMD="python3"
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
      while [[ $# -gt 0 && "$1" != --* ]]; do
        TPIDS+=("$1"); shift
      done
      AUTO=0 ;;
    --auto)
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
    --force)
      FORCE=1; shift ;;
    --debug)
      DEBUG=1; shift ;;
    --py)
      PYTHON_CMD="${2:-python3}"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    --*)
      echo "Error: Unknown option: $1" >&2; usage; exit 1 ;;
    *)
      POS+=("$1"); shift ;;
  esac
done

if [[ ${#POS[@]} -lt 2 ]]; then
  echo "Error: Missing required positional arguments." >&2; usage; exit 1
fi
BIDS_ROOT="${POS[0]}"
OUTPUT_DIR="${POS[1]}"

############################################
# Basic Validations
############################################
[[ -z "${CONFIG}" ]]      && { echo "Error: Config file not specified. Use -c <config.json>" >&2; exit 1; }
[[ ! -f "${CONFIG}" ]]    && { echo "Error: Config file '${CONFIG}' not found." >&2; exit 1; }
! command -v jq >/dev/null 2>&1          && { echo "Error: 'jq' is required but not found in PATH." >&2; exit 1; }
! command -v singularity >/dev/null 2>&1 && { echo "Error: 'singularity' not found in PATH." >&2; exit 1; }
[[ ! -d "${BIDS_ROOT}" ]]  && { echo "Error: BIDS root '${BIDS_ROOT}' does not exist." >&2; exit 1; }
if [[ ! -d "${OUTPUT_DIR}" ]]; then
  echo "[INFO] Creating output directory: ${OUTPUT_DIR}"
  mkdir -p "${OUTPUT_DIR}"
fi

if [[ $AUTO -eq 0 ]]; then
  [[ -z "${TEMPLATE_SUBJECT}" ]] && { echo "Error: --tpids used without --tid <subject>." >&2; exit 1; }
  [[ "${TEMPLATE_SUBJECT}" != sub-* ]] && TEMPLATE_SUBJECT="sub-${TEMPLATE_SUBJECT}"
  [[ ${#TPIDS[@]} -eq 0 ]] && { echo "Error: Manual mode requires at least one --tpids entry (sub-XXX_ses-YYY)." >&2; exit 1; }
fi

[[ $PILOT -eq 1 && $AUTO -eq 0 ]] && { echo "Error: --pilot can only be used in auto mode (omit --tid/--tpids)." >&2; exit 1; }

if [[ -n "$BATCH_SIZE" && $NOHUP -eq 0 ]]; then
  echo "Error: --batch_size requires --nohup." >&2; exit 1
fi

# Sequential default to protect single-GPU setups
if [[ $NOHUP -eq 1 && -z "$BATCH_SIZE" ]]; then
  echo "[INFO] --nohup without --batch_size: defaulting to sequential (batch_size=1) to protect GPU memory"
  BATCH_SIZE=1
fi

if [[ -n "${RERUN_FILE}" ]]; then
  [[ ! -f "${RERUN_FILE}" ]] && { echo "Error: Re-run file '${RERUN_FILE}' does not exist." >&2; exit 1; }
  [[ $AUTO -eq 0 ]] && { echo "Error: --re-run cannot be used with --tid/--tpids." >&2; exit 1; }
  RERUN_MODE=1
fi

############################################
# Extract Config Values
############################################
SIF_FILE=$(jq -r '.sif_file // empty' "${CONFIG}")
FS_LICENSE=$(jq -r '.fs_license // empty' "${CONFIG}")

[[ -z "${SIF_FILE}" || ! -f "${SIF_FILE}" ]] && {
  echo "Error: SIF file '${SIF_FILE}' not found (sif_file in config)." >&2; exit 1; }
[[ -z "${FS_LICENSE}" || ! -f "${FS_LICENSE}" ]] && {
  echo "Error: License file '${FS_LICENSE}' not found (fs_license in config)." >&2; exit 1; }

LICENSE_DIR=$(dirname "${FS_LICENSE}")

############################################
# GPU Detection — only pass --nv when a GPU is available
############################################
declare -a NV_OPTS=()
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
  NV_OPTS=("--nv")
  [[ $DEBUG -eq 1 ]] && echo "[DEBUG] GPU detected; Singularity will use --nv"
else
  echo "[INFO] No GPU detected; running in CPU mode (--nv omitted)"
fi

############################################
# Collect Long Options from Config
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
    printf '%s\n' "${VALID_LONG_KEYS[@]}" | grep -q "^${k}$" || echo "[DEBUG] Ignoring unknown long key '$k'"
  done
fi

############################################
# Manual Mode
############################################
if [[ $AUTO -eq 0 ]]; then
  declare -a T1_PATHS=()
  for tpid in "${TPIDS[@]}"; do
    [[ ! "${tpid}" =~ ^sub-[^_]+_ses-[^_]+$ ]] && {
      echo "Error: TPID '${tpid}' does not match pattern sub-XXX_ses-YYY." >&2; exit 1; }
    subj_part="${tpid%%_ses-*}"
    ses_part="ses-${tpid##*_ses-}"
    anat_dir="${BIDS_ROOT%/}/${subj_part}/${ses_part}/anat"
    [[ ! -d "${anat_dir}" ]] && {
      echo "Error: anat directory missing for '${tpid}' at '${anat_dir}'." >&2; exit 1; }
    t1=$(ls -1 "${anat_dir}"/*_T1w.nii.gz 2>/dev/null | head -n1 || true)
    [[ -z "$t1" ]] && t1=$(ls -1 "${anat_dir}"/*_T1w.nii 2>/dev/null | head -n1 || true)
    [[ -z "$t1" ]] && {
      echo "Error: No T1w image found for '${tpid}' in '${anat_dir}'." >&2; exit 1; }
    rel="${t1#${BIDS_ROOT%/}/}"
    T1_PATHS+=("/data/${rel}")
    [[ $DEBUG -eq 1 ]] && echo "[DEBUG] TPID=${tpid}  Host=${t1}  Container=/data/${rel}"
  done

  # Check if already fully processed
  subject_fully_processed=1
  for tpid in "${TPIDS[@]}"; do
    long_dir="${OUTPUT_DIR%/}/${tpid}.long.${TEMPLATE_SUBJECT}"
    [[ ! -d "$long_dir" || ! -f "${long_dir}/stats/aseg.stats" ]] && {
      subject_fully_processed=0; break; }
  done

  if [[ $subject_fully_processed -eq 1 ]]; then
    if [[ $FORCE -eq 1 ]]; then
      force_delete_long_outputs "${TEMPLATE_SUBJECT}" "${TPIDS[@]}"
      [[ $DRY_RUN -eq 1 ]] && exit 0
    else
      echo "[SKIP] ${TEMPLATE_SUBJECT} already fully processed (${#TPIDS[@]} sessions) - skipping"
      exit 0
    fi
  fi

  cmd=( singularity exec "${NV_OPTS[@]}" --no-home
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
    --py "${PYTHON_CMD}" )
  [[ ${#LONG_OPTS[@]} -gt 0 ]] && cmd+=("${LONG_OPTS[@]}")

  if [[ $DEBUG -eq 1 ]]; then
    echo "[DEBUG] TEMPLATE_SUBJECT: ${TEMPLATE_SUBJECT}"
    echo "[DEBUG] TPIDS: ${TPIDS[*]}"
    echo "[DEBUG] T1_PATHS (container): ${T1_PATHS[*]}"
    echo "[DEBUG] LONG_OPTS: ${LONG_OPTS[*]:-none}"
    echo "[DEBUG] NV_OPTS: ${NV_OPTS[*]:-none}"
  fi

  echo "Running longitudinal FastSurfer for '${TEMPLATE_SUBJECT}' with ${#TPIDS[@]} timepoints."
  printf 'Command:'; printf ' %q' "${cmd[@]}"; echo

  if [[ $DRY_RUN -eq 1 ]]; then
    [[ $NOHUP -eq 1 ]] && echo "[DRY RUN] Would run with nohup -> ${OUTPUT_DIR%/}/long_fastsurfer_${TEMPLATE_SUBJECT}.log"
    echo "[DRY RUN] Not executing."
    exit 0
  fi

  ERROR_LOG="${OUTPUT_DIR%/}/fastsurfer_errors.log"
  log_file="${OUTPUT_DIR%/}/long_fastsurfer_${TEMPLATE_SUBJECT}.log"

  if [[ $NOHUP -eq 1 ]]; then
    echo "Running with nohup -> $log_file"
    ( if nohup "${cmd[@]}" > "$log_file" 2>&1; then
        create_long_symlinks "${OUTPUT_DIR}" "${TEMPLATE_SUBJECT}" "${TPIDS[@]}"
      else
        rc=$?
        echo "$(date) [ERROR] ${TEMPLATE_SUBJECT} exited with code $rc (nohup)" >> "$ERROR_LOG"
        echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
        echo "--- tail of ${log_file} ---" >> "$ERROR_LOG"
        tail -n 300 "$log_file" >> "$ERROR_LOG" 2>/dev/null
        echo "--- end tail ---" >> "$ERROR_LOG"
      fi ) &
    pid=$!
    echo "Started PID: $pid  (monitor: tail -f $log_file)"
  else
    if "${cmd[@]}"; then
      create_long_symlinks "${OUTPUT_DIR}" "${TEMPLATE_SUBJECT}" "${TPIDS[@]}"
    else
      rc=$?
      echo "$(date) [ERROR] ${TEMPLATE_SUBJECT} exited with code $rc (foreground)" >> "$ERROR_LOG"
      echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
      exit $rc
    fi
  fi

  exit 0
fi

############################################
# Auto / Re-run Mode: build subject list
############################################
declare -a subjects=()

if [[ $RERUN_MODE -eq 1 ]]; then
  echo "[RE-RUN] Reading subjects from ${RERUN_FILE}"
  while IFS= read -r subject; do
    subjects+=("${BIDS_ROOT%/}/${subject}")
  done < <(jq -r '.subjects[]' "${RERUN_FILE}")
  [[ ${#subjects[@]} -eq 0 ]] && { echo "Error: No subjects found in ${RERUN_FILE}" >&2; exit 1; }
  echo "[RE-RUN] Found ${#subjects[@]} subject(s) to re-run"
else
  shopt -s nullglob
  all_subj=("${BIDS_ROOT%/}"/sub-*)
  [[ ${#all_subj[@]} -eq 0 ]] && { echo "[AUTO] No subjects found in ${BIDS_ROOT}"; exit 1; }
  for sp in "${all_subj[@]}"; do
    [[ -d "$sp" ]] || continue
    sbase=$(basename "$sp")
    mapfile -t ses_list < <(find "$sp" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sort)
    if [[ ${#ses_list[@]} -ge 2 ]]; then
      subjects+=("$sp")
    else
      [[ $DEBUG -eq 1 ]] && echo "[AUTO][SKIP] $sbase has <2 sessions"
    fi
  done
  [[ ${#subjects[@]} -eq 0 ]] && { echo "[AUTO] No longitudinal subjects (>=2 sessions) found"; exit 1; }
  echo "[AUTO] Found ${#subjects[@]} eligible longitudinal subject(s)"
fi

# Apply --pilot
if [[ $PILOT -eq 1 ]]; then
  pick_idx=$(( RANDOM % ${#subjects[@]} ))
  echo "[PILOT] Selected $(basename "${subjects[$pick_idx]}") from ${#subjects[@]} subject(s)"
  subjects=("${subjects[$pick_idx]}")
fi

############################################
# Batch wait helper — waits for current batch, logs errors
############################################
ERROR_LOG="${OUTPUT_DIR%/}/fastsurfer_errors.log"
declare -a batch_pids=()
declare -a batch_subs=()
declare -a batch_logs=()
batch_count=0

wait_for_batch() {
  [[ ${#batch_pids[@]} -eq 0 ]] && return 0
  echo "[BATCH] Waiting for ${#batch_pids[@]} subject(s) to complete..."
  for i in "${!batch_pids[@]}"; do
    pid="${batch_pids[$i]}"
    sbase="${batch_subs[$i]}"
    log="${batch_logs[$i]}"
    if wait "$pid"; then
      echo "  [OK] ${sbase} (PID $pid)"
    else
      rc=$?
      echo "  [ERROR] ${sbase} (PID $pid) exited with code $rc"
      echo "$(date) [ERROR] ${sbase} exited with code $rc PID=$pid" >> "$ERROR_LOG"
      echo "--- tail of ${log} ---" >> "$ERROR_LOG"
      tail -n 300 "$log" >> "$ERROR_LOG" 2>/dev/null
      echo "--- end tail ---" >> "$ERROR_LOG"
    fi
  done
  batch_pids=()
  batch_subs=()
  batch_logs=()
  batch_count=0
  echo "[BATCH] Batch complete."
}

############################################
# Process subjects
############################################
total_processed=0

for sp in "${subjects[@]}"; do
  sbase=$(basename "$sp")
  mapfile -t ses_list < <(find "$sp" -maxdepth 1 -type d -name 'ses-*' -exec basename {} \; | sort)

  TPIDS_LOCAL=()
  T1_PATHS_LOCAL=()
  skip=0
  for ses in "${ses_list[@]}"; do
    anat_dir="$sp/$ses/anat"
    if [[ ! -d "$anat_dir" ]]; then
      echo "[WARN] Missing anat for $sbase/$ses -> skipping subject"; skip=1; break
    fi
    t1=$(ls -1 "$anat_dir"/*_T1w.nii.gz 2>/dev/null | head -n1 || true)
    [[ -z "$t1" ]] && t1=$(ls -1 "$anat_dir"/*_T1w.nii 2>/dev/null | head -n1 || true)
    if [[ -z "$t1" ]]; then
      echo "[WARN] No T1w for $sbase/$ses -> skipping subject"; skip=1; break
    fi
    rel="${t1#${BIDS_ROOT%/}/}"
    TPIDS_LOCAL+=("${sbase}_${ses}")
    T1_PATHS_LOCAL+=("/data/${rel}")
  done
  [[ $skip -eq 1 ]] && continue
  if [[ ${#TPIDS_LOCAL[@]} -lt 2 ]]; then
    [[ $DEBUG -eq 1 ]] && echo "[SKIP] $sbase: insufficient valid sessions"; continue
  fi

  # Skip / force check
  subject_fully_processed=1
  for tpid in "${TPIDS_LOCAL[@]}"; do
    long_dir="${OUTPUT_DIR%/}/${tpid}.long.${sbase}"
    [[ ! -d "$long_dir" || ! -f "${long_dir}/stats/aseg.stats" ]] && {
      subject_fully_processed=0; break; }
  done

  if [[ $subject_fully_processed -eq 1 ]]; then
    if [[ $FORCE -eq 1 ]]; then
      force_delete_long_outputs "${sbase}" "${TPIDS_LOCAL[@]}"
      [[ $DRY_RUN -eq 1 ]] && continue
    else
      echo "[SKIP] $sbase already fully processed (${#TPIDS_LOCAL[@]} sessions)"
      continue
    fi
  fi

  cmd=( singularity exec "${NV_OPTS[@]}" --no-home
    -B "${BIDS_ROOT%/}":/data
    -B "${OUTPUT_DIR%/}":/output
    -B "${LICENSE_DIR}":/fs_license
    "${SIF_FILE}"
    /fastsurfer/long_fastsurfer.sh
    --tid "$sbase"
    --t1s "${T1_PATHS_LOCAL[@]}"
    --tpids "${TPIDS_LOCAL[@]}"
    --sd /output
    --fs_license /fs_license/license.txt
    --py "${PYTHON_CMD}" )
  [[ ${#LONG_OPTS[@]} -gt 0 ]] && cmd+=("${LONG_OPTS[@]}")

  echo "[AUTO] Subject $sbase (${#TPIDS_LOCAL[@]} sessions)"
  printf '  CMD:'; printf ' %q' "${cmd[@]}"; echo

  if [[ $DRY_RUN -eq 1 ]]; then
    [[ $NOHUP -eq 1 ]] && echo "  [DRY RUN] Would run with nohup -> ${OUTPUT_DIR%/}/long_fastsurfer_${sbase}.log"
    [[ $NOHUP -eq 0 ]] && echo "  [DRY RUN] Would run directly"
    continue
  fi

  if [[ $NOHUP -eq 1 ]]; then
    log_file="${OUTPUT_DIR%/}/long_fastsurfer_${sbase}.log"
    echo "  Nohup -> $log_file"
    nohup "${cmd[@]}" > "$log_file" 2>&1 &
    pid=$!
    echo "  PID: $pid"
    batch_pids+=("$pid")
    batch_subs+=("$sbase")
    batch_logs+=("$log_file")
    total_processed=$((total_processed + 1))
    batch_count=$((batch_count + 1))

    # When batch is full, wait before launching more
    if [[ -n "$BATCH_SIZE" && $batch_count -ge $BATCH_SIZE ]]; then
      wait_for_batch
    fi
  else
    if "${cmd[@]}"; then
      total_processed=$((total_processed + 1))
      create_long_symlinks "${OUTPUT_DIR}" "${sbase}" "${TPIDS_LOCAL[@]}"
    else
      rc=$?
      echo "$(date) [ERROR] ${sbase} exited with code $rc (foreground)" >> "$ERROR_LOG"
      echo "CMD: ${cmd[*]}" >> "$ERROR_LOG"
    fi
  fi
done

# Drain any remaining subjects in the final partial batch
if [[ $NOHUP -eq 1 && ${#batch_pids[@]} -gt 0 ]]; then
  wait_for_batch
fi

############################################
# Summary
############################################
if [[ $DRY_RUN -eq 1 ]]; then
  echo "[DONE] Dry run complete."
elif [[ $NOHUP -eq 1 ]]; then
  echo "[DONE] Processed $total_processed subject(s) (batch_size=${BATCH_SIZE})."
  if [[ $total_processed -gt 0 ]]; then
    echo "[INFO] .long symlinks were not created automatically for nohup jobs."
    echo "       Run: bash scripts/create_missing_long_symlinks.sh ${OUTPUT_DIR}"
  fi
else
  echo "[DONE] Processed $total_processed subject(s)."
fi

exit 0
