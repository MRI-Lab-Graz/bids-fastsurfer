#!/usr/bin/env bash
set -euo pipefail

# run_meta_pipeline.sh
# End-to-end, reproducible meta script to:
#   1) Generate QDEC (and aseg.long.table) from participants.tsv and SUBJECTS_DIR
#   2) Run one or more univariate LME models from configs/
#   3) Produce per-effect summaries and an HTML report per model
#
# Defaults align with this repo's configs and folder layout.
#
# Example:
#   bash scripts/run_meta_pipeline.sh \
#     --subjects-dir /path/to/derivatives/fastsurfer \
#     --bids-root /path/to/BIDS \
#     --participants configs/participants.tsv \
#     --out-root results/running_intervention \
#     --models group5_factor_time,group5_decomposed \
#     --link-long
#

# Load shared functions
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/common_functions.sh"

ROOT_DIR="$(get_project_root "$SCRIPT_DIR")"
cd "$ROOT_DIR"

# Get Python command
PY="$(get_python_cmd)"
export PYTHON_CMD="$PY"

PARTICIPANTS="configs/participants.tsv"
SUBJECTS_DIR=""
BIDS_ROOT=""
OUT_ROOT="results/running_intervention"
MODELS="group5_factor_time,group5_decomposed"
VERIFY_MODE="verify"   # "verify" or "link"
DO_ASEG=1               # 1=run asegstats2table via generator, 0=skip
DO_SUMMARIES=1          # 1=run summarization, 0=skip
CONFIG_JSON=""
DRY_RUN=0               # 1=print planned actions only
DO_QC=0                 # 1=run QC outlier check vs eTIV, 0=skip
SKIP_SUBJECTS=""       # comma-separated list of fsid_base to skip
SKIP_FILE=""           # path to file with one fsid_base per line to skip

usage() {
  cat <<EOF
Usage: $0 [options]

Required:
  --subjects-dir PATH         FastSurfer/FreeSurfer SUBJECTS_DIR with sub-*/sub_*_ses-* dirs

Optional:
  --participants PATH         participants.tsv (default: configs/participants.tsv)
  --bids-root PATH            BIDS root to cross-check subjects/sessions (optional)
  --out-root DIR              Results root (default: results/running_intervention)
  --models LIST               Comma-separated models to run (default: group5_factor_time,group5_decomposed)
                              or "all" to run all known models
                              Known keys map to configs:
                                group5_factor_time -> configs/fslmer_group5_factor_time.json
                                group5_decomposed  -> configs/fslmer_group5_decomposed.json
                                group_beh          -> configs/fslmer_group_beh.json
  --link-long                 Create .long symlinks when missing (otherwise just verify)
  --skip-aseg                 Do not run asegstats2table (requires FreeSurfer if omitted)
  --no-summaries              Skip the summarization step (CSV/HTML reports)
  --qc                        Run QC: flag ROI vs eTIV outliers at baseline and emit plots
  --skip-sub LIST             Comma-separated fsid_base IDs to exclude from QDEC and downstream
  --skip-file PATH            File with fsid_base IDs (one per line) to exclude
  --config PATH               JSON config with the same fields; CLI flags override JSON
  --dry-run                   Print planned actions and commands; do not execute
  -h, --help                  Show this help and exit

The script regenerates QDEC and aseg.long.table under --out-root, then runs selected models
and produces summaries/plots per key effect into each model's results directory.

JSON schema (minimal):
{
  "subjects_dir": "/path/to/derivatives/fastsurfer",
  "participants": "configs/participants.tsv",
  "bids_root": "/path/to/BIDS",             # optional
  "out_root": "results/running_intervention",
  "models": ["group5_factor_time", "group5_decomposed"],
  "link_long": true,                          # default: false (verify only)
  "skip_aseg": false,                         # default: false
  "no_summaries": false,                      # default: false
  "qc": false,                                # optional: run QC step
  "skip_sub": ["sub-0001", "sub-0002"],      # optional: subjects (fsid_base) to exclude
  "skip_file": "configs/skip_subjects.txt",   # optional: file with fsid_base per line to exclude
  "effects": {                                # optional: override effect lists per model key
    "group5_factor_time": ["sexM", "age"],
    "group5_decomposed": ["factor(tp)3:weeks4"]
  },
  "model_map": {                              # optional: override config file paths
    "group5_factor_time": "configs/fslmer_group5_factor_time.json"
  }
}

Inline multi-model format (alternative):
{
  "subjects_dir": "/path/to/derivatives/fastsurfer",
  "participants": "configs/participants.tsv",
  "out_root": "results/running_intervention",
  "qdec": "results/running_intervention/qdec.table.dat",
  "aseg": "results/running_intervention/aseg.long.table",
  "time_col": "tp",
  "id_col": "Measure.volume",
  "models": [
    { "id": "first_model",  "formula": "~ tp * group_5 + age + sex", "zcols": "RIRS" },
    { "id": "second_model", "formula": "~ tp * group_5 + age",      "zcols": "RI" }
  ]
}
This will be expanded into per-model configs under <out_root>/.expanded_models/<id>.json with model_id=<id>
and outdir defaulting to <out_root>/fslmer_<id> unless overridden per model.

Human-friendly effect names:
- You can write effects as human-friendly tokens and the pipeline will translate them:
  tp2 -> factor(tp)2
  tp3 -> factor(tp)3
  sex -> sexM
  Interactions: tp3:smallgroup_2w -> factor(tp)3:group_5smallgroup_2w
                tp2:alone_4w      -> factor(tp)2:group_5alone_4w
  Age remains 'age'. Unknown tokens are passed through unchanged.
Structured effects:
- In the meta JSON, effects can also be objects that will be converted to tokens:
  {"tp": 3}                       -> "tp3"
  {"tp": 2, "group_5": "smallgroup_2w"} -> "tp2:smallgroup_2w"
  {"age": true}                  -> "age"
  {"sex": true}                  -> "sex"
Validation of effects:
- Before summarization, each effect is checked against the model's lme_coefficients.csv.
- If an effect isn't present (e.g., tp4 when only tp2/tp3 exist), it's skipped with a clear note.
- For time effects, the script will print which time levels were detected.
EOF
}

# maybe_run is now loaded from common_functions.sh

print_input_summary() {
  echo "[0/4] Input summary" >&2
  echo "subjects_dir: $SUBJECTS_DIR" >&2
  echo "participants: $PARTICIPANTS" >&2
  if [[ -n "$BIDS_ROOT" ]]; then
    echo "bids_root:   $BIDS_ROOT" >&2
  else
    echo "bids_root:   (none)" >&2
  fi
  echo "out_root:    $OUT_ROOT" >&2
  echo "models:      $MODELS" >&2
  echo "long mode:   $VERIFY_MODE" >&2
  echo "aseg:        $([[ $DO_ASEG -eq 1 ]] && echo on || echo off)" >&2
  echo "summaries:   $([[ $DO_SUMMARIES -eq 1 ]] && echo on || echo off)" >&2
  echo "qc:          $([[ $DO_QC -eq 1 ]] && echo on || echo off)" >&2
  if [[ -n "$SKIP_SUBJECTS" || -n "$SKIP_FILE" ]]; then
    echo "skip:       subjects=[$SKIP_SUBJECTS] file=${SKIP_FILE:-none}" >&2
  fi
  if [[ -n "$CONFIG_JSON" ]]; then echo "config:      $CONFIG_JSON" >&2; fi
  if [[ $DRY_RUN -eq 1 ]]; then echo "mode:        DRY-RUN (no changes will be made)" >&2; fi
}

validate_inputs() {
  local ok=1
  validate_path "$SUBJECTS_DIR" dir "subjects_dir" || ok=0
  validate_path "$PARTICIPANTS" file "participants.tsv" || ok=0
  if [[ -n "$BIDS_ROOT" ]]; then
    validate_path "$BIDS_ROOT" dir "bids_root" || ok=0
  fi

  # Validate model configs exist
  local models_to_check="$MODELS"
  if [[ "$models_to_check" == "all" || "$models_to_check" == "ALL" ]]; then
    models_to_check="group5_factor_time,group5_decomposed,group_beh"
  fi
  IFS=',' read -r -a _models <<< "$models_to_check"
  for m in "${_models[@]}"; do
    local cfg=""
    if [[ -f "$m" || ( "$m" == *.json && -f "$m" ) ]]; then
      cfg="$m"
    else
      case "$m" in
        group5_factor_time) cfg="${CFG_GROUP5_FACTOR_TIME:-configs/fslmer_group5_factor_time.json}" ;;
        group5_decomposed)  cfg="${CFG_GROUP5_DECOMP:-configs/fslmer_group5_decomposed.json}"  ;;
        group_beh)          cfg="${CFG_GROUP_BEH:-configs/fslmer_group_beh.json}"          ;;
        *) echo "ERROR: unknown model spec in --models: $m (neither known key nor existing .json file)" >&2; ok=0; continue ;;
      esac
    fi
    [[ -f "$cfg" ]] || { echo "ERROR: model config not found for $m: $cfg" >&2; ok=0; }
  done

  # If aseg is requested, check tool availability
  if [[ $DO_ASEG -eq 1 && $DRY_RUN -eq 0 ]]; then
    if ! check_freesurfer_tools "asegstats2table"; then
      echo "ERROR: asegstats2table not found in PATH. Source FreeSurfer or use --skip-aseg" >&2
      ok=0
    fi
  fi

  if [[ $ok -eq 0 ]]; then
    echo "Validation failed. Fix the above issues and retry." >&2
    exit 2
  fi
}

# json_get and json_get_array are now loaded from common_functions.sh

load_config() {
  local f="$1"
  validate_path "$f" file "Config" || exit 2
  local v
  v=$(json_get "$f" subjects_dir);      [[ -n "$v" ]] && SUBJECTS_DIR="$v"
  v=$(json_get "$f" participants);      [[ -n "$v" ]] && PARTICIPANTS="$v"
  v=$(json_get "$f" bids_root);         [[ -n "$v" ]] && BIDS_ROOT="$v"
  v=$(json_get "$f" out_root);          [[ -n "$v" ]] && OUT_ROOT="$v"
  # models: either array or comma string
  local arr
  arr=$(json_get_array "$f" models || true)
  if [[ -n "$arr" ]]; then
    MODELS=$(echo "$arr" | paste -sd, -)
  else
    v=$(json_get "$f" models)
    [[ -n "$v" ]] && MODELS="$v"
  fi
  v=$(json_get "$f" link_long)
  if [[ "$v" == "true" ]]; then VERIFY_MODE="link"; fi
  v=$(json_get "$f" skip_aseg)
  if [[ "$v" == "true" ]]; then DO_ASEG=0; fi
  v=$(json_get "$f" no_summaries)
  if [[ "$v" == "true" ]]; then DO_SUMMARIES=0; fi
  v=$(json_get "$f" qc)
  if [[ "$v" == "true" ]]; then DO_QC=1; fi
  # Skip subjects: accept array and/or file
  local arr
  arr=$(json_get_array "$f" skip_sub || true)
  if [[ -n "$arr" ]]; then SKIP_SUBJECTS=$(echo "$arr" | paste -sd, -); fi
  v=$(json_get "$f" skip_file); [[ -n "$v" ]] && SKIP_FILE="$v"

  # Optional: override effect lists per model key
  local eff
  eff=$(json_get_array "$f" effects.group5_factor_time || true)
  if [[ -n "$eff" ]]; then mapfile -t EFFECTS_GROUP5_FACTOR_TIME < <(printf "%s\n" "$eff"); fi
  eff=$(json_get_array "$f" effects.group5_decomposed || true)
  if [[ -n "$eff" ]]; then mapfile -t EFFECTS_GROUP5_DECOMP < <(printf "%s\n" "$eff"); fi
  eff=$(json_get_array "$f" effects.group_beh || true)
  if [[ -n "$eff" ]]; then mapfile -t EFFECTS_GROUP_BEH < <(printf "%s\n" "$eff"); fi

  # Optional: override model config paths
  v=$(json_get "$f" model_map.group5_factor_time); [[ -n "$v" ]] && CFG_GROUP5_FACTOR_TIME="$v" || true
  v=$(json_get "$f" model_map.group5_decomposed);  [[ -n "$v" ]] && CFG_GROUP5_DECOMP="$v" || true
  v=$(json_get "$f" model_map.group_beh);          [[ -n "$v" ]] && CFG_GROUP_BEH="$v" || true
}

# Expand inline multi-model meta config to per-model JSONs; updates MODELS to JSON paths
expand_inline_models() {
  local f="$1"
  local out_root="$OUT_ROOT"
  local py_out
  py_out=$("$PY" - "$f" "$OUT_ROOT" <<'PY'
import json,sys,os
fn=sys.argv[1]
out_root=sys.argv[2]
with open(fn) as fh:
  d=json.load(fh)
models=d.get('models')
if not models:
  sys.exit(0)
# Determine if models is list of strings (keys/paths) -> then do nothing
if isinstance(models, list) and all(isinstance(x, str) for x in models):
  sys.exit(0)

# Prepare base fields for fslmer configs
base=dict(d)
for k in ['models','effects','model_map','subjects_dir','participants','bids_root','out_root','link_long','skip_aseg','no_summaries']:
  base.pop(k, None)

expanded_dir=os.path.join(out_root, '.expanded_models')
os.makedirs(expanded_dir, exist_ok=True)

def write_cfg(mid, spec):
  cfg=dict(base)
  cfg.update(spec)
  # Ensure model_id
  cfg['model_id']=spec.get('id', mid)
  # Compute outdir default
  if 'outdir' not in cfg or not cfg['outdir']:
    cfg['outdir']=os.path.join(out_root, f"fslmer_{cfg['model_id']}")
  path=os.path.join(expanded_dir, f"{cfg['model_id']}.json")
  with open(path,'w') as g:
    json.dump(cfg,g,indent=2)
  print(path)

if isinstance(models, dict):
  for mid, spec in models.items():
    if not isinstance(spec, dict):
      continue
    spec=dict(spec)
    spec.setdefault('id', mid)
    write_cfg(mid, spec)
elif isinstance(models, list):
  for i, spec in enumerate(models):
    if not isinstance(spec, dict):
      continue
    mid=spec.get('id') or f"model{i+1}"
    write_cfg(mid, spec)
PY
  )
  if [[ -n "$py_out" ]]; then
  # Convert newline separated to comma list
  local joined
  joined=$(echo "$py_out" | paste -sd, -)
  MODELS="$joined"
  echo "Expanded inline models to: $MODELS" >&2
  fi
}

# human_to_effect, effect_exists, list_time_effects are now loaded from common_functions.sh

# Read effects for a given tag from meta CONFIG_JSON, supporting strings and simple objects
get_effects_for_tag() {
  local config="$1" tag="$2"
  "$PY" - "$config" "$tag" <<'PY'
import json,sys
fn,tag=sys.argv[1],sys.argv[2]
try:
  with open(fn) as f:
    d=json.load(f)
except Exception:
  d={}
effs=(d.get('effects') or {}).get(tag)
if not effs:
  sys.exit(0)

def emit_token(x):
  if isinstance(x,str):
    print(x)
    return
  if not isinstance(x,dict):
    return
  tp=None
  grp=None
  # time
  if isinstance(x.get('tp'), int):
    tp=f"tp{x['tp']}"
  elif isinstance(x.get('time'), int):
    tp=f"tp{x['time']}"
  # group 5
  g5=x.get('group_5')
  if isinstance(g5,str):
    grp=g5
  else:
    for g in ['alone_4w','smallgroup_2w','smallgroup_4w','control']:
      if x.get(g) is True:
        grp=g
        break
  # simple boolean main effects
  for simple in ['sex','age','weeks4','smallgroup']:
    if x.get(simple) is True:
      print(simple)
  # combined tokens
  if tp and grp:
    print(f"{tp}:{grp}")
  elif tp:
    print(tp)
  elif grp:
    print(grp)

for item in effs:
  emit_token(item)
PY
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)       CONFIG_JSON="$2"; load_config "$2"; shift 2 ;;
    --subjects-dir) SUBJECTS_DIR="$2"; shift 2 ;;
    --participants) PARTICIPANTS="$2"; shift 2 ;;
    --bids-root)    BIDS_ROOT="$2"; shift 2 ;;
    --out-root)     OUT_ROOT="$2"; shift 2 ;;
    --models)       MODELS="$2"; shift 2 ;;
    --link-long)    VERIFY_MODE="link"; shift ;;
    --skip-aseg)    DO_ASEG=0; shift ;;
    --no-summaries) DO_SUMMARIES=0; shift ;;
  --qc)          DO_QC=1; shift ;;
  --skip-sub)    SKIP_SUBJECTS="$2"; shift 2 ;;
  --skip-file)   SKIP_FILE="$2"; shift 2 ;;
    --dry-run)      DRY_RUN=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$SUBJECTS_DIR" ]]; then
  ENV_SUBJECTS_DIR="$(printenv SUBJECTS_DIR 2>/dev/null || true)"
  if [[ -n "$ENV_SUBJECTS_DIR" ]]; then
    SUBJECTS_DIR="$ENV_SUBJECTS_DIR"
    echo "Using SUBJECTS_DIR from environment: $SUBJECTS_DIR" >&2
  fi
fi
[[ -n "$SUBJECTS_DIR" ]] || { echo "--subjects-dir is required (or export SUBJECTS_DIR)" >&2; exit 2; }
[[ -f "$PARTICIPANTS" ]] || { echo "participants.tsv not found: $PARTICIPANTS" >&2; exit 2; }

print_input_summary
# If a meta config with inline models was provided, expand it now
if [[ -n "$CONFIG_JSON" ]]; then
  expand_inline_models "$CONFIG_JSON"
fi
validate_inputs

if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY-RUN: would create directory $OUT_ROOT" >&2
else
  mkdir -p "$OUT_ROOT"
fi

QDEC_PATH="$OUT_ROOT/qdec.table.dat"
ASEG_PATH="$OUT_ROOT/aseg.long.table"   # generator writes this next to QDEC when --aseg is used

echo "[1/4] Generating QDEC (and aseg table) into $OUT_ROOT" >&2
[[ -n "$BIDS_ROOT" ]] && BIDS_ARG=(--bids "$BIDS_ROOT") || BIDS_ARG=()

LINK_ARGS=()
if [[ "$VERIFY_MODE" == "link" ]]; then
  LINK_ARGS+=(--link-long)
else
  LINK_ARGS+=(--verify-long)
fi

ASEG_FLAG=()
if [[ $DO_ASEG -eq 1 ]]; then
  ASEG_FLAG+=(--aseg)
fi

SKIP_ARGS=()
if [[ -n "$SKIP_SUBJECTS" ]]; then
  SKIP_ARGS+=(--skip-sub "$SKIP_SUBJECTS")
fi
if [[ -n "$SKIP_FILE" ]]; then
  SKIP_ARGS+=(--skip-file "$SKIP_FILE")
fi

maybe_run "$PY" scripts/generate_qdec.py \
  --participants "$PARTICIPANTS" \
  --subjects-dir "$SUBJECTS_DIR" \
  --output "$QDEC_PATH" \
  "${BIDS_ARG[@]}" \
  "${LINK_ARGS[@]}" \
  "${ASEG_FLAG[@]}" \
  "${SKIP_ARGS[@]}"

echo "QDEC: $QDEC_PATH" >&2
if [[ $DRY_RUN -eq 0 ]]; then
  [[ -f "$QDEC_PATH" ]] || { echo "QDEC not found after generation: $QDEC_PATH" >&2; exit 3; }
  if [[ $DO_ASEG -eq 1 ]]; then
    [[ -f "$ASEG_PATH" ]] || { echo "aseg.long.table not found: $ASEG_PATH (source FreeSurfer or use --skip-aseg)" >&2; exit 3; }
  fi
else
  echo "DRY-RUN: would validate existence of $QDEC_PATH and $ASEG_PATH (if requested)" >&2
fi

# Optional QC step
if [[ $DO_QC -eq 1 ]]; then
  echo "[1b/4] QC: ROI vs eTIV outlier check (baseline)" >&2
  if [[ $DRY_RUN -eq 1 ]]; then
    echo "DRY-RUN: Rscript scripts/qc_etiv_outliers.R --qdec $QDEC_PATH --aseg $ASEG_PATH --outdir $OUT_ROOT/qc --time-col tp --id-col Measure.volume" >&2
  else
    maybe_run Rscript scripts/qc_etiv_outliers.R --qdec "$QDEC_PATH" --aseg "$ASEG_PATH" --outdir "$OUT_ROOT/qc" --time-col tp --id-col Measure.volume || true
  fi
fi

if [[ "$MODELS" == "all" || "$MODELS" == "ALL" ]]; then
  MODELS="group5_factor_time,group5_decomposed,group_beh"
fi

echo "[2/4] Running selected models: $MODELS" >&2
IFS=',' read -r -a MODELS_ARR <<< "$MODELS"

# Keep track for summarization: arrays of result dirs and tags
declare -a __RES_DIRS=()
declare -a __TAGS=()

run_model() {
  local spec="$1" cfg="" tag=""
  # Resolve spec to config path and a tag (key or model_id)
  if [[ -f "$spec" || ( "$spec" == *.json && -f "$spec" ) ]]; then
    cfg="$spec"
    tag="$(json_get "$cfg" model_id)"
    if [[ -z "$tag" ]]; then
      tag="$(basename "$cfg" .json)"
    fi
  else
    # Allow override via JSON model_map
    local cfg_override=""
    case "$spec" in
      group5_factor_time) cfg_override="${CFG_GROUP5_FACTOR_TIME:-}" ;;
      group5_decomposed)  cfg_override="${CFG_GROUP5_DECOMP:-}" ;;
      group_beh)          cfg_override="${CFG_GROUP_BEH:-}" ;;
    esac
    case "$spec" in
      group5_factor_time) cfg="${cfg_override:-configs/fslmer_group5_factor_time.json}" ;;
      group5_decomposed)  cfg="${cfg_override:-configs/fslmer_group5_decomposed.json}"  ;;
      group_beh)          cfg="${cfg_override:-configs/fslmer_group_beh.json}"          ;;
      *) echo "Unknown model key: $spec" >&2; return 2 ;;
    esac
    tag="$spec"
  fi

  echo "-- $tag -> $cfg" >&2
  # Ensure the config uses the OUT_ROOT paths (our configs already do). If you customize OUT_ROOT,
  # make sure config JSONs point to relative paths under $OUT_ROOT.
  # Pass through region selection if present; default to --all-regions for multi-region analyses
  local rp
  rp="$(json_get "$cfg" region_pattern)"; if [[ -z "$rp" ]]; then rp="$(json_get "$cfg" region-pattern)"; fi
  local ar
  ar="$(json_get "$cfg" all_regions)"; if [[ -z "$ar" ]]; then ar="$(json_get "$cfg" all-regions)"; fi
  if [[ -n "$rp" ]]; then
    maybe_run bash scripts/run_fslmer_univariate.sh --config "$cfg" --region-pattern "$rp"
  else
    if [[ "$ar" == "true" || -z "$ar" ]]; then
      maybe_run bash scripts/run_fslmer_univariate.sh --config "$cfg" --all-regions
    else
      maybe_run bash scripts/run_fslmer_univariate.sh --config "$cfg"
    fi
  fi

  # Collect result dir for summarization
  local res_dir
  res_dir="$(json_get "$cfg" outdir)"
  if [[ -z "$res_dir" ]]; then
    # Fallback to OUT_ROOT/tag if not specified (legacy)
    res_dir="$OUT_ROOT/fslmer_${tag}"
  fi
  __RES_DIRS+=("$res_dir")
  __TAGS+=("$tag")
}

for m in "${MODELS_ARR[@]}"; do
  run_model "$m"
done

if [[ $DO_SUMMARIES -eq 1 ]]; then
  echo "[3/4] Summarizing key effects per model" >&2

summarize_effects() {
  local res_dir="$1"; shift
  local effects=("$@")
  for eff in "${effects[@]}"; do
    local eff_resolved
    eff_resolved="$(human_to_effect "$eff")"
    if effect_exists "$res_dir" "$eff_resolved"; then
      maybe_run Rscript scripts/fslmer_summarize.R --results-dir "$res_dir" --effect "$eff_resolved" --verbose || true
    else
      echo "Note: Skipping effect '$eff' -> '$eff_resolved' (not found in coefficients)." >&2
      if [[ "$eff_resolved" =~ ^factor\(tp\)[0-9]+(|:.*)$ ]]; then
        echo -n "       Available time effects: " >&2
        list_time_effects "$res_dir" >&2 || true
      fi
    fi
  done
}

# Effects to summarize for each known model
EFFECTS_GROUP5_FACTOR_TIME=(
  # Baseline group differences vs reference (alone_2w)
  "group_5alone_4w" "group_5smallgroup_2w" "group_5smallgroup_4w" "group_5control"
  # Time main effects (change from baseline in reference arm)
  "factor(tp)2" "factor(tp)3"
  # Interactions (differential change vs reference per group and time)
  "factor(tp)2:group_5alone_4w" "factor(tp)3:group_5alone_4w"
  "factor(tp)2:group_5smallgroup_2w" "factor(tp)3:group_5smallgroup_2w"
  "factor(tp)2:group_5smallgroup_4w" "factor(tp)3:group_5smallgroup_4w"
  "factor(tp)2:group_5control" "factor(tp)3:group_5control"
  # Covariates of interest
  "sexM" "age"
)

EFFECTS_GROUP5_DECOMP=(
  # Time and interactions for decomposed design
  "factor(tp)2" "factor(tp)3"
  "smallgroup" "weeks4" "sexM" "age"
  "factor(tp)2:smallgroup" "factor(tp)3:smallgroup"
  "factor(tp)2:weeks4" "factor(tp)3:weeks4"
)

EFFECTS_GROUP_BEH=(
  # Example; adjust to your factors used in configs/fslmer_group_beh.json
  "factor(tp)2" "factor(tp)3" "sexM" "age"
)

  for i in "${!__RES_DIRS[@]}"; do
    res_dir="${__RES_DIRS[$i]}"
    tag="${__TAGS[$i]}"

    # Try per-config overrides first if a meta CONFIG_JSON was provided
    effects_override=()
    if [[ -n "$CONFIG_JSON" ]]; then
      mapfile -t effects_override < <(get_effects_for_tag "$CONFIG_JSON" "$tag" || true)
    fi

    if [[ ${#effects_override[@]} -gt 0 ]]; then
      summarize_effects "$res_dir" "${effects_override[@]}"
      continue
    fi

    case "$tag" in
      group5_factor_time)
        summarize_effects "$res_dir" "${EFFECTS_GROUP5_FACTOR_TIME[@]}"
        ;;
      group5_decomposed)
        summarize_effects "$res_dir" "${EFFECTS_GROUP5_DECOMP[@]}"
        ;;
      group_beh)
        summarize_effects "$res_dir" "${EFFECTS_GROUP_BEH[@]}"
        ;;
      *)
        echo "Note: No built-in effect list for '$tag' and no overrides in meta JSON. Add effects overrides under 'effects' to summarize this model." >&2
        ;;
    esac
  done
else
  echo "[3/4] Skipping summarization (per --no-summaries)" >&2
fi

echo "[4/4] Done. Artifacts under: $OUT_ROOT" >&2
echo "- QDEC: $QDEC_PATH" >&2
if [[ -f "$ASEG_PATH" ]]; then echo "- aseg.long.table: $ASEG_PATH" >&2; fi
for m in "${MODELS_ARR[@]}"; do
  case "$m" in
    group5_factor_time) echo "- Model outputs: $OUT_ROOT/fslmer_group5_factor_time" >&2 ;;
    group5_decomposed)  echo "- Model outputs: $OUT_ROOT/fslmer_group5_decomposed" >&2 ;;
    group_beh)          echo "- Model outputs: $OUT_ROOT/fslmer_group_beh" >&2 ;;
  esac
done
