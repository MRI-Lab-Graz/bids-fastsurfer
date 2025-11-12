#!/usr/bin/env bash
# Wrapper to run analyse_qdec.py with nohup so you can detach from a remote session.
# Usage: scripts/run_analyse_qdec_nohup.sh [args-for-analyse_qdec.py]
# Example:
# ./scripts/run_analyse_qdec_nohup.sh --qdec /tmp/qdec.table.dat --subjects-dir /data/local/129_PK01/derivatives/fastsurfer/ --output results/qdec_analysis --aseg --aparc --surf --surf-measures thickness --force

set -euo pipefail

# Path to the micromamba/python environment used for this project
PYTHON=~/.local/micromamba/envs/fastsurfer-r/bin/python
ANALYSE_SCRIPT="$(pwd)/scripts/analyse_qdec.py"

if [ ! -x "$PYTHON" ] && [ ! -f "$PYTHON" ]; then
  echo "Warning: python not found at $PYTHON. Using 'python' from PATH."
  PYTHON="python"
fi

# Default output dir if not provided
OUT_DIR="results/qdec_analysis"

# Parse args quickly to find --output value (if given)
ARGS=("$@")
for ((i=0;i<${#ARGS[@]};i++)); do
  arg="${ARGS[$i]}"
  case "$arg" in
    --output)
      if [ $((i+1)) -lt ${#ARGS[@]} ]; then
        OUT_DIR="${ARGS[$((i+1))]}"
      fi
      ;;
    --output=*)
      OUT_DIR="${arg#--output=}"
      ;;
  esac
done

# Ensure output dir exists so we can place a log there
mkdir -p "$OUT_DIR"

# Create a timestamped log file inside output dir
TS=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$OUT_DIR/analyse_qdec.nohup.$TS.log"

# Build the command to run
CMD=("$PYTHON" "$ANALYSE_SCRIPT" "$@")

echo "Launching analyse_qdec.py with nohup"
echo "Log: $LOG_FILE"

echo "nohup ${CMD[*]} >'$LOG_FILE' 2>&1 &"
nohup "${CMD[@]}" >"$LOG_FILE" 2>&1 &
PID=$!

echo "Started analyse_qdec.py (PID=$PID)."
echo "To follow logs: tail -f $LOG_FILE"
echo "To disown: disown $PID (if using a shell that supports it)"

# Print PID and log path in a machine-parsable form too
printf '{"pid": %d, "log": "%s"}\n' "$PID" "$LOG_FILE"

exit 0
