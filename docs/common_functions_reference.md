# Common Functions Quick Reference

This is a quick reference guide for developers using `scripts/common_functions.sh`.

## How to Use

Add this to the top of your bash script:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Load shared functions
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1090
source "$SCRIPT_DIR/common_functions.sh"
```

## Available Functions

### Error Handling

```bash
# Display error and exit with code 1
show_error "Something went wrong"

# Display warning (doesn't exit)
show_warning "This might be a problem"
```

### Path Utilities

```bash
# Get absolute path to project root
ROOT_DIR="$(get_project_root "$SCRIPT_DIR")"
cd "$ROOT_DIR"

# Validate file exists
validate_path "/path/to/file.txt" file "Config file" || exit 1

# Validate directory exists
validate_path "/path/to/dir" dir "Data directory" || exit 1
```

### Python Detection

```bash
# Get python or python3 command
PY="$(get_python_cmd)"
export PYTHON_CMD="$PY"

# Use it
$PY script.py --arg value
```

### JSON Configuration

```bash
# Get single value from JSON
value=$(json_get "config.json" "key")
nested=$(json_get "config.json" "parent.child.key")
boolean=$(json_get "config.json" "enabled")  # returns "true" or "false"

# Get array from JSON (prints one value per line)
json_get_array "config.json" "models" | while read -r model; do
  echo "Model: $model"
done

# Or into array
mapfile -t models < <(json_get_array "config.json" "models")
```

### R Environment Setup

```bash
# Detect R environment (tries micromamba, falls back to system R)
# Validates required packages are installed
RSCRIPT_CMD="$(setup_r_environment "$SCRIPT_DIR")"

# Use it
$RSCRIPT_CMD my_analysis.R --data input.csv
```

### FreeSurfer Tools

```bash
# Check if FreeSurfer tool is available
if check_freesurfer_tools "asegstats2table"; then
  echo "FreeSurfer is available"
else
  echo "FreeSurfer not found"
fi
```

### Dry Run Support

```bash
# Set dry run mode
DRY_RUN=1

# Use maybe_run for conditional execution
maybe_run echo "This will only print in dry-run mode"
maybe_run rm -rf /important/data  # Safe! Only prints in dry-run

# Normal execution
DRY_RUN=0
maybe_run echo "This actually runs"
```

### Statistical Analysis Helpers

```bash
# Convert human-friendly effect names to R terms
effect=$(human_to_effect "tp2")              # "factor(tp)2"
effect=$(human_to_effect "sex")              # "sexM"
effect=$(human_to_effect "tp3:smallgroup_2w") # "factor(tp)3:group_5smallgroup_2w"

# Check if effect exists in model
if effect_exists "/path/to/results" "factor(tp)2"; then
  echo "Effect found in model"
fi

# List available time effects
time_effects=$(list_time_effects "/path/to/results")
echo "Available: $time_effects"  # e.g., "tp2, tp3"
```

## Real-World Examples

### Example 1: Simple Analysis Script

```bash
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/common_functions.sh"

# Setup
ROOT_DIR="$(get_project_root "$SCRIPT_DIR")"
PY="$(get_python_cmd)"

# Validate inputs
CONFIG="$1"
validate_path "$CONFIG" file "Configuration" || exit 1

# Parse config
DATA_DIR=$(json_get "$CONFIG" "data_dir")
validate_path "$DATA_DIR" dir "Data directory" || exit 1

# Run analysis
echo "Running analysis..."
$PY analyze.py --config "$CONFIG"
```

### Example 2: R Analysis with Dry Run

```bash
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/common_functions.sh"

# Parse arguments
DRY_RUN=0
INPUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    *) show_error "Unknown argument: $1" ;;
  esac
done

# Validate
[[ -n "$INPUT" ]] || show_error "Missing --input"
validate_path "$INPUT" file "Input file" || exit 1

# Setup R
RSCRIPT_CMD="$(setup_r_environment "$SCRIPT_DIR")"

# Run (respects dry-run)
maybe_run $RSCRIPT_CMD analysis.R --input "$INPUT"
```

### Example 3: Multi-Model Pipeline

```bash
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/common_functions.sh"

CONFIG="configs/meta.json"
PY="$(get_python_cmd)"

# Get models from config
mapfile -t MODELS < <(json_get_array "$CONFIG" "models")

# Run each model
for model in "${MODELS[@]}"; do
  echo "Processing model: $model"
  
  # Get model-specific settings
  formula=$(json_get "configs/${model}.json" "formula")
  outdir=$(json_get "configs/${model}.json" "outdir")
  
  # Run analysis
  bash scripts/run_analysis.sh \
    --model "$model" \
    --formula "$formula" \
    --outdir "$outdir"
done
```

## Environment Variables

These are set/used by common functions:

- `PYTHON_CMD` - Python executable (set by `get_python_cmd`)
- `DRY_RUN` - Set to `1` for dry-run mode (used by `maybe_run`)
- `MAMBA_ROOT_PREFIX` - Micromamba installation path (used by `setup_r_environment`)

## Tips and Best Practices

### 1. Always validate inputs early
```bash
validate_path "$INPUT_FILE" file "Input" || exit 1
validate_path "$OUTPUT_DIR" dir "Output" || exit 1
# Now you can safely use these paths
```

### 2. Use meaningful names in validation
```bash
# Bad - unclear what failed
validate_path "$f" file || exit 1

# Good - clear error message
validate_path "$config_file" file "Configuration file" || exit 1
```

### 3. Export PYTHON_CMD for consistency
```bash
PY="$(get_python_cmd)"
export PYTHON_CMD="$PY"  # Makes it available to json_get/etc
```

### 4. Check dry-run mode for destructive operations
```bash
# Good - respects dry-run
maybe_run rm -rf "$output_dir"

# Bad - always runs
rm -rf "$output_dir"
```

### 5. Handle missing JSON values
```bash
# JSON returns empty string if key missing
value=$(json_get "config.json" "optional_key")
if [[ -z "$value" ]]; then
  value="default_value"
fi
```

## Troubleshooting

### Error: "Config not found"
- Check file path is correct
- Use absolute paths or make sure you're in the right directory

### Error: "Missing required R packages"
- Run `scripts/install.sh` to set up environment
- Or install manually: `R -e "install.packages(c('optparse', 'jsonlite'))"`

### Error: "Python not found"
- Install Python 3
- Ensure it's in your PATH

### Functions not found
- Make sure you sourced the file: `source "$SCRIPT_DIR/common_functions.sh"`
- Check the path to common_functions.sh is correct

## Testing Your Script

Always test with these scenarios:

```bash
# 1. Normal execution
bash your_script.sh --input data.csv

# 2. Dry-run mode
bash your_script.sh --input data.csv --dry-run

# 3. Missing required files
bash your_script.sh --input nonexistent.csv  # Should show clear error

# 4. Invalid arguments
bash your_script.sh --invalid-flag  # Should show error

# 5. Help flag
bash your_script.sh --help  # Should show usage
```

## See Also

- `REFACTORING_SUMMARY.md` - Overview of the refactoring
- `run_fslmer_univariate.sh` - Example usage in simple script
- `run_meta_pipeline.sh` - Example usage in complex pipeline
