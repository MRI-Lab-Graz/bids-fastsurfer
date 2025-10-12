#!/usr/bin/env bash
# common_functions.sh
# Shared functions and utilities for FastSurfer longitudinal analysis pipeline scripts
# Source this file in other scripts: source "$(dirname "$0")/common_functions.sh"

# Error handling and display functions
show_error() {
    echo "ERROR: $1" >&2
    echo "Use --help for usage examples and argument details." >&2
    exit 1
}

show_warning() {
    echo "WARNING: $1" >&2
}

# Detect Python command (python or python3)
get_python_cmd() {
    local py=python
    command -v "$py" >/dev/null 2>&1 || py=python3
    if ! command -v "$py" >/dev/null 2>&1; then
        show_error "Python not found. Install Python 3 to continue."
    fi
    echo "$py"
}

# JSON parsing functions
json_get() {
    local file="$1" key="$2"
    local py="${PYTHON_CMD:-$(get_python_cmd)}"
    "$py" - "$file" "$key" <<'PY'
import json,sys
fn,key=sys.argv[1],sys.argv[2]
with open(fn) as f: d=json.load(f)
cur=d
for part in key.split('.'):
    if isinstance(cur, dict) and part in cur:
        cur=cur[part]
    else:
        cur=None; break
if cur is None:
    print("")
elif isinstance(cur, bool):
    print("true" if cur else "false")
else:
    print(cur)
PY
}

json_get_array() {
    local file="$1" key="$2"
    local py="${PYTHON_CMD:-$(get_python_cmd)}"
    "$py" - "$file" "$key" <<'PY'
import json,sys
fn,key=sys.argv[1],sys.argv[2]
with open(fn) as f: d=json.load(f)
cur=d
for part in key.split('.'):
    if isinstance(cur, dict) and part in cur:
        cur=cur[part]
    else:
        cur=None; break
if isinstance(cur, list):
    for v in cur:
        print(v)
PY
}

# Environment setup and R detection
setup_r_environment() {
    # Try micromamba environment first
    local script_dir="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
    local env_record="$script_dir/.mamba_env"
    local env_name="fastsurfer-r"
    local mamba_root_prefix="${MAMBA_ROOT_PREFIX:-$HOME/.local/micromamba}"
    local mamba_bin="$mamba_root_prefix/bin/micromamba"
    
    # Load environment record if it exists
    if [[ -f "$env_record" ]]; then
        # shellcheck disable=SC1090
        source "$env_record"
        env_name="${ENV_NAME:-fastsurfer-r}"
        mamba_root_prefix="${MAMBA_ROOT_PREFIX:-$HOME/.local/micromamba}"
        mamba_bin="$mamba_root_prefix/bin/micromamba"
    fi
    
    # Check for required R packages function
    check_r_packages() {
        local rscript_cmd="$1"
        echo "Checking R package dependencies..." >&2
        
        if ! $rscript_cmd -e "
            required_pkgs <- c('optparse', 'jsonlite', 'stats')
            missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]
            if (length(missing) > 0) {
                cat('Missing required R packages:', paste(missing, collapse=', '), '\n', file=stderr())
                cat('Run scripts/install.sh to install required packages.\n', file=stderr())
                quit(status=1)
            }
            cat('All required R packages available.\n', file=stderr())
        " 2>/dev/null; then
            echo "ERROR: Missing required R packages. Run scripts/install.sh to install dependencies." >&2
            return 1
        fi
    }
    
    # Try micromamba environment first
    if [[ -x "$mamba_bin" && -d "$mamba_root_prefix/envs/$env_name" ]]; then
        echo "Using micromamba environment: $env_name" >&2
        local rscript_cmd="$mamba_bin run -p $mamba_root_prefix/envs/$env_name Rscript"
        
        # Verify environment works
        if check_r_packages "$rscript_cmd"; then
            echo "$rscript_cmd"
            return 0
        else
            echo "WARNING: micromamba environment has missing packages, falling back to system R" >&2
        fi
    fi
    
    # Fallback: system Rscript
    if ! command -v Rscript >/dev/null 2>&1; then
        show_error "Rscript not found in PATH and micromamba environment unavailable. Install R or run scripts/install.sh to set up environment."
    fi
    
    echo "Using system R installation" >&2
    local rscript_cmd="Rscript"
    
    # Check system R packages
    if ! check_r_packages "$rscript_cmd"; then
        show_error "System R missing required packages. Run scripts/install.sh or install packages manually: install.packages(c('optparse', 'jsonlite'))"
    fi
    
    echo "$rscript_cmd"
}

# Conditional execution for dry-run mode
maybe_run() {
    if [[ "${DRY_RUN:-0}" -eq 1 ]]; then
        # Print the command as it would be executed
        printf 'DRY-RUN: '
        printf '%q ' "$@"
        printf '\n'
        return 0
    else
        "$@"
    fi
}

# Validate file/directory existence
validate_path() {
    local path="$1"
    local type="${2:-file}"  # file or dir
    local name="${3:-Path}"
    
    case "$type" in
        file)
            [[ -f "$path" ]] || { echo "ERROR: $name not found: $path" >&2; return 1; }
            ;;
        dir)
            [[ -d "$path" ]] || { echo "ERROR: $name directory not found: $path" >&2; return 1; }
            ;;
        *)
            echo "ERROR: Invalid validation type: $type" >&2
            return 1
            ;;
    esac
    return 0
}

# Check if FreeSurfer tools are available
check_freesurfer_tools() {
    local tool="${1:-asegstats2table}"
    if command -v "$tool" >/dev/null 2>&1; then
        return 0
    elif [[ -x "/usr/local/freesurfer/bin/$tool" ]]; then
        return 0
    else
        return 1
    fi
}

# Get absolute path to project root (assuming scripts/ subdirectory)
get_project_root() {
    local script_dir="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
    cd "$script_dir/.." && pwd
}

# Parse comma-separated list into array
parse_csv_to_array() {
    local csv="$1"
    local -n arr_ref="$2"
    IFS=',' read -r -a arr_ref <<< "$csv"
}

# Convert human-friendly effect tokens to internal effect terms
human_to_effect() {
    local e="$1"
    # Lowercase copy for matching but preserve original for other cases
    local el="${e,,}"
    # Interaction?
    if [[ "$el" == *:* ]]; then
        local lhs="${el%%:*}" rhs="${el#*:}"
        lhs="$(human_to_effect "$lhs")"
        # Map group tokens on RHS
        case "$rhs" in
            smallgroup_2w|smallgroup-2w) rhs="group_5smallgroup_2w" ;;
            smallgroup_4w|smallgroup-4w) rhs="group_5smallgroup_4w" ;;
            alone_4w|alone-4w)           rhs="group_5alone_4w" ;;
            control)                     rhs="group_5control" ;;
            *)                           rhs="$(human_to_effect "$rhs")" ;;
        esac
        echo "${lhs}:${rhs}"
        return 0
    fi
    case "$el" in
        tp2|time2|t2) echo "factor(tp)2" ;;
        tp3|time3|t3) echo "factor(tp)3" ;;
        sex|gender)   echo "sexM" ;;
        age)          echo "age" ;;
        *)            echo "$e" ;;
    esac
}

# Check if an effect term exists in lme_coefficients.csv
effect_exists() {
    local res_dir="$1" eff="$2"
    local py="${PYTHON_CMD:-$(get_python_cmd)}"
    "$py" - "$res_dir" "$eff" <<'PY'
import sys,os,csv
res_dir,eff=sys.argv[1],sys.argv[2]
fp=os.path.join(res_dir,'lme_coefficients.csv')
if not os.path.isfile(fp):
    sys.exit(2)
with open(fp,newline='') as f:
    r=csv.DictReader(f)
    for row in r:
        if row.get('coef')==eff:
            sys.exit(0)
sys.exit(1)
PY
    case $? in
        0) return 0;;
        1) return 1;;
        2) echo "WARN: coefficients file not found: $res_dir/lme_coefficients.csv" >&2; return 1;;
    esac
}

# List available time effect levels from lme_coefficients.csv
list_time_effects() {
    local res_dir="$1"
    local py="${PYTHON_CMD:-$(get_python_cmd)}"
    "$py" - "$res_dir" <<'PY'
import sys,os,csv,re
res_dir=sys.argv[1]
fp=os.path.join(res_dir,'lme_coefficients.csv')
levels=set()
if os.path.isfile(fp):
    with open(fp,newline='') as f:
        for row in csv.DictReader(f):
            c=row.get('coef','')
            m=re.match(r'^factor\(tp\)(\d+)', c)
            if m:
                levels.add(int(m.group(1)))
if levels:
    print(', '.join(f'tp{n}' for n in sorted(levels)))
else:
    print('(none)')
PY
}

# Export functions for use in scripts that source this file
export -f show_error show_warning
export -f json_get json_get_array
export -f maybe_run validate_path
export -f human_to_effect effect_exists list_time_effects
