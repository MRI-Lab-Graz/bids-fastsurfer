# Script Refactoring Summary

## Overview

The `run_fslmer_univariate.sh` and `run_meta_pipeline.sh` scripts have been refactored to eliminate code duplication by extracting common functionality into a shared library.

## Changes Made

### New File Created

**`scripts/common_functions.sh`** - Shared library containing:

#### Core Utilities
- `show_error()` - Standardized error reporting
- `show_warning()` - Standardized warning messages
- `get_python_cmd()` - Detect Python 2 or 3 executable
- `get_project_root()` - Get absolute path to project root
- `maybe_run()` - Conditional execution for dry-run mode

#### Path Validation
- `validate_path()` - Generic file/directory existence checker
- `check_freesurfer_tools()` - Check for FreeSurfer commands

#### JSON Parsing
- `json_get()` - Extract single values from JSON configs
- `json_get_array()` - Extract arrays from JSON configs

#### R Environment Setup
- `setup_r_environment()` - Detect and configure R environment
  - Tries micromamba environment first
  - Falls back to system R if needed
  - Validates required R packages

#### Statistical Analysis Helpers
- `human_to_effect()` - Convert human-readable effect names to R terms
  - `tp2` → `factor(tp)2`
  - `sex` → `sexM`
  - `tp3:smallgroup_2w` → `factor(tp)3:group_5smallgroup_2w`
- `effect_exists()` - Check if effect exists in model coefficients
- `list_time_effects()` - List available time effect levels

### Modified Files

#### `scripts/run_fslmer_univariate.sh`

**Before:** ~160 lines with embedded R environment setup logic

**After:** ~40 lines, delegates to common functions

**Changes:**
- Sources `common_functions.sh` at the beginning
- Removed duplicate environment setup code (~100 lines)
- Uses `setup_r_environment()` for R detection
- Uses `validate_path()` for file validation
- Uses `get_project_root()` for directory navigation

#### `scripts/run_meta_pipeline.sh`

**Before:** ~750 lines with many duplicate functions

**After:** ~640 lines, cleaner separation of concerns

**Changes:**
- Sources `common_functions.sh` at the beginning
- Removed duplicate `json_get()` and `json_get_array()` (~40 lines)
- Removed duplicate `maybe_run()` (~10 lines)
- Removed duplicate `human_to_effect()` (~25 lines)
- Removed duplicate `effect_exists()` (~20 lines)
- Removed duplicate `list_time_effects()` (~20 lines)
- Uses `get_python_cmd()` for Python detection
- Uses `validate_path()` for path validation
- Uses `check_freesurfer_tools()` for FreeSurfer tool detection

## Benefits

### 1. **Reduced Code Duplication**
- ~200 lines of duplicate code eliminated
- Single source of truth for common functionality
- Easier to maintain and update

### 2. **Improved Consistency**
- Both scripts now use identical error handling
- Standardized JSON parsing behavior
- Consistent environment detection logic

### 3. **Better Testability**
- Common functions can be tested independently
- Changes to shared logic automatically propagate to both scripts

### 4. **Enhanced Maintainability**
- Bug fixes in common code benefit both scripts
- New features can be added to shared library
- Clearer separation between script-specific and shared logic

### 5. **Easier Extension**
- New scripts can easily reuse common functionality
- Simply source `common_functions.sh` and access all utilities

## Backward Compatibility

✅ **Fully backward compatible** - No changes to:
- Command-line interfaces
- Configuration file formats
- Input/output behavior
- Environment variables
- Error messages or exit codes

## Usage

### For New Scripts

To use the shared functions in a new script:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Load shared functions
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "$SCRIPT_DIR/common_functions.sh"

# Now you can use all common functions
ROOT_DIR="$(get_project_root "$SCRIPT_DIR")"
PY="$(get_python_cmd)"
RSCRIPT_CMD="$(setup_r_environment "$SCRIPT_DIR")"

# Use validation
validate_path "/some/file.txt" file "My file" || exit 1

# Use JSON parsing
value=$(json_get "config.json" "some.nested.key")
```

### For Existing Scripts

No changes needed! Both refactored scripts work exactly as before:

```bash
# Same usage as always
bash scripts/run_fslmer_univariate.sh --help
bash scripts/run_meta_pipeline.sh --subjects-dir /path/to/data
```

## Testing

Both scripts have been tested and verified to work:
- ✅ Help flags display correctly
- ✅ Environment setup works (micromamba and system R)
- ✅ All error handling preserved
- ✅ No functional changes to script behavior

## Future Improvements

Potential enhancements to consider:

1. **Unit Tests**: Create test suite for `common_functions.sh`
2. **More Shared Functions**: Extract additional common patterns
3. **Configuration Management**: Centralize default config values
4. **Logging Framework**: Add standardized logging utilities
5. **Progress Indicators**: Shared progress bar/spinner functions

## Architecture Diagram

```
┌──────────────────────────────────────────────────────┐
│                  User Scripts                         │
├──────────────────────────────────────────────────────┤
│  run_fslmer_univariate.sh  │  run_meta_pipeline.sh   │
│         (40 lines)          │      (640 lines)        │
└─────────────┬──────────────┴──────────┬──────────────┘
              │                         │
              │    source (imports)     │
              └─────────┬───────────────┘
                        │
              ┌─────────▼──────────────┐
              │  common_functions.sh    │
              │     (250 lines)         │
              ├─────────────────────────┤
              │ • Error handling        │
              │ • Path validation       │
              │ • JSON parsing          │
              │ • R environment setup   │
              │ • Python detection      │
              │ • Effect name mapping   │
              │ • FreeSurfer tools      │
              └─────────────────────────┘
```

## Migration Notes

If you have custom scripts based on the old versions:

1. **No immediate action required** - Old scripts still work
2. **To adopt shared functions**:
   - Add `source "$(dirname "$0")/common_functions.sh"` to your script
   - Replace duplicate functions with calls to shared functions
   - Test thoroughly before committing

## Questions?

For questions about this refactoring:
- Review the inline documentation in `common_functions.sh`
- Compare old vs. new versions using `git diff`
- Test in a development environment first
