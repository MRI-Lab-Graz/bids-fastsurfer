# Validation Checks Documentation

## Overview
Both `bids_long_fastsurfer.sh` and `bids_fastsurfer.sh` now include comprehensive pre-flight validation checks that run before any processing begins.

## Validation Categories

### 1. **System Requirements**
- ✓ Checks for required commands: `singularity`, `jq`, `find`, `grep`
- ✓ Reports Singularity/Apptainer version
- ⚠️ Checks for NVIDIA GPU availability via `nvidia-smi`
- ⚠️ Reports GPU name, driver version, and memory
- ⚠️ Warns if GPU unavailable (will fall back to CPU)

### 2. **FastSurfer Container**
- ✓ Verifies SIF file exists
- ✓ Reports container file size
- ❌ Exits if SIF file missing

### 3. **FreeSurfer License**
- ✓ Verifies license file exists
- ✓ Basic validation of license format
- ⚠️ Warns if license format appears invalid
- ❌ Exits if license file missing

### 4. **JSON Configuration**
- ✓ Validates JSON syntax
- ✓ Checks for required keys:
  - `sif_file`
  - `fs_license`
  - `long` (for longitudinal script)
  - `cross` (for cross-sectional script)
- ❌ Exits if JSON invalid or missing required keys

### 5. **BIDS Directory Structure**
- ✓ Verifies BIDS root directory exists
- ✓ Counts and reports number of subjects found
- ❌ Exits if no subjects found (no `sub-*` directories)

### 6. **Output Directory**
- ✓ Verifies output directory exists
- ✓ Checks write permissions
- ✓ Reports available disk space
- ❌ Exits if directory missing or not writable

### 7. **Input Data Validation**
- ✓ Searches for T1w images in BIDS directory
- ✓ Reports count of T1w images found
- ✓ Supports multiple formats: `*_T1w.nii.gz`, `*_T1w.nii`, `*_desc-preproc_T1w.nii.gz`
- ❌ Exits if no T1w images found

## Error Handling

### Exit Codes
- **0**: All validations passed
- **1**: One or more validation errors detected

### Error Types
- ❌ **Critical Errors**: Stop execution immediately
- ⚠️ **Warnings**: Report issue but continue (e.g., no GPU available)

## Example Output

```bash
[VALIDATION] Checking system requirements...
  [OK] Found: singularity
  [OK] Found: jq
  [OK] Found: find
  [OK] Singularity version: apptainer version 1.4.3
  [OK] GPU detected: NVIDIA A2, 580.95.05, 15356 MiB
[VALIDATION] Checking FastSurfer container...
  [OK] SIF file exists: /path/to/fastsurfer.sif (3.8G)
[VALIDATION] Checking FreeSurfer license...
  [OK] License file exists: /path/to/license.txt
  [OK] License file appears valid
[VALIDATION] Checking JSON config...
  [OK] JSON config is valid
[VALIDATION] Checking BIDS directory structure...
  [OK] Found 150 subject(s) in BIDS directory
[VALIDATION] Checking output directory...
  [OK] Output directory writable: /output (Available: 1.2T)
[VALIDATION] Checking for T1w images...
  [OK] Found 432 T1w image(s)
[VALIDATION] ==========================================
[VALIDATION] ✓ All checks passed!
```

## Benefits

1. **Early Error Detection**: Catches configuration issues before processing starts
2. **Clear Feedback**: Detailed error messages help users fix problems quickly
3. **System Health**: Reports GPU status and disk space availability
4. **Time Saving**: Avoids wasting hours on processing that will fail due to setup issues
5. **Reproducibility**: Ensures all requirements are met before analysis begins

## When Validation Runs

Validation occurs:
- After argument parsing
- Before any subject processing
- Before background jobs are spawned
- Before batch processing triggers

This ensures all prerequisites are met regardless of processing mode (pilot, batch, nohup, etc.).
