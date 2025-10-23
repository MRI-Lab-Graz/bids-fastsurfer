#!/usr/bin/env bash
# validate_fastsurfer.sh
# Shared validation function for FastSurfer processing scripts
# Usage: source this file and call validate_requirements()

validate_requirements() {
  local bids_data="$1"
  local output_dir="$2"
  local config="$3"
  local sif_file="$4"
  local fs_license="$5"
  
  echo "[VALIDATION] Checking system requirements..."
  local validation_errors=0
  
  # 1. Check required commands
  for cmd in singularity jq find grep; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      echo "  [ERROR] Required command '$cmd' not found in PATH"
      validation_errors=$((validation_errors + 1))
    else
      echo "  [OK] Found: $cmd"
    fi
  done
  
  # 2. Check Singularity version
  if command -v singularity >/dev/null 2>&1; then
    sing_version=$(singularity --version 2>&1 | head -n1)
    echo "  [OK] Singularity version: $sing_version"
  fi
  
  # 3. Check NVIDIA/CUDA availability
  if command -v nvidia-smi >/dev/null 2>&1; then
    if nvidia-smi >/dev/null 2>&1; then
      gpu_info=$(nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader | head -n1)
      echo "  [OK] GPU detected: $gpu_info"
    else
      echo "  [WARN] nvidia-smi found but failed (driver issue?)"
      echo "  [WARN] FastSurfer will fall back to CPU (much slower!)"
    fi
  else
    echo "  [WARN] nvidia-smi not found - no GPU acceleration available"
    echo "  [WARN] FastSurfer will run on CPU (much slower!)"
  fi
  
  # 4. Check SIF file
  echo "[VALIDATION] Checking FastSurfer container..."
  if [[ ! -f "$sif_file" ]]; then
    echo "  [ERROR] Singularity image file not found: $sif_file"
    validation_errors=$((validation_errors + 1))
  else
    sif_size=$(du -h "$sif_file" | cut -f1)
    echo "  [OK] SIF file exists: $sif_file ($sif_size)"
  fi
  
  # 5. Check FreeSurfer license
  echo "[VALIDATION] Checking FreeSurfer license..."
  if [[ ! -f "$fs_license" ]]; then
    echo "  [ERROR] FreeSurfer license file not found: $fs_license"
    validation_errors=$((validation_errors + 1))
  else
    echo "  [OK] License file exists: $fs_license"
    # Check license content
    if grep -q "^[[:space:]]*#" "$fs_license" && grep -q "[A-Za-z0-9]" "$fs_license"; then
      echo "  [OK] License file appears valid"
    else
      echo "  [WARN] License file format may be invalid"
    fi
  fi
  
  # 6. Check JSON config validity
  echo "[VALIDATION] Checking JSON config..."
  if ! jq empty "$config" 2>/dev/null; then
    echo "  [ERROR] Invalid JSON in config file: $config"
    validation_errors=$((validation_errors + 1))
  else
    echo "  [OK] JSON config is valid"
    # Check required keys (check for either 'cross' or 'long')
    if ! jq -e '.sif_file' "$config" >/dev/null 2>&1; then
      echo "  [ERROR] Missing required key 'sif_file' in config"
      validation_errors=$((validation_errors + 1))
    fi
    if ! jq -e '.fs_license' "$config" >/dev/null 2>&1; then
      echo "  [ERROR] Missing required key 'fs_license' in config"
      validation_errors=$((validation_errors + 1))
    fi
    if ! jq -e '.cross' "$config" >/dev/null 2>&1 && ! jq -e '.long' "$config" >/dev/null 2>&1; then
      echo "  [ERROR] Missing required key 'cross' or 'long' in config"
      validation_errors=$((validation_errors + 1))
    fi
  fi
  
  # 7. Check BIDS directory structure
  echo "[VALIDATION] Checking BIDS directory structure..."
  if [[ ! -d "$bids_data" ]]; then
    echo "  [ERROR] BIDS input directory not found: $bids_data"
    validation_errors=$((validation_errors + 1))
  else
    subject_count=$(find "$bids_data" -maxdepth 1 -type d -name 'sub-*' 2>/dev/null | wc -l)
    if [[ $subject_count -eq 0 ]]; then
      echo "  [ERROR] No subjects found (no sub-* directories in $bids_data)"
      validation_errors=$((validation_errors + 1))
    else
      echo "  [OK] Found $subject_count subject(s) in BIDS directory"
    fi
  fi
  
  # 8. Check output directory
  echo "[VALIDATION] Checking output directory..."
  if [[ ! -d "$output_dir" ]]; then
    echo "  [ERROR] Output directory not found: $output_dir"
    validation_errors=$((validation_errors + 1))
  elif [[ ! -w "$output_dir" ]]; then
    echo "  [ERROR] Output directory not writable: $output_dir"
    validation_errors=$((validation_errors + 1))
  else
    available_space=$(df -h "$output_dir" | awk 'NR==2 {print $4}')
    echo "  [OK] Output directory writable: $output_dir (Available: $available_space)"
  fi
  
  # 9. Check for T1w images in BIDS directory
  echo "[VALIDATION] Checking for T1w images..."
  t1w_count=$(find "$bids_data" -type f \( -name "*_T1w.nii.gz" -o -name "*_T1w.nii" -o -name "*_desc-preproc_T1w.nii.gz" \) 2>/dev/null | wc -l)
  if [[ $t1w_count -eq 0 ]]; then
    echo "  [ERROR] No T1w images found in BIDS directory"
    validation_errors=$((validation_errors + 1))
  else
    echo "  [OK] Found $t1w_count T1w image(s)"
  fi
  
  # Summary
  echo "[VALIDATION] =========================================="
  if [[ $validation_errors -eq 0 ]]; then
    echo "[VALIDATION] ✓ All checks passed!"
    return 0
  else
    echo "[VALIDATION] ✗ Found $validation_errors error(s)"
    echo "[VALIDATION] Please fix the errors above before proceeding."
    return 1
  fi
}
