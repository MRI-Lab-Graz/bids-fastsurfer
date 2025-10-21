#!/usr/bin/env bash
#
# create_missing_long_symlinks.sh
# 
# Creates missing .long symlinks for completed FastSurfer longitudinal processing.
# FastSurfer creates template directories (sub-XXX) and timepoint directories (sub-XXX_ses-Y),
# but the .long symlinks (sub-XXX_ses-Y.long.sub-XXX) need to be created for FreeSurfer
# tools like asegstats2table --qdec-long to work properly.
#
# Usage:
#   bash create_missing_long_symlinks.sh <OUTPUT_DIR>
#
# Example:
#   bash create_missing_long_symlinks.sh /data/local/129_PK01/derivatives/fastsurfer
#

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <OUTPUT_DIR>"
  echo ""
  echo "Creates .long symlinks for completed FastSurfer longitudinal processing."
  exit 1
fi

OUTPUT_DIR="$1"

if [[ ! -d "$OUTPUT_DIR" ]]; then
  echo "Error: Output directory '$OUTPUT_DIR' does not exist."
  exit 1
fi

cd "$OUTPUT_DIR"

echo "Scanning for template directories in: $OUTPUT_DIR"
echo ""

created_count=0
skipped_count=0
total_templates=0

# Find all template directories (directories without _ses- in name and with at least 2 session subdirectories)
for base_dir in sub-*; do
  [[ ! -d "$base_dir" ]] && continue
  
  # Skip if this is a timepoint directory (contains _ses-)
  if [[ "$base_dir" =~ _ses- ]]; then
    continue
  fi
  
  # Check if this has a base-tps.fastsurfer file (indicates it's a template)
  if [[ ! -f "$base_dir/base-tps.fastsurfer" ]]; then
    continue
  fi
  
  total_templates=$((total_templates + 1))
  
  # Find all session directories for this template
  mapfile -t session_dirs < <(find . -maxdepth 1 -type d -name "${base_dir}_ses-*" -exec basename {} \; | sort)
  
  if [[ ${#session_dirs[@]} -lt 2 ]]; then
    echo "[SKIP] $base_dir: only ${#session_dirs[@]} session(s) found"
    continue
  fi
  
  echo "[TEMPLATE] $base_dir (${#session_dirs[@]} sessions)"
  
  # Create symlinks for each session
  for sess_dir in "${session_dirs[@]}"; do
    long_link="${sess_dir}.long.${base_dir}"
    
    if [[ -e "$long_link" ]]; then
      echo "  ✓ $long_link (already exists)"
      skipped_count=$((skipped_count + 1))
    else
      if ln -s "$sess_dir" "$long_link"; then
        echo "  ✓ Created: $long_link -> $sess_dir"
        created_count=$((created_count + 1))
      else
        echo "  ✗ Failed to create: $long_link"
      fi
    fi
  done
  echo ""
done

echo "=========================================="
echo "Summary:"
echo "  Templates found: $total_templates"
echo "  Symlinks created: $created_count"
echo "  Symlinks already existed: $skipped_count"
echo "=========================================="
