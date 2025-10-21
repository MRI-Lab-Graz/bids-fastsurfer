# Changelog: .long Symlink Creation Fix

## Issue
After FastSurfer longitudinal processing completed successfully, subjects were still showing as "missing" because the `.long` symlinks were not being created. These symlinks are required by FreeSurfer tools like `asegstats2table --qdec-long`.

## Root Cause
The `bids_long_fastsurfer.sh` script checked for `.long` directories to determine if processing was complete (line 339), but never actually created these symlinks after processing finished.

## Solution

### 1. Added `create_long_symlinks()` Helper Function (Line 105-146)
- Creates symlinks for all timepoints of a template subject
- Format: `sub-XXX_ses-Y.long.sub-XXX` → `sub-XXX_ses-Y`
- Reports status: [OK] for created, [SKIP] for existing, [WARN] for missing dirs
- Only runs when not in dry-run mode

### 2. Integrated into Manual Mode (Line 418)
- Calls `create_long_symlinks()` after successful processing (exit code 0)
- Automatically creates symlinks for all timepoints of the specified template subject

### 3. Integrated into Auto Mode (Line 536)
- Calls `create_long_symlinks()` after successful processing for each subject
- Automatically creates symlinks during batch processing (non-nohup)

### 4. Added Warning for Nohup Mode (Line 549-551)
- Background jobs cannot automatically create symlinks (job returns immediately)
- Users are instructed to run the repair script after jobs complete:
  ```bash
  bash scripts/create_missing_long_symlinks.sh <OUTPUT_DIR>
  ```

### 5. Created Standalone Repair Script
- **Location**: `scripts/create_missing_long_symlinks.sh`
- **Purpose**: Creates missing symlinks for already-processed subjects
- **Usage**: `bash scripts/create_missing_long_symlinks.sh <OUTPUT_DIR>`
- **Features**:
  - Scans for template directories (those with `base-tps.fastsurfer`)
  - Finds all session directories for each template
  - Creates symlinks only for directories that exist
  - Reports summary statistics

## Verification
Tested on existing processed data:
- 69 template subjects scanned
- 203 existing symlinks verified
- 0 new symlinks created (all already existed from Oct 1 processing)

## Usage

### For New Processing (Non-Background)
Symlinks are created automatically:
```bash
bash bids_long_fastsurfer.sh --auto
```

### For New Processing (Background with --nohup)
Run the repair script after jobs complete:
```bash
bash bids_long_fastsurfer.sh --auto --nohup
# Wait for jobs to complete, then:
bash scripts/create_missing_long_symlinks.sh /path/to/OUTPUT_DIR
```

### For Existing Processed Data
Run the repair script once:
```bash
bash scripts/create_missing_long_symlinks.sh /path/to/OUTPUT_DIR
```

## Files Modified
1. `bids_long_fastsurfer.sh` - Added automatic symlink creation
2. `scripts/create_missing_long_symlinks.sh` - New standalone repair script

## Testing Recommendations
1. Test manual mode with a small subject to verify automatic creation
2. Run repair script on existing processed data to fix historical subjects
3. For nohup users: add repair script to workflow documentation
