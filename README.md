# FastSurfer BIDS Wrappers

This repository contains two Bash wrappers to run the FastSurfer cross-sectional and longitudinal pipelines on BIDS datasets via Singularity containers. They add BIDS-aware discovery, JSON-based configuration, and quality-of-life flags like dry-run, debug, and pilot sampling.

- `bids_faststurfer.sh` — Cross-sectional wrapper (note the file name spelling)
- `bids_long_fastsurfer.sh` — Longitudinal wrapper (auto-detects longitudinal subjects by default)

The scripts are designed to behave like lightweight BIDS apps: you pass the input dataset, an output directory, and a config JSON file; the scripts discover inputs and build the correct Singularity commands for you.

## Requirements

- A Linux environment with Singularity installed (or Apptainer; if using Apptainer, you can alias `singularity` to `apptainer`). On macOS, run within a Linux VM or environment that provides Singularity/Apptainer.
- `jq` for JSON parsing
- FastSurfer GPU Singularity image (`.sif`)
- A valid FreeSurfer license file

Quick checks:

```zsh
# should print a version
singularity --version

# should print a path after install
jq --version
```

## Repository layout

- `bids_faststurfer.sh` — Cross-sectional wrapper
- `bids_long_fastsurfer.sh` — Longitudinal wrapper
- `fastsurfer_options.json` — Example configuration file with cross/long sections
- `license.txt` — Placeholder filename referenced by scripts (replace with your actual FreeSurfer license file path in the config)
- `scripts/install.sh` — One-command installer for the R analysis environment (micromamba)
- `scripts/generate_qdec.py` — Build Qdec from BIDS participants.tsv and subjects_dir; optional .long symlink helper
- `scripts/fslmer_univariate.R` — Flexible univariate LME/GAM/GLM helper (fslmer/mgcv/stats)
- `configs/fslmer_univariate.example.json` — Example JSON config for the R script

## Example BIDS datasets

If you need public BIDS datasets to try the wrappers, you can use these OpenNeuro examples (or any other BIDS-compatible dataset):

- Longitudinal: ds004937 v1.0.1 — <https://openneuro.org/datasets/ds004937/versions/1.0.1>
- Cross-sectional: ds004965 v1.0.1 — <https://openneuro.org/datasets/ds004965/versions/1.0.1>

Download using your preferred tool (OpenNeuro CLI, openneuro-py, DataLad, or web download) and point `<BIDS_ROOT>` to the dataset root.

### Download examples

You can download these datasets with either the OpenNeuro CLI or the Python package `openneuro-py`.

- OpenNeuro CLI (requires Node.js):

```zsh
# Install once (optional)
npm install -g @openneuro/cli

# Longitudinal dataset (ds004937 v1.0.1)
openneuro download --snapshot 1.0.1 --dataset ds004937 --destination /path/to/BIDS-ds004937

# Cross-sectional dataset (ds004965 v1.0.1)
openneuro download --snapshot 1.0.1 --dataset ds004965 --destination /path/to/BIDS-ds004965
```

- Python: openneuro-py (<https://pypi.org/project/openneuro-py/>)

```zsh
# Create/activate a virtual environment (optional)
python3 -m venv ~/.venvs/openneuro && source ~/.venvs/openneuro/bin/activate
pip install --upgrade pip openneuro-py

# Download longitudinal dataset ds004937 v1.0.1
python - <<'PY'
from openneuro import download
download(dataset='ds004937', target_dir='/path/to/BIDS-ds004937', include='*', tag='1.0.1')
PY

# Download cross-sectional dataset ds004965 v1.0.1
python - <<'PY'
from openneuro import download
download(dataset='ds004965', target_dir='/path/to/BIDS-ds004965', include='*', tag='1.0.1')
PY
```

After downloading, use the dataset root directories (`/path/to/BIDS-ds004937` or `/path/to/BIDS-ds004965`) as `<BIDS_ROOT>` in the examples below.

## Configuration: `fastsurfer_options.json`

Both scripts read a single JSON file with top-level container paths and per-pipeline option sections. Replace placeholder paths below with your actual locations.

```json
{
  "fs_license": "/abs/path/to/license.txt",
  "sif_file": "/abs/path/to/fastsurfer-gpu.sif",
  "cross": {
    "vox_size": "min",
    "seg_only": true,
    "3T": true,
    "reg_mode": "coreg",
    "qc_snap": false
    // ... more cross-sectional options (see below)
  },
  "long": {
    "parallel": null,
    "parallel_seg": null,
    "parallel_surf": null,
    "reg_mode": "coreg",
    "3T": true,
    "qc_snap": false
    // ... more longitudinal options (see below)
  }
}
```

- Notes:

- Top-level keys:
  - `fs_license` — Absolute path to your FreeSurfer license file.
  - `sif_file` — Absolute path to the FastSurfer GPU Singularity image.

- The `cross` and `long` sections mirror FastSurfer flags. Values are translated as follows:
  - Boolean `true` → adds `--flag`
  - Boolean `false` or `null` → omitted
  - Strings/numbers → adds `--flag value`

- Unknown key behavior differs:
  - Longitudinal (`long`): unknown keys are ignored; run with `--debug` to see warnings.
  - Cross-sectional (`cross`): keys are forwarded as flags to FastSurfer; unknown ones may cause errors. Keep the `cross` section minimal and only include supported flags.

Common option keys you can use in `cross`/`long`:
- Parallelization and threading: `parallel`, `parallel_seg`, `parallel_surf`, `threads`, `threads_seg`, `threads_surf`, `batch`
- Device selection: `device`, `viewagg_device` (GPU index if applicable)
- MRI specifics and modes: `3T`, `reg_mode`, `surf_only`, `no_fs_T1`, `no_surfreg`
- Additional toggles: `qc_snap`, `ignore_fs_version`, `fstess`, `fsqsphere`, `fsaparc`, `allow_root`, `base`

Tip: Keep only the options you actually need; remove or set unneeded keys to `null`.

## Configuration reference

This section lists commonly used keys you can set in `fastsurfer_options.json`. Not every key is required—omit what you don’t need.

- Top-level (required):
	- `fs_license` (string): Absolute path to your FreeSurfer license file (license.txt).
	- `sif_file` (string): Absolute path to your FastSurfer GPU `.sif` container.

- Cross (`cross` section) — forwarded to `/fastsurfer/run_fastsurfer.sh`:
	- `3T` (bool): Enable 3T-specific heuristics.
	- `seg_only` (bool): Run segmentation only (no surfaces).
	- `surf_only` (bool): Run surface pipeline only.
	- `reg_mode` (string): Registration mode, e.g. `coreg`.
	- `qc_snap` (bool): Generate QC snapshots.
	- `device`, `viewagg_device` (string/int): GPU selection if multiple GPUs.
	- `threads`, `threads_seg`, `threads_surf` (int): Thread counts.
	- `parallel`, `parallel_seg`, `parallel_surf` (bool): Enable parallel execution where supported.
	- `batch` (int): Batch size or job parallelism (depends on stage).
	- `ignore_fs_version` (bool): Skip FreeSurfer version checks.
	- `fstess`, `fsqsphere`, `fsaparc` (bool): Advanced FreeSurfer toggles.
	- `no_fs_T1`, `no_surfreg` (bool): Disable specific steps.
	- `vox_size` (string): Voxel size handling (e.g., `min`).
	- `no_biasfield`, `tal_reg`, `no_asegdkt`, `no_cereb`, `no_hypothal` (bool): Advanced toggles.
	- `t2` (string): Path to a T2 image (container path). Note: the wrapper also tries to auto-detect `--t2`; avoid setting both to prevent duplicates.
	- Avoid including runtime flags the wrapper sets: `t1`, `sid`, `sd`, `fs_license`, `sif_file`, `py`.

- Longitudinal (`long` section) — filtered and mapped to `/fastsurfer/long_fastsurfer.sh`:
	- `parallel`, `parallel_seg`, `parallel_surf` (bool): Enable parallel execution.
	- `reg_mode` (string): Registration mode (e.g., `coreg`).
	- `qc_snap` (bool): Generate QC snapshots.
	- `surf_only` (bool): Surface pipeline only.
	- `3T` (bool): 3T-specific heuristics.
	- `device`, `viewagg_device` (string/int): GPU selection.
	- `threads`, `threads_seg`, `threads_surf` (int): Thread counts.
	- `batch` (int): Batch size or job parallelism.
	- `ignore_fs_version` (bool): Skip FreeSurfer version checks.
	- `fstess`, `fsqsphere`, `fsaparc` (bool): Advanced FreeSurfer toggles.
	- `no_fs_T1`, `no_surfreg` (bool): Disable specific steps.
	- `allow_root` (bool): Allow running as root inside the container.
	- `base` (string): Advanced: base name control (used by some workflows).

For any option not listed here, consult FastSurfer’s upstream documentation. If an option isn’t recognized by the underlying scripts, it will be ignored in the longitudinal wrapper and may error in the cross wrapper.

## Cross-sectional wrapper: `bids_faststurfer.sh`

Run FastSurfer cross-sectionally across a BIDS dataset, a specific subject, or a specific session across all subjects.

Usage:

```zsh
bash bids_faststurfer.sh <BIDS_ROOT> <OUTPUT_DIR> -c <config.json> [--sub <sub-XXX>] [--ses <ses-YYY>] [--pilot] [--dry_run] [--debug]
```

- `<BIDS_ROOT>`: Path to your BIDS dataset root.
- `<OUTPUT_DIR>`: Destination directory for FastSurfer outputs (bind-mounted at `/output`).
- `-c/--config`: Path to your JSON config file (see above).
- `--sub sub-XXX`: Restrict to a single subject.
- `--ses ses-YYY`: Restrict to a session; can be used without `--sub` to process that session for all subjects.
- `--pilot`: Randomly picks a single T1w to test the pipeline quickly.
- `--dry_run`: Prints the constructed Singularity command(s) without executing.
- `--debug`: Prints extra debug information (resolved paths, options, etc.).

What it does:
- Scans for T1w images (`*_T1w.nii.gz`, `*_T1w.nii`, or `*_desc-preproc_T1w.nii.gz`) respecting `--sub/--ses` limits if provided.
- Optionally attempts to find a matching T2w (`*_T2w.nii.gz`) per T1w and adds `--t2` if found.
- Builds and executes the FastSurfer cross-sectional pipeline with the appropriate `--t1`, `--sid`, and other options derived from your BIDS dataset and config.
- The wrapper takes care of container bindings internally; you only provide host paths for `<BIDS_ROOT>` and `<OUTPUT_DIR>`.

Examples:

```zsh
# Dry run across whole BIDS dataset
bash bids_faststurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer \
	-c fastsurfer_options.json --dry_run

# Only subject sub-001 session ses-01
bash bids_faststurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer \
	-c fastsurfer_options.json --sub sub-001 --ses ses-01

# Session-only across all subjects
bash bids_faststurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer \
	-c fastsurfer_options.json --ses ses-01 --dry_run

# Pilot run (pick one random T1w)
bash bids_faststurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer \
	-c fastsurfer_options.json --pilot --dry_run --debug
```

Notes:
- Use a proper BIDS tree (e.g., `sub-XXX/ses-YYY/anat/`). The OpenNeuro examples above are already BIDS-compliant.

## Longitudinal wrapper: `bids_long_fastsurfer.sh`

Run FastSurfer’s longitudinal pipeline either automatically (default) or manually.

Usage (auto mode, default):

```zsh
bash bids_long_fastsurfer.sh <BIDS_ROOT> <OUTPUT_DIR> -c <config.json> [--pilot] [--dry_run] [--debug]
```

Usage (manual mode):

```zsh
bash bids_long_fastsurfer.sh <BIDS_ROOT> <OUTPUT_DIR> -c <config.json> \
	--tid <subject> --tpids <sub-XXX_ses-YYY> [<sub-XXX_ses-ZZZ> ...] [--py <python>] [--dry_run] [--debug]
```

- Auto mode:
	- Scans each `sub-*` for `ses-*` and requires at least 2 sessions with valid T1w in `anat/`.
	- Builds timepoint lists and runs `long_fastsurfer.sh` per eligible subject.
	- `--pilot` randomly selects one eligible subject to process.
- Manual mode:
	- `--tid` is the template subject (with or without `sub-` prefix).
	- `--tpids` are timepoint IDs like `sub-001_ses-01` (one or more).
	- The script locates the T1w for each TPID under `sub-XXX/ses-YYY/anat/`.
- `--py` lets you choose the Python executable inside the container (default: `python3`).

What it does:
- Builds and executes the FastSurfer longitudinal pipeline with the appropriate `--tid`, `--t1s`, `--tpids`, and options based on your dataset and config. The wrapper handles container bindings internally.

Examples:

```zsh
# Auto-detect all longitudinal subjects (>= 2 sessions)
bash bids_long_fastsurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer_long \
	-c fastsurfer_options.json --dry_run

# Auto + pilot (random one subject)
bash bids_long_fastsurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer_long \
	-c fastsurfer_options.json --pilot --dry_run

# Manual subject with 2 sessions
bash bids_long_fastsurfer.sh /path/to/BIDS /path/to/derivatives/fastsurfer_long \
	-c fastsurfer_options.json --tid sub-001 \
	--tpids sub-001_ses-01 sub-001_ses-02 --dry_run --debug
```

 Notes:
- `--pilot` is only valid in auto mode.
- Unknown keys in the `long` section are ignored, with warnings printed in `--debug` mode.

## Outputs

- All outputs go to `<OUTPUT_DIR>`.
- Cross-sectional runs create one FreeSurfer subject directory per input (named by `--sid`).
- Longitudinal runs create base and TP directories under the same `<OUTPUT_DIR>` according to FastSurfer’s conventions.

## Troubleshooting

- “singularity: command not found”
	- Install Singularity/Apptainer or run inside a suitable environment; on Apptainer, you can create an alias: `alias singularity=apptainer`.
- “jq: command not found”
	- Install jq (e.g., `sudo apt-get install jq`).
- “Config or license not found”
	- Ensure `fastsurfer_options.json` exists and contains valid absolute paths for `fs_license` and `sif_file`.
- “No T1w found”
	- Verify your BIDS tree and that `anat/*_T1w.nii[.gz]` files exist for the targeted subjects/sessions.
- Permission issues writing to output
	- Use an output directory you own; avoid binding network or read-only locations.
- GPU access inside Singularity
	- These wrappers pass `--nv`. Ensure the host has NVIDIA drivers and the container supports GPU.

### Headless servers and fsqc surfaces

- On servers without a graphical display (`DISPLAY` not set), the fsqc "surfaces" module requires OpenGL and will fail with a GLFW error. Our prep script automatically disables fsqc surfaces when it detects no DISPLAY and will log an info message. You can still generate headless-friendly outputs by enabling screenshots/metrics (e.g., `--qc --qc-screenshots --qc-html`).

## Tips

- Start with `--dry_run --debug` to inspect the constructed commands.
- Use `--pilot` to quickly validate the environment on a small subset before scaling up.
- Keep your config minimal and explicit; remove unused keys to reduce ambiguity.

## Longitudinal statistics: generate Qdec

To run FreeSurfer-style longitudinal statistics (e.g., with the R package fslmer), generate a Qdec file from your BIDS `participants.tsv` and your FastSurfer/FreeSurfer subjects directory.

Script: `scripts/generate_qdec.py`

Inputs:
- `--participants`: Path to `participants.tsv` with at least `participant_id` and optionally `session_id`, plus covariates (`age`, `sex`, `group`, ...).
- `--subjects-dir`: Path to the FastSurfer/FreeSurfer subjects directory. Expect directories like `sub-1291003` (base) and `sub-1291003_ses-1` (timepoint).
- `--output`: Output path for the Qdec TSV (default: `qdec.table.dat`).
- `--include-columns`: Optional explicit list of covariates to include from `participants.tsv`. If omitted, all columns except participant/session are included.
- `--strict`: If set, the script fails when a timepoint has no matching participants row; otherwise fills `n/a`.
- `--inspect`: Print `participants.tsv` column names and exit (useful for column selection).
- `--bids`: Optional BIDS root; if provided, the script prints a consistency summary comparing `participants.tsv`, `subjects_dir`, and the BIDS tree.
- `--list-limit`: Maximum number of subject IDs to print per list in the summary (default: 20). The summary avoids printing the same list twice when it’s identical across sources.

Output columns:

- `fsid` — timepoint subject id (e.g., `sub-1291003_ses-1`)
- `fsid-base` — base/template id (e.g., `sub-1291003`)
- `tp` — numeric timepoint derived from `ses-<number>` (or `n/a` if not parseable)
- Covariates — selected from `participants.tsv`

Example:

```zsh
python scripts/generate_qdec.py \
  --participants /path/to/BIDS/participants.tsv \
  --subjects-dir /path/to/derivatives/fastsurfer_long \
  --output /path/to/qdec.table.dat \
  --include-columns age sex group
```

You can also inspect columns or print a consistency summary:

```zsh
# See columns in participants.tsv
python scripts/generate_qdec.py --participants /path/to/participants.tsv --subjects-dir /path/to/subjects --inspect

# Generate Qdec and print a summary of overlaps/missing subjects across sources
python scripts/generate_qdec.py --participants /path/to/participants.tsv --subjects-dir /path/to/subjects --bids /path/to/BIDS
```

Use the resulting `qdec.table.dat` with FreeSurfer longitudinal statistics and tools like `fslmer`:

- FreeSurfer Longitudinal: <https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalStatistics>
- fslmer: <https://github.com/Deep-MI/fslmer>

## Univariate stats in R (fslmer / GAM / GLM) 🎓

A flexible R helper is included to run univariate models on aseg/aparc tables using:
- Mixed-effects via `fslmer`
- Generalized Additive Models via `mgcv`
- Generalized Linear Models via base `stats`

Script: `scripts/fslmer_univariate.R`

Basic usage:

```zsh
Rscript scripts/fslmer_univariate.R \
  --qdec /path/to/qdec.table.dat \
  --aseg /path/to/aseg.long.table \
  --roi Left.Hippocampus \
  --formula "~ tp*group" \
  --zcols 1,2 \
  --contrast 0,0,0,1 \
  --outdir results_univar --save-merged
```

Flags:

- `--qdec`: Qdec TSV generated by this repo.
- `--aseg`: aseg/aparc table made with `asegstats2table --qdec-long` or `aparcstats2table --qdec-long`.
- `--roi`: ROI column name (R uses `.` instead of `-`, e.g., `Left.Hippocampus`).
- `--formula`: R fixed-effect formula (default `~ tp`). Add covariates/interaction as needed.
- `--zcols`: Random-effect columns (in design matrix order). `1,2` means random intercept and random slope for time.
- `--contrast`: Comma-separated vector matching `ncol(X)` to test a specific effect (optional).
- `--time-col`: If your time variable is not `tp` or `Time.From.Baseline`, set it here.
- `--id-col`: If the ID column in aseg is not `Measure:volume`, specify it.
- `--print-cols`: Print column names of qdec and aseg and exit.
- `--outdir`: Output directory (default `fslmer_out`).
- `--save-merged`: Save merged design/response as CSV.

Outputs in `--outdir`:

- `fit.rds` — model fit from `lme_fit_FS` (RDS).
- `Bhat.txt` — fixed-effect coefficients.
- `F_test.txt` — if `--contrast` given, F-values, p-values, sign, and df.
- `merged_data.csv` — optional merged data dump.

This maps closely to the fslmer README tutorial while automating the data loading/merging and making the model specification configurable via flags.

Config-driven usage and reproducibility:

- Instead of passing many flags, you can provide a JSON config via `--config` (see `configs/fslmer_univariate.example.json`). Any CLI flag overrides the same field in the JSON. The script writes the merged settings to `<outdir>/used_config.json`.
- Example JSON fields: `qdec`, `aseg`, `roi`, `formula`, `zcols` (array of integers), `contrast` (array), `outdir`, `time_col`, `id_col`.
- To run: point `--config` to your JSON file; optionally add `--outdir` or other flags to override.

Additional flag:

- `--config`: Path to a JSON file with the same keys as the CLI flags. Useful to version-control analyses.

## Installation (one command)

Set up everything needed for the R analysis helper with micromamba:

```zsh
bash scripts/install.sh
```

What it does:
- Bootstraps micromamba under `~/.local/micromamba` if missing
- Creates/updates an env (default `fastsurfer-r`) from `scripts/environment.yml` in strict mode
- Installs/validates R packages in the right order (checkmate/backports → bettermc → fslmer; plus mgcv, optparse, jsonlite)

Activate, verify, and use:

```zsh
# Activate the env
source scripts/mamba_activate.sh

# Verify packages
Rscript -e 'pkgs <- c("optparse","jsonlite","mgcv","checkmate","bettermc","fslmer"); print(sapply(pkgs, requireNamespace, quietly=TRUE))'

# See helper usage
Rscript scripts/fslmer_univariate.R --help

# Deactivate later
source scripts/mamba_deactivate.sh
```

Advanced install options:

```zsh
# Use spec-based installation (alternative to YAML), skip compilers on macOS
bash scripts/install.sh --use-specs --no-compilers

# Choose a custom env name and R version
bash scripts/install.sh --env my-fast-env --r 4.4
```

### Using FreeSurfer tools with FastSurfer longitudinal outputs

Some FreeSurfer utilities (e.g., `asegstats2table --qdec-long`) expect timepoints to be arranged as `<fsid>.long.<fsid-base>/stats/aseg.stats`. FastSurfer’s longitudinal pipeline may not create these `.long.*` directories by default.

To make FreeSurfer tools work without re-running with FreeSurfer, this repo provides a safe symlink workflow:

```zsh
# Preview/verify which links would be needed (no changes)
python scripts/generate_qdec.py \
  --participants /path/to/participants.tsv \
  --subjects-dir /path/to/fastsurfer_subjects \
  --verify-long --list-limit 10

# Create .long.<base> symlinks that point to each timepoint directory (dry-run)
python scripts/generate_qdec.py \
  --participants /path/to/participants.tsv \
  --subjects-dir /path/to/fastsurfer_subjects \
  --link-long --link-dry-run

# Actually create/update links (careful):
python scripts/generate_qdec.py \
  --participants /path/to/participants.tsv \
  --subjects-dir /path/to/fastsurfer_subjects \
  --link-long

# If an existing symlink points elsewhere, allow updating it:
python scripts/generate_qdec.py \
  --participants /path/to/participants.tsv \
  --subjects-dir /path/to/fastsurfer_subjects \
  --link-long --link-force
```

After creating the links, `asegstats2table --qdec-long` will be able to find `stats/aseg.stats` under the expected names.

## Credits and License

- FastSurfer is developed by the FastSurfer team; see their documentation for details and licensing.
- This repository’s scripts were created for the 2025 workshop. See `license.txt` for terms applying to this repo.
