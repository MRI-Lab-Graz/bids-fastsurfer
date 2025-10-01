#!/usr/bin/env python3
"""
Prepare longitudinal analysis inputs:

- Generate a FreeSurfer longitudinal Qdec file from a BIDS participants.tsv and a
  FastSurfer/FreeSurfer subjects directory.
- Optionally verify/create .long symlinks for FreeSurfer compatibility.
- Optionally run asegstats2table and/or aparcstats2table in longitudinal mode.

Usage examples:

  1. Basic usage - generate QDEC table only:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --output results/

  2. Generate QDEC + aseg and aparc tables:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --output results/ \
       --aseg --aparc

  3. Include surface preprocessing with custom smoothing:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --output results/ \
       --surf --smooth 5,10,15 \
       --surf-measures thickness,area

  4. Enable quality control (QC) with screenshots:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --output results/ \
       --qc --qc-screenshots --qc-html

  5. Use config file for settings:
     python scripts/prep_long.py \
       --config configs/prep_long.example.json

  6. Create FreeSurfer-compatible .long symlinks:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --verify-long --link-long

  7. Skip specific subjects:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --skip-sub "sub-01,sub-02" \
       --output results/

  8. Include only specific columns from participants.tsv:
     python scripts/prep_long.py \
       --participants data/participants.tsv \
       --subjects-dir data/subjects \
       --include-columns age sex group \
       --output results/

Notes:
- SUBJECTS_DIR is set for the aseg/aparc table commands to the provided subjects-dir.
- By default, aparc uses parc=aparc, hemis=lh,rh, measures=thickness,area,volume.
- Surface preprocessing (--surf) is enabled by default; use --no-surf to disable.
- Aseg and aparc table generation are enabled by default; use --no-aseg/--no-aparc to disable.
"""

from __future__ import annotations

import argparse
import json
import os
import csv
import re
import sys
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set


SUBJECT_DIR_PATTERN = re.compile(r"^(?P<base>sub-[^/]+?)(?:_(?P<ses>ses-[^/]+))?$")
SES_NUM_PATTERN = re.compile(r"^ses-(?P<num>\d+)$")


def _coerce_list(val) -> Optional[List[str]]:
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        return [str(x) for x in val]
    # split on comma or whitespace
    s = str(val).strip()
    if not s:
        return []
    if "," in s:
        return [x.strip() for x in s.split(",") if x.strip()]
    return [x for x in s.split() if x]


def _coerce_int_list(val) -> Optional[List[int]]:
    lst = _coerce_list(val)
    if lst is None:
        return None
    out: List[int] = []
    for x in lst:
        try:
            out.append(int(x))
        except Exception:
            continue
    return out


def check_dependencies(args: argparse.Namespace) -> List[str]:
    """Check for required dependencies and return list of missing tools/packages.
    
    Returns:
        List of missing dependencies with installation instructions.
    """
    missing = []
    
    # Check FreeSurfer tools
    if args.aseg and not args.link_dry_run:
        if shutil.which("asegstats2table") is None:
            missing.append("asegstats2table (FreeSurfer) - ensure FreeSurfer is sourced")
    
    if args.aparc and not args.link_dry_run:
        if shutil.which("aparcstats2table") is None:
            missing.append("aparcstats2table (FreeSurfer) - ensure FreeSurfer is sourced")
    
    if args.surf:
        if shutil.which("mris_preproc") is None:
            missing.append("mris_preproc (FreeSurfer) - ensure FreeSurfer is sourced")
        if shutil.which("mri_surf2surf") is None:
            missing.append("mri_surf2surf (FreeSurfer) - ensure FreeSurfer is sourced")
    
    # Check Python packages
    if args.qc:
        if shutil.which("run_fsqc") is None:
            missing.append("fsqc (Python package) - install with: bash scripts/install.sh (or pip install fsqc)")
    
    return missing


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    # Phase 1: parse only --config to optionally load defaults from JSON
    p0 = argparse.ArgumentParser(add_help=False)
    p0.add_argument("--config", type=Path, default=None, help="Optional JSON config file with settings; CLI overrides it")
    ns0, _ = p0.parse_known_args(argv)

    cfg: Dict[str, object] = {}
    if ns0.config is not None:
        cfg_path = Path(ns0.config)
        if not cfg_path.exists():
            raise FileNotFoundError(f"Config file not found: {cfg_path}")
        with cfg_path.open("r") as fh:
            cfg = json.load(fh)

    p = argparse.ArgumentParser(description="Prepare longitudinal QDEC and optional aseg/aparc tables", parents=[p0])
    
    # Required arguments
    required = p.add_argument_group('Required arguments')
    required.add_argument("--participants", required=True, type=Path, help="Path to BIDS participants.tsv")
    required.add_argument(
        "--subjects-dir",
        required=True,
        type=Path,
        help="Path to FastSurfer/FreeSurfer subjects directory",
    )
    
    # Input/Output configuration
    io_group = p.add_argument_group('Input/Output configuration')
    io_group.add_argument(
        "--output",
        type=Path,
        default=Path("results"),
        help=(
            "Output directory (default: results). The QDEC file 'qdec.table.dat' and other outputs"
            " (aseg/aparc tables, surf, fsqc) will be written under this folder. For backwards"
            " compatibility, a file path ending in .dat/.tsv/.table is also accepted."
        ),
    )
    io_group.add_argument("--participant-column", default="participant_id", help="Column name for participant id (default: participant_id)")
    io_group.add_argument("--session-column", default="session_id", help="Column name for session id if present (default: session_id)")
    io_group.add_argument(
        "--include-columns",
        nargs="*",
        default=None,
        help="Optional explicit list of covariate columns to include from participants.tsv. "
             "If omitted, include all columns except participant and session columns.",
    )
    io_group.add_argument("--strict", action="store_true", help="Fail if a subjects_dir timepoint has no matching participants row")
    io_group.add_argument("--inspect", action="store_true", help="Print participants.tsv columns and exit")
    io_group.add_argument("--bids", type=Path, default=None, help="Optional BIDS root to cross-check subjects/sessions consistency")
    io_group.add_argument("--list-limit", type=int, default=20, help="Max number of IDs to show when listing missing subjects (default: 20)")
    io_group.add_argument("--force", action="store_true", help="Overwrite/replace existing outputs where applicable (surf, qc, tables)")
    
    # FreeSurfer .long compatibility
    long_group = p.add_argument_group('FreeSurfer .long compatibility')
    long_group.add_argument("--verify-long", action="store_true", help="Verify presence of <fsid>.long.<base>/stats/aseg.stats for each timepoint")
    long_group.add_argument("--link-long", action="store_true", help="Create <fsid>.long.<base> symlinks pointing to the timepoint directories when missing")
    long_group.add_argument("--link-dry-run", action="store_true", help="Print the symlink actions without making changes")
    long_group.add_argument("--link-force", action="store_true", help="If an existing symlink points elsewhere, replace it (does not delete real directories)")
    
    # Statistical tables (aseg/aparc)
    tables_group = p.add_argument_group('Statistical tables (aseg/aparc)')
    tables_group.add_argument("--aseg", dest="aseg", action="store_true", help="Enable asegstats2table (default: enabled)")
    tables_group.add_argument("--no-aseg", dest="aseg", action="store_false", help="Disable asegstats2table")
    tables_group.add_argument("--aparc", dest="aparc", action="store_true", help="Enable aparcstats2table (default: enabled; auto-detects parcellation)")
    tables_group.add_argument("--no-aparc", dest="aparc", action="store_false", help="Disable aparcstats2table")
    p.set_defaults(aseg=True, aparc=True)
    tables_group.add_argument("--aparc-parc", default="aparc", help="Aparc parcellation preference (default: aparc; auto-fallback tries aparc.DKTatlas.mapped, aparc, aparc.a2009s)")
    tables_group.add_argument("--aparc-measures", nargs="*", default=["thickness", "area", "volume"], help="Measures for aparcstats2table (default: thickness area volume)")
    tables_group.add_argument("--aparc-hemis", nargs="*", default=["lh", "rh"], help="Hemispheres for aparcstats2table (default: lh rh)")
    
    # Subject filtering
    filter_group = p.add_argument_group('Subject filtering')
    filter_group.add_argument("--skip-sub", default=None, help="Comma-separated fsid_base IDs to skip (exclude) from QDEC")
    filter_group.add_argument("--skip-file", type=Path, default=None, help="File with fsid_base IDs (one per line) to skip from QDEC")

    # Surface preprocessing
    surf_group = p.add_argument_group('Surface preprocessing')
    surf_group.add_argument("--surf", dest="surf", action="store_true", help="Enable mass-univariate surface prep (mris_preproc + mri_surf2surf); default enabled")
    surf_group.add_argument("--no-surf", dest="surf", action="store_false", help="Disable mass-univariate surface prep")
    surf_group.add_argument("--surf-target", default="fsaverage", help="Surface template target for mris_preproc/mri_surf2surf (default: fsaverage)")
    surf_group.add_argument("--surf-measures", nargs="*", default=["thickness"], help="Surface measures to prepare (default: thickness)")
    surf_group.add_argument("--surf-hemis", nargs="*", default=["lh", "rh"], help="Hemispheres to prepare (default: lh rh)")
    surf_group.add_argument("--smooth", default=None, help="Comma- or space-separated smoothing kernels in mm (e.g., '5,10,15'). Overrides --surf-fwhm")
    surf_group.add_argument("--surf-fwhm", type=int, default=10, help="[Deprecated] Single smoothing FWHM (mm) if --smooth not provided (default: 10)")
    surf_group.add_argument("--surf-outdir", type=Path, default=None, help="Output directory for surface files (default: alongside QDEC under 'surf/')")
    p.set_defaults(surf=True)
    
    # Quality control (QC)
    qc_group = p.add_argument_group('Quality control (QC)')
    qc_group.add_argument("--qc", dest="qc", action="store_true", help="Run fsqc on selected subjects (requires 'run_fsqc' in PATH)")
    qc_group.add_argument("--no-qc", dest="qc", action="store_false", help="Disable fsqc step")
    p.set_defaults(qc=False)
    qc_group.add_argument("--qc-output", type=Path, default=None, help="fsqc output directory (default: alongside QDEC under 'fsqc/')")
    qc_group.add_argument("--qc-from", choices=["fsid", "base"], default="fsid", help="Select subject IDs for fsqc from QDEC: timepoints (fsid) or bases (fsid-base)")
    qc_group.add_argument("--qc-fastsurfer", action="store_true", help="Tell fsqc to use --fastsurfer (default: on)")
    qc_group.add_argument("--qc-no-fastsurfer", dest="qc_fastsurfer", action="store_false")
    p.set_defaults(qc_fastsurfer=True)
    qc_group.add_argument("--qc-screenshots", action="store_true", help="Enable fsqc screenshots module (and HTML if --qc-html)")
    qc_group.add_argument("--qc-surfaces", action="store_true", help="Enable fsqc surfaces module (and HTML if --qc-html)")
    qc_group.add_argument("--qc-skullstrip", action="store_true", help="Enable fsqc skullstrip module (and HTML if --qc-html)")
    qc_group.add_argument("--qc-outlier", action="store_true", help="Enable fsqc outlier detection module")
    qc_group.add_argument("--qc-html", action="store_true", help="For enabled modules, also produce HTML summary pages")
    qc_group.add_argument("--qc-skip-existing", action="store_true", help="Pass --skip-existing to fsqc to avoid recomputation")
    # Apply config defaults if provided
    if cfg:
        # Map config keys directly to argparse dests; allow hyphen or underscore forms
        def get_cfg(key: str, default=None):
            if key in cfg:
                return cfg[key]
            alt = key.replace("_", "-")
            if alt in cfg:
                return cfg[alt]
            alt2 = key.replace("-", "_")
            if alt2 in cfg:
                return cfg[alt2]
            return default

        # Pre-seed defaults from config; lists may be given as list or comma strings
        if get_cfg("participants"):
            p.set_defaults(participants=Path(get_cfg("participants")))
        if get_cfg("subjects_dir") or get_cfg("subjects-dir"):
            p.set_defaults(subjects_dir=Path(get_cfg("subjects_dir") or get_cfg("subjects-dir")))
        if get_cfg("output"):
            p.set_defaults(output=Path(get_cfg("output")))
        if get_cfg("participant_column"):
            p.set_defaults(participant_column=str(get_cfg("participant_column")))
        if get_cfg("session_column"):
            p.set_defaults(session_column=str(get_cfg("session_column")))
        if get_cfg("include_columns") is not None:
            p.set_defaults(include_columns=_coerce_list(get_cfg("include_columns")))
        if get_cfg("strict") is not None:
            p.set_defaults(strict=bool(get_cfg("strict")))
        if get_cfg("bids"):
            p.set_defaults(bids=Path(get_cfg("bids")))
        if get_cfg("list_limit") is not None:
            p.set_defaults(list_limit=int(get_cfg("list_limit")))
        if get_cfg("force") is not None:
            p.set_defaults(force=bool(get_cfg("force")))
        if get_cfg("verify_long") is not None:
            p.set_defaults(verify_long=bool(get_cfg("verify_long")))
        if get_cfg("link_long") is not None:
            p.set_defaults(link_long=bool(get_cfg("link_long")))
        if get_cfg("link_dry_run") is not None:
            p.set_defaults(link_dry_run=bool(get_cfg("link_dry_run")))
        if get_cfg("link_force") is not None:
            p.set_defaults(link_force=bool(get_cfg("link_force")))
        if get_cfg("aseg") is not None:
            p.set_defaults(aseg=bool(get_cfg("aseg")))
        if get_cfg("aparc") is not None:
            p.set_defaults(aparc=bool(get_cfg("aparc")))
        if get_cfg("aparc_parc"):
            p.set_defaults(aparc_parc=str(get_cfg("aparc_parc")))
        if get_cfg("aparc_measures") is not None:
            p.set_defaults(aparc_measures=_coerce_list(get_cfg("aparc_measures")) or ["thickness", "area", "volume"])
        if get_cfg("aparc_hemis") is not None:
            p.set_defaults(aparc_hemis=_coerce_list(get_cfg("aparc_hemis")) or ["lh", "rh"])
        if get_cfg("skip_sub") is not None:
            p.set_defaults(skip_sub=str(get_cfg("skip_sub")))
        if get_cfg("skip_file") is not None:
            p.set_defaults(skip_file=Path(get_cfg("skip_file")))
        if get_cfg("surf") is not None:
            p.set_defaults(surf=bool(get_cfg("surf")))
        if get_cfg("surf_target"):
            p.set_defaults(surf_target=str(get_cfg("surf_target")))
        if get_cfg("surf_measures") is not None:
            p.set_defaults(surf_measures=_coerce_list(get_cfg("surf_measures")) or ["thickness"])
        if get_cfg("surf_hemis") is not None:
            p.set_defaults(surf_hemis=_coerce_list(get_cfg("surf_hemis")) or ["lh", "rh"])
        if get_cfg("smooth") is not None:
            p.set_defaults(smooth=",".join(map(str, _coerce_int_list(get_cfg("smooth")) or [])))
        elif get_cfg("surf_fwhm") is not None:
            p.set_defaults(surf_fwhm=int(get_cfg("surf_fwhm")))
        if get_cfg("surf_outdir") is not None:
            p.set_defaults(surf_outdir=Path(get_cfg("surf_outdir")))
        # fsqc
        if get_cfg("qc") is not None:
            p.set_defaults(qc=bool(get_cfg("qc")))
        if get_cfg("qc_output") is not None:
            p.set_defaults(qc_output=Path(get_cfg("qc_output")))
        if get_cfg("qc_from") is not None:
            p.set_defaults(qc_from=str(get_cfg("qc_from")))
        if get_cfg("qc_fastsurfer") is not None:
            p.set_defaults(qc_fastsurfer=bool(get_cfg("qc_fastsurfer")))
        if get_cfg("qc_screenshots") is not None:
            p.set_defaults(qc_screenshots=bool(get_cfg("qc_screenshots")))
        if get_cfg("qc_surfaces") is not None:
            p.set_defaults(qc_surfaces=bool(get_cfg("qc_surfaces")))
        if get_cfg("qc_skullstrip") is not None:
            p.set_defaults(qc_skullstrip=bool(get_cfg("qc_skullstrip")))
        if get_cfg("qc_outlier") is not None:
            p.set_defaults(qc_outlier=bool(get_cfg("qc_outlier")))
        if get_cfg("qc_html") is not None:
            p.set_defaults(qc_html=bool(get_cfg("qc_html")))
        if get_cfg("qc_skip_existing") is not None:
            p.set_defaults(qc_skip_existing=bool(get_cfg("qc_skip_existing")))

    return p.parse_args(argv)


def prepare_output_directory(output_path: Path, force: bool = False) -> bool:
    """
    Prepare output directory, creating it if needed and handling overwrite confirmation.
    
    Args:
        output_path: The output directory path
        force: If True, skip overwrite confirmation
        
    Returns:
        True if directory is ready to use, False if user declined overwrite
    """
    # Create directory if it doesn't exist
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)
        return True
    
    # Directory exists - check if it's empty or if we should ask for confirmation
    if not force:
        try:
            # Check if directory has any files (ignore hidden files starting with .)
            files = [f for f in output_path.iterdir() if not f.name.startswith('.')]
            if files:
                print(f"[WARN] Output directory '{output_path}' is not empty and contains {len(files)} items.")
                print("Files/directories found:")
                for i, f in enumerate(files[:5]):  # Show first 5 items
                    print(f"  - {f.name}")
                if len(files) > 5:
                    print(f"  ... and {len(files) - 5} more items")
                
                while True:
                    response = input("Do you want to overwrite/continue? [y/N]: ").strip().lower()
                    if response in ('y', 'yes'):
                        print("[INFO] Continuing with existing output directory.")
                        return True
                    elif response in ('n', 'no', ''):
                        print("[INFO] Operation cancelled by user.")
                        return False
                    else:
                        print("Please enter 'y' for yes or 'n' for no.")
        except Exception as e:
            print(f"[WARN] Could not check directory contents: {e}", file=sys.stderr)
    
    return True


def read_participants(
    tsv_path: Path,
    participant_col: str,
    session_col: str,
) -> Tuple[List[str], List[Dict[str, str]], str, str]:
    if not tsv_path.exists():
        raise FileNotFoundError(f"participants.tsv not found: {tsv_path}")

    with tsv_path.open("r", newline="") as f:
        sniffer = csv.Sniffer()
        sample = f.read(2048)
        f.seek(0)
        dialect = csv.excel_tab
        if sniffer.has_header(sample):
            pass
        reader = csv.DictReader(f, dialect=dialect)
        # Normalize headers to their raw form but we will lookup case-insensitively
        fieldnames = reader.fieldnames or []
        rows = [dict(row) for row in reader]

    # Case-insensitive mapping for column names
    lower_map = {fn.lower(): fn for fn in fieldnames}
    # Allow common alternates for participant/session
    if participant_col.lower() not in lower_map:
        for alt in ("participant", "sub", "subject_id", "subject"):
            if alt in lower_map:
                participant_col = lower_map[alt]
                break
    else:
        participant_col = lower_map[participant_col.lower()]

    if session_col.lower() not in lower_map:
        for alt in ("session", "ses", "visit"):
            if alt in lower_map:
                session_col = lower_map[alt]
                break
    else:
        session_col = lower_map[session_col.lower()]

    return fieldnames, rows, participant_col, session_col


def scan_subjects_dir(subjects_dir: Path) -> List[Tuple[str, str, Optional[str]]]:
    """Return a list of (fsid, fsid_base, session_label) for each longitudinal timepoint.

    Skips base-only directories (those without a _ses-* suffix).
    """
    if not subjects_dir.exists():
        raise FileNotFoundError(f"subjects_dir not found: {subjects_dir}")
    if not subjects_dir.is_dir():
        raise NotADirectoryError(f"subjects_dir is not a directory: {subjects_dir}")

    entries: List[Tuple[str, str, Optional[str]]] = []
    for child in sorted(subjects_dir.iterdir()):
        if not child.is_dir():
            continue
        if ".long." in child.name:
            # Skip longitudinal derivative directories to avoid treating them as timepoints
            continue
        m = SUBJECT_DIR_PATTERN.match(child.name)
        if not m:
            continue
        base = m.group("base")
        ses = m.group("ses")
        if ses:  # this is a timepoint directory
            fsid = child.name
            entries.append((fsid, base, ses))
        # else: base-only directory, skip
    return entries


def session_to_tp(ses_label: Optional[str]) -> Optional[int]:
    if ses_label is None:
        return None
    m = SES_NUM_PATTERN.match(ses_label)
    if m:
        try:
            return int(m.group("num"))
        except ValueError:
            return None
    return None


def build_qdec_rows(
    timepoints: List[Tuple[str, str, Optional[str]]],
    participants_rows: List[Dict[str, str]],
    participant_col: str,
    session_col: Optional[str],
    include_columns: Optional[List[str]],
    strict: bool,
    skip_set: Optional[Set[str]] = None,
) -> Tuple[List[str], List[List[str]]]:
    # Normalize include columns
    available_cols = set(participants_rows[0].keys()) if participants_rows else set()
    cols_to_include: List[str]
    if include_columns:
        # Keep only those that exist
        cols_to_include = [c for c in include_columns if c in available_cols]
    else:
        cols_to_include = [c for c in available_cols if c not in {participant_col, session_col}]

    header = ["fsid", "fsid-base", "tp"] + cols_to_include

    def find_row(base: str, ses: Optional[str]) -> Optional[Dict[str, str]]:
        # exact match on base and session (if column exists)
        if session_col and session_col in available_cols and ses is not None:
            # prefer exact match
            for r in participants_rows:
                if r.get(participant_col) == base and r.get(session_col) == ses:
                    return r
        # fallback: match by participant only
        for r in participants_rows:
            if r.get(participant_col) == base:
                return r
        return None

    rows: List[List[str]] = []
    skipped_missing_sex: List[str] = []
    missing_tokens = {"", "na", "n/a", "nan", "null"}
    sex_col_idx: Optional[int] = None
    if "sex" in cols_to_include:
        sex_col_idx = cols_to_include.index("sex")

    for fsid, base, ses in timepoints:
        if skip_set and base in skip_set:
            continue
        r = find_row(base, ses)
        if r is None:
            if strict:
                raise ValueError(
                    f"No participants.tsv row found for subject {base} session {ses!r}"
                )
            # fill NA values when not strict
            values = ["n/a" for _ in cols_to_include]
        else:
            values = [r.get(c, "n/a") for c in cols_to_include]

        if sex_col_idx is not None:
            sex_value = values[sex_col_idx]
            norm_sex = str(sex_value).strip().lower() if sex_value is not None else None
            if norm_sex is None or norm_sex in missing_tokens:
                skipped_missing_sex.append(fsid)
                continue

        tp = session_to_tp(ses)
        tp_str = str(tp) if tp is not None else "n/a"
        rows.append([fsid, base, tp_str] + values)

    # sort by base, then numeric tp if possible
    def sort_key(row: List[str]):
        base = row[1]
        try:
            tp_val = int(row[2])
        except Exception:
            tp_val = 10**9
        return (base, tp_val, row[0])

    rows.sort(key=sort_key)
    if skipped_missing_sex:
        limit = 10
        sample = ", ".join(skipped_missing_sex[:limit])
        more = " ..." if len(skipped_missing_sex) > limit else ""
        print(
            f"Skipped {len(skipped_missing_sex)} timepoints due to missing/invalid sex values: {sample}{more}"
        )
    return header, rows


def scan_bids_subjects(bids_root: Path) -> Tuple[Set[str], Set[Tuple[str, str]]]:
    """Scan a BIDS root for participants and (participant, session) pairs.

    Returns:
      - subjects: set of participant ids like 'sub-001'
      - sessions: set of (participant, session) like ('sub-001', 'ses-01')
    """
    if not bids_root.exists():
        raise FileNotFoundError(f"BIDS root not found: {bids_root}")
    if not bids_root.is_dir():
        raise NotADirectoryError(f"BIDS root is not a directory: {bids_root}")

    subs: Set[str] = set()
    sess: Set[Tuple[str, str]] = set()
    for child in bids_root.iterdir():
        if not child.is_dir() or not child.name.startswith("sub-"):
            continue
        sub = child.name
        subs.add(sub)
        # look for ses-* under subject
        for sesdir in child.iterdir():
            if sesdir.is_dir() and sesdir.name.startswith("ses-"):
                sess.add((sub, sesdir.name))
    return subs, sess


def summarize_consistency(
    bids_root: Optional[Path],
    subjects_dir: Path,
    participants_rows: List[Dict[str, str]],
    participant_col: str,
    session_col: Optional[str],
    timepoints: List[Tuple[str, str, Optional[str]]],
) -> None:
    """Print a summary comparing participants.tsv, subjects_dir, and optional BIDS tree.

    Reports:
      - counts of subjects/timepoints found in subjects_dir
      - subjects present in participants.tsv but missing in subjects_dir
      - subjects present in subjects_dir but missing in participants.tsv
      - if BIDS is provided: subjects/sessions present in BIDS but missing elsewhere
    """
    # Participants sets
    parts_subjects: Set[str] = set(r.get(participant_col, "") for r in participants_rows if r.get(participant_col))
    parts_pairs: Set[Tuple[str, str]] = set()
    if session_col:
        for r in participants_rows:
            sub = r.get(participant_col)
            ses = r.get(session_col)
            if sub and ses:
                parts_pairs.add((sub, ses))

    # Subjects_dir sets
    sd_subjects: Set[str] = set()
    sd_pairs: Set[Tuple[str, str]] = set()
    for fsid, base, ses in timepoints:
        sd_subjects.add(base)
        if ses:
            sd_pairs.add((base, ses))

    print("=== Qdec/Subjects summary ===")
    print(f"subjects_dir: {subjects_dir}")
    print(f"participants.tsv subjects: {len(parts_subjects)}")
    print(f"subjects_dir subjects (with any timepoints): {len(sd_subjects)}")
    print(f"subjects_dir timepoints: {len(timepoints)}")

    only_in_participants = sorted(parts_subjects - sd_subjects)
    only_in_subjects_dir = sorted(sd_subjects - parts_subjects)
    if only_in_participants:
        print(f"Subjects in participants.tsv but missing in subjects_dir: {len(only_in_participants)}")
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", 20)
        print(", ".join(only_in_participants[:limit]) + (" ..." if len(only_in_participants) > limit else ""))
    if only_in_subjects_dir:
        print(f"Subjects in subjects_dir but missing in participants.tsv: {len(only_in_subjects_dir)}")
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", 20)
        print(", ".join(only_in_subjects_dir[:limit]) + (" ..." if len(only_in_subjects_dir) > limit else ""))

    if bids_root:
        bids_subjects, bids_pairs = scan_bids_subjects(bids_root)
        print(f"BIDS subjects: {len(bids_subjects)}")
        missing_in_sd = sorted(bids_subjects - sd_subjects)
        missing_in_parts = sorted(bids_subjects - parts_subjects)
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", 20)
        if missing_in_sd:
            print(f"BIDS subjects missing in subjects_dir: {len(missing_in_sd)}")
            if missing_in_sd != only_in_participants:
                print(", ".join(missing_in_sd[:limit]) + (" ..." if len(missing_in_sd) > limit else ""))
        if missing_in_parts:
            print(f"BIDS subjects missing in participants.tsv: {len(missing_in_parts)}")
            if missing_in_parts != only_in_subjects_dir:
                print(", ".join(missing_in_parts[:limit]) + (" ..." if len(missing_in_parts) > limit else ""))


def write_qdec(output_path: Path, header: List[str], rows: List[List[str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, dialect=csv.excel_tab)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def _ensure_symlink(link_path: Path, target_path: Path, dry_run: bool = True, force: bool = False) -> Tuple[bool, str]:
    """Ensure link_path is a symlink to target_path.

    Returns (changed, message)
    - changed True if a new symlink was created or updated.
    - message contains a short description of the action taken or why it was skipped.
    """
    # If link exists and is a symlink
    if link_path.is_symlink():
        try:
            current = link_path.resolve()
        except FileNotFoundError:
            current = None
        if current and current == target_path.resolve():
            return False, f"exists (correct symlink): {link_path} -> {target_path}"
        if not force:
            return False, f"exists (symlink to different target, use --link-force to update): {link_path}"
        if not dry_run:
            link_path.unlink()
            link_path.symlink_to(target_path, target_is_directory=True)
        return True, f"updated symlink: {link_path} -> {target_path}"
    # If link path exists but is not a symlink, do not touch
    if link_path.exists():
        return False, f"exists (not a symlink, skipping): {link_path}"
    # Create new
    if not dry_run:
        link_path.symlink_to(target_path, target_is_directory=True)
    return True, f"created symlink: {link_path} -> {target_path}"


def verify_and_link_long(
    subjects_dir: Path,
    timepoints: List[Tuple[str, str, Optional[str]]],
    link: bool = False,
    dry_run: bool = True,
    force: bool = False,
    require_stats: bool = True,
) -> None:
    """Verify presence of .long directories and optionally create symlinks.

    For each timepoint (fsid, base, ses):
      - expected long dir: <fsid>.long.<base>
      - if missing, and stats exist in <fsid>/stats/aseg.stats, optionally create a symlink
        <fsid>.long.<base> -> <fsid>

    Prints a short summary at the end.
    """
    created = 0
    updated = 0
    skipped = 0
    missing_stats: List[str] = []
    present = 0

    def has_any_evidence(tp_dir: Path) -> bool:
        """Return True if tp_dir shows evidence of a completed run.

        Evidence includes any of:
          - stats/aseg.stats
          - stats/<hemi>.aparc*.stats (classic or DKT mapped)
          - surf/<hemi>.thickness (surface measures exist)
        """
        stats_dir = tp_dir / "stats"
        surf_dir = tp_dir / "surf"
        if (stats_dir / "aseg.stats").exists():
            return True
        for hemi in ("lh", "rh"):
            # aparc variants
            for parc in ("aparc.DKTatlas.mapped", "aparc", "aparc.a2009s"):
                if (stats_dir / f"{hemi}.{parc}.stats").exists():
                    return True
            # surface measures
            if (surf_dir / f"{hemi}.thickness").exists():
                return True
        return False

    for fsid, base, ses in timepoints:
        if ".long." in fsid:
            skipped += 1
            print(f"skipping: {fsid} (already a .long entry)")
            continue
        tp_dir = subjects_dir / fsid
        long_dir = subjects_dir / f"{fsid}.long.{base}"
        stats_path = tp_dir / "stats" / "aseg.stats"

        if long_dir.exists() and long_dir.is_dir():
            present += 1
            continue

        # If long_dir is missing, optionally require some evidence of processing before linking
        if require_stats and not has_any_evidence(tp_dir):
            missing_stats.append(fsid)
            skipped += 1
            continue

        if link:
            changed, msg = _ensure_symlink(long_dir, tp_dir, dry_run=dry_run, force=force)
            if "created" in msg:
                created += 1
            elif "updated" in msg:
                updated += 1
            else:
                skipped += 1
            print(msg)
        else:
            note_missing = " [NO-EVIDENCE]" if require_stats and not stats_path.exists() else ""
            print(f"would link: {long_dir} -> {tp_dir} (use --link-long to create){note_missing}")
            skipped += 1

    print("=== Long symlink verification ===")
    print(f"Existing long dirs: {present}")
    print(f"Created: {created}, Updated: {updated}, Skipped: {skipped}")
    if missing_stats:
        print(f"Timepoints missing stats/aseg.stats in {subjects_dir}: {len(missing_stats)}")
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", 20)
        sample = ", ".join(sorted(missing_stats)[:limit])
        print(sample + (" ..." if len(missing_stats) > limit else ""))


def run_asegstats2table(qdec_path: Path, subjects_dir: Path) -> int:
    """Run asegstats2table with SUBJECTS_DIR pointing to subjects_dir."""

    aseg_bin = shutil.which("asegstats2table")
    if not aseg_bin:
        print(
            "asegstats2table not found in PATH. Source FreeSurfer before using --aseg.",
            file=sys.stderr,
        )
        return 4

    aseg_out = qdec_path.parent / "aseg.long.table"
    aseg_out.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["SUBJECTS_DIR"] = str(subjects_dir.resolve())

    cmd = [
        aseg_bin,
        "--qdec-long",
        str(qdec_path),
        "-t",
        str(aseg_out),
        "--skip",
    ]
    print(f"Running: {' '.join(cmd)} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})")

    try:
        result = subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        error_output = exc.stderr or exc.stdout or ""
        if "IndexError: list index out of range" in error_output or "list index out of range" in error_output:
            print(
                f"asegstats2table failed because no valid longitudinal data was found. "
                f"This likely means .long directories are missing or don't contain proper stats files. "
                f"Try using --link-long to create the required symlinks first, or check that FastSurfer/FreeSurfer "
                f"processing completed successfully for the timepoints.",
                file=sys.stderr,
            )
        else:
            print(
                f"asegstats2table failed with exit code {exc.returncode}. Command: {' '.join(cmd)}",
                file=sys.stderr,
            )
            if error_output:
                print(f"Error output: {error_output}", file=sys.stderr)
        return exc.returncode or 5

    print(f"Wrote asegstats2table output: {aseg_out}")
    return 0


def run_surf_mass_univariate(
    qdec_path: Path,
    subjects_dir: Path,
    target: str,
    measures: List[str],
    hemis: List[str],
    smooth_kernels: List[int],
    outdir: Optional[Path] = None,
    force: bool = False,
    dry_run: bool = False,
) -> int:
    """Prepare mass-univariate surface data using mris_preproc and mri_surf2surf.

        For each hemi and measure, creates:
            - <hemi>.<measure>.mgh
            - For each smoothing kernel k in smooth_kernels: <hemi>.<measure>_sm{k}.mgh
    under outdir (defaults to qdec_dir/surf).
    """
    mris_preproc_bin = shutil.which("mris_preproc")
    surf2surf_bin = shutil.which("mri_surf2surf")
    if not mris_preproc_bin or not surf2surf_bin:
        missing = [n for n, b in [("mris_preproc", mris_preproc_bin), ("mri_surf2surf", surf2surf_bin)] if not b]
        print(f"[WARN] Missing FreeSurfer binaries: {', '.join(missing)}. Skipping surface prep.", file=sys.stderr)
        return 8

    # Verify target surface template exists under subjects_dir (eg, subjects_dir/fsaverage)
    target_dir = subjects_dir / str(target)
    if not target_dir.exists() or not target_dir.is_dir():
        print(
            f"[WARN] Surface target '{target}' not found under {subjects_dir}. Expected directory: {target_dir}. Skipping surface prep.",
            file=sys.stderr,
        )
        return 0

    out_root = outdir if outdir is not None else (qdec_path.parent / "surf")
    # If forcing, we may remove existing files per pair; otherwise just ensure dir exists
    out_root.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["SUBJECTS_DIR"] = str(subjects_dir.resolve())

    # Ensure .long symlinks exist so that mris_preproc can resolve <fsid>.long.<base> paths
    # Skip auto-linking if in dry-run mode
    if not dry_run:
        try:
            tps = scan_subjects_dir(subjects_dir)
            verify_and_link_long(subjects_dir, tps, link=True, dry_run=False, force=False, require_stats=False)
        except Exception as e:
            print(f"[WARN] Failed to auto-link .long symlinks before surface prep: {e}", file=sys.stderr)
    else:
        print("[INFO] Skipping automatic .long symlink creation due to dry-run mode.")

    # Helper: filter QDEC rows for which the surf measure exists; return filtered qdec path
    def build_filtered_qdec_for(
        hemi: str, meas: str
    ) -> Tuple[Path, int, int, List[Tuple[str, str]]]:
        """Create a QDEC subset keeping only rows with existing surf files for (hemi, meas).

        Returns (qdec_filtered_path, kept_count, dropped_count, dropped_pairs[(fsid, base)]).
        If no rows are dropped, returns the original qdec_path.
        """
        kept_rows: List[List[str]] = []
        dropped = 0
        dropped_pairs: List[Tuple[str, str]] = []
        # Read QDEC (tab-separated) as generic CSV
        with qdec_path.open("r", newline="") as fh:
            reader = csv.reader(fh, dialect=csv.excel_tab)
            rows = list(reader)
        if not rows:
            return qdec_path, 0, 0, []
        header = rows[0]
        # Expect at least fsid and fsid-base
        try:
            fsid_idx = header.index("fsid")
            base_idx = header.index("fsid-base")
        except ValueError:
            # Unexpected format; fallback to original
            return qdec_path, len(rows) - 1, 0, []
        for row in rows[1:]:
            if not row or len(row) <= max(fsid_idx, base_idx):
                continue
            fsid = row[fsid_idx]
            base = row[base_idx]
            link_dir = subjects_dir / f"{fsid}.long.{base}"
            surf_file = link_dir / "surf" / f"{hemi}.{meas}"
            if surf_file.exists():
                kept_rows.append(row)
            else:
                dropped += 1
                dropped_pairs.append((fsid, base))
        # If nothing dropped, reuse original QDEC
        if dropped == 0:
            return qdec_path, len(kept_rows), 0, []
        # If everything dropped, skip gracefully by returning a path with no rows
        # but we'll detect 0 kept later and skip the computation
        filt_path = qdec_path.parent / f"qdec.{hemi}.{meas}.filtered.dat"
        with filt_path.open("w", newline="") as fh:
            writer = csv.writer(fh, dialect=csv.excel_tab)
            writer.writerow(header)
            for r in kept_rows:
                writer.writerow(r)
        print(f"[INFO] Filtered QDEC for {hemi}/{meas}: kept={len(kept_rows)}, dropped={dropped} -> {filt_path}")
        return filt_path, len(kept_rows), dropped, dropped_pairs

    # QC summary rows
    qc_rows: List[List[str]] = [["hemi", "measure", "kept", "dropped", "filtered_qdec", "missing_list"]]

    for hemi in hemis:
        for meas in measures:
            pre_path = out_root / f"{hemi}.{meas}.mgh"
            # If not forcing and file exists, we still rebuild base pre_path to reflect kept set.
            if force:
                try:
                    if pre_path.exists():
                        pre_path.unlink()
                except Exception:
                    pass
            # Build filtered QDEC (drop rows missing the required surf file)
            qdec_for_pair, kept, dropped, dropped_pairs = build_filtered_qdec_for(hemi, meas)
            if kept == 0:
                print(f"[WARN] Skipping surface prep for {hemi}/{meas}: no subjects with existing surf files.", file=sys.stderr)
                # record QC row with zero kept
                qc_rows.append([hemi, meas, str(kept), str(dropped), str(qdec_for_pair), ""])
                continue

            # Write missing list if any dropped
            missing_path = ""
            if dropped > 0:
                miss_file = out_root / f"{hemi}.{meas}.missing.tsv"
                if force and miss_file.exists():
                    try:
                        miss_file.unlink()
                    except Exception:
                        pass
                with miss_file.open("w", newline="") as fh:
                    w = csv.writer(fh, dialect=csv.excel_tab)
                    w.writerow(["fsid", "fsid-base"])  # header
                    for fsid, base in dropped_pairs:
                        w.writerow([fsid, base])
                missing_path = str(miss_file)

            # mris_preproc
            cmd1 = [
                mris_preproc_bin,
                "--qdec-long", str(qdec_for_pair),
                "--target", target,
                "--hemi", hemi,
                "--meas", meas,
                "--out", str(pre_path),
            ]
            print(f"Running: {' '.join(cmd1)} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})")
            if not dry_run:
                try:
                    subprocess.run(cmd1, check=True, env=env)
                except subprocess.CalledProcessError as exc:
                    print(f"mris_preproc failed (hemi={hemi}, meas={meas}) with code {exc.returncode}", file=sys.stderr)
                    return exc.returncode or 9
            else:
                print("[DRY-RUN] Would execute mris_preproc command above")
            # mri_surf2surf smoothing for each kernel
            for fwhm in smooth_kernels:
                sm_path = out_root / f"{hemi}.{meas}_sm{fwhm}.mgh"
                if force and sm_path.exists():
                    try:
                        sm_path.unlink()
                    except Exception:
                        pass
                cmd2 = [
                    surf2surf_bin,
                    "--hemi", hemi,
                    "--s", target,
                    "--sval", str(pre_path),
                    "--tval", str(sm_path),
                    "--fwhm-trg", str(fwhm),
                    "--cortex",
                    "--noreshape",
                ]
                print(f"Running: {' '.join(cmd2)} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})")
                if not dry_run:
                    try:
                        subprocess.run(cmd2, check=True, env=env)
                    except subprocess.CalledProcessError as exc:
                        print(f"mri_surf2surf failed (hemi={hemi}, meas={meas}, fwhm={fwhm}) with code {exc.returncode}", file=sys.stderr)
                        return exc.returncode or 10
                    print(f"Wrote: {pre_path}\nWrote: {sm_path}")
                else:
                    print("[DRY-RUN] Would execute mri_surf2surf command above")
                    print(f"[DRY-RUN] Would write: {pre_path}")
                    print(f"[DRY-RUN] Would write: {sm_path}")

            # record QC summary
            qc_rows.append([hemi, meas, str(kept), str(dropped), str(qdec_for_pair), missing_path])
    # Write QC summary TSV
    try:
        qc_path = out_root / "qc_summary.tsv"
        if force and qc_path.exists():
            try:
                qc_path.unlink()
            except Exception:
                pass
        if not dry_run:
            with qc_path.open("w", newline="") as fh:
                w = csv.writer(fh, dialect=csv.excel_tab)
                for r in qc_rows:
                    w.writerow(r)
            print(f"Wrote surface QC summary: {qc_path}")
        else:
            print(f"[DRY-RUN] Would write surface QC summary: {qc_path}")
    except Exception as e:
        print(f"[WARN] Failed to write surface QC summary: {e}", file=sys.stderr)

    return 0


def run_fsqc(
    qdec_path: Path,
    subjects_dir: Path,
    outdir: Optional[Path] = None,
    pick_from: str = "fsid",
    fastsurfer: bool = True,
    screenshots: bool = False,
    surfaces: bool = False,
    skullstrip: bool = False,
    outlier: bool = False,
    html: bool = False,
    skip_existing: bool = False,
    force: bool = False,
) -> int:
    """Run fsqc via run_fsqc CLI if available.

    Selects subjects from the QDEC table (fsid or fsid-base). If pick_from=base, we pass unique fsid-base.
    Returns 0 on success or when fsqc is unavailable.
    """
    fsqc_bin = shutil.which("run_fsqc")
    if not fsqc_bin:
        print("[WARN] fsqc not found (run_fsqc). Skipping --qc step. Install with: bash scripts/install.sh (or pip install fsqc)", file=sys.stderr)
        return 0

    out_root = outdir if outdir is not None else (qdec_path.parent / "fsqc")
    out_root.mkdir(parents=True, exist_ok=True)

    # parse QDEC and collect ids
    with qdec_path.open("r", newline="") as fh:
        reader = csv.reader(fh, dialect=csv.excel_tab)
        rows = list(reader)
    if not rows:
        print("[WARN] QDEC empty; skipping fsqc", file=sys.stderr)
        return 0
    header = rows[0]
    id_col = "fsid" if pick_from == "fsid" else "fsid-base"
    try:
        idx = header.index(id_col)
    except ValueError:
        print(f"[WARN] Column '{id_col}' not found in QDEC; skipping fsqc", file=sys.stderr)
        return 0
    values = [r[idx] for r in rows[1:] if len(r) > idx and r[idx]]
    # de-duplicate, preserve order
    seen = set()
    subjects = []
    for v in values:
        if v not in seen:
            seen.add(v)
            subjects.append(v)
    if not subjects:
        print("[WARN] No subjects found to run fsqc on; skipping", file=sys.stderr)
        return 0

    # Detect headless environment (no DISPLAY) and auto-disable surfaces to avoid OpenGL/GLFW errors
    try:
        disp = os.environ.get("DISPLAY", "").strip()
        headless = (disp == "")
    except Exception:
        headless = True
    if surfaces and headless:
        print("[INFO] No DISPLAY detected; disabling fsqc surfaces module to avoid OpenGL errors.", file=sys.stderr)
        surfaces = False

    cmd = [
        fsqc_bin,
        "--subjects_dir", str(subjects_dir),
        "--output_dir", str(out_root),
    ]
    # subject list; run_fsqc expects subject IDs (timepoints or bases)
    cmd += ["--subjects", *subjects]
    if fastsurfer:
        cmd.append("--fastsurfer")
    if screenshots:
        cmd.append("--screenshots")
        if html:
            cmd.append("--screenshots-html")
    if surfaces:
        cmd.append("--surfaces")
        if html:
            cmd.append("--surfaces-html")
    if skullstrip:
        cmd.append("--skullstrip")
        if html:
            cmd.append("--skullstrip-html")
    if outlier:
        cmd.append("--outlier")
    if skip_existing and not force:
        cmd.append("--skip-existing")

    env = os.environ.copy()
    env["SUBJECTS_DIR"] = str(subjects_dir.resolve())
    print(f"Running fsqc: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, env=env)
        print(f"Wrote fsqc outputs to: {out_root}")
    except subprocess.CalledProcessError as exc:
        print(f"[WARN] fsqc failed with exit code {exc.returncode}; continuing. Command: {' '.join(cmd)}", file=sys.stderr)
        return 0
    return 0


def run_aparcstats2table(
    qdec_path: Path,
    subjects_dir: Path,
    parc: str = "aparc",
    measures: Optional[List[str]] = None,
    hemis: Optional[List[str]] = None,
) -> int:
    """Run aparcstats2table in longitudinal mode for provided hemis and measures.

    Returns 0 on success, non-zero on first failure.
    """
    aparc_bin = shutil.which("aparcstats2table")
    if not aparc_bin:
        print(
            "aparcstats2table not found in PATH. Source FreeSurfer before using --aparc.",
            file=sys.stderr,
        )
        return 6

    if measures is None:
        measures = ["thickness", "area", "volume"]
    if hemis is None:
        hemis = ["lh", "rh"]

    out_root = qdec_path.parent / "aparc_tables"
    out_root.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["SUBJECTS_DIR"] = str(subjects_dir.resolve())

    # Preflight: auto-detect available parcellation stats. Try user-provided parc first,
    # then FastSurfer's DKT-mapped, then classic aparc, then aparc.a2009s.
    candidate_parcs = [parc]
    if parc != "aparc.DKTatlas.mapped":
        candidate_parcs.append("aparc.DKTatlas.mapped")
    if parc != "aparc":
        candidate_parcs.append("aparc")
    if parc != "aparc.a2009s":
        candidate_parcs.append("aparc.a2009s")

    chosen_parc: Optional[str] = None
    for p in candidate_parcs:
        found = False
        for hemi in hemis:
            pattern = f"**/*.long.*/stats/{hemi}.{p}.stats"
            if list(subjects_dir.glob(pattern)):
                found = True
                break
        if found:
            chosen_parc = p
            break

    if not chosen_parc:
        print(
            f"[WARN] No aparc stats files found for any of parcs {candidate_parcs} and hemis={hemis} under {subjects_dir}. Skipping aparc tables.")
        return 0
    if chosen_parc != parc:
        print(f"[INFO] Using detected parcellation '{chosen_parc}' for aparc tables (requested '{parc}').")
    parc = chosen_parc

    for hemi in hemis:
        for meas in measures:
            out_path = out_root / f"{hemi}.{parc}.{meas}.long.table"
            cmd = [
                aparc_bin,
                "--qdec-long",
                str(qdec_path),
                "--hemi",
                hemi,
                "--meas",
                meas,
                "--parc",
                parc,
                "-t",
                str(out_path),
                "--skip",
            ]
            print(f"Running: {' '.join(cmd)} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})")
            try:
                subprocess.run(cmd, check=True, env=env)
            except subprocess.CalledProcessError as exc:
                print(
                    f"aparcstats2table failed with exit code {exc.returncode}. Command: {' '.join(cmd)}",
                    file=sys.stderr,
                )
                return exc.returncode or 7
            print(f"Wrote aparcstats2table output: {out_path}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    # If no arguments provided, show help
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        parse_args(["-h"])
        return 0
    
    args = parse_args(argv)
    
    # Early dependency check
    missing_deps = check_dependencies(args)
    if missing_deps:
        print("ERROR: Missing required dependencies:", file=sys.stderr)
        for dep in missing_deps:
            print(f"  - {dep}", file=sys.stderr)
        print("\nPlease install missing dependencies and ensure FreeSurfer is properly sourced.", file=sys.stderr)
        return 1
    
    # Existence checks
    if not args.participants.exists():
        print(f"ERROR: participants.tsv not found: {args.participants}", file=sys.stderr)
        return 2
    subj_dir: Path = args.subjects_dir
    if not subj_dir.exists() or not subj_dir.is_dir():
        print(f"ERROR: subjects_dir not found or not a directory: {subj_dir}", file=sys.stderr)
        return 2

    fieldnames, participants_rows, participant_col, session_col = read_participants(
        args.participants, args.participant_column, args.session_column
    )
    if args.inspect:
        print("participants.tsv columns:")
        for fn in fieldnames:
            print(f"- {fn}")
        return 0
    timepoints = scan_subjects_dir(subj_dir)
    # Quick overview for the user: number of bases and timepoints
    bases: Set[str] = set(tp[1] for tp in timepoints)
    print(f"[INFO] Subjects overview: bases={len(bases)}, timepoints={len(timepoints)} in {subj_dir}")

    # Build skip set from CLI options
    skip_set: Set[str] = set()
    if args.skip_sub:
        for tok in str(args.skip_sub).split(","):
            tok = tok.strip()
            if tok:
                skip_set.add(tok)
    if args.skip_file and args.skip_file.exists():
        try:
            with args.skip_file.open("r") as fh:
                for line in fh:
                    tok = line.strip()
                    if tok and not tok.startswith("#"):
                        skip_set.add(tok)
        except Exception as e:
            print(f"[WARN] Failed reading --skip-file {args.skip_file}: {e}", file=sys.stderr)

    header, rows = build_qdec_rows(
        timepoints,
        participants_rows,
        participant_col,
        session_col,
        args.include_columns,
        args.strict,
        skip_set=skip_set or None,
    )
    if skip_set:
        print(f"[INFO] Skipped subjects (fsid-base) provided: {len(skip_set)}")

    # set list limit globally for summary printing
    setattr(sys.modules[__name__], "_LIST_LIMIT", max(0, int(args.list_limit)))

    # Resolve --output: treat as directory by default; support legacy file paths ending with known extensions
    out_root = args.output
    qdec_filename = "qdec.table.dat"
    out_path: Path
    try:
        is_file_like = any(str(out_root).lower().endswith(ext) for ext in (".dat", ".tsv", ".table"))
        if is_file_like:
            # Backwards-compat: user provided a file path
            out_path = out_root
            out_root = out_path.parent if out_path.parent != Path("") else Path(".")
            print(f"[INFO] --output looks like a file; will write QDEC to: {out_path}")
            # Ensure parent directory exists
            if not prepare_output_directory(out_root, args.force):
                print("ERROR: Output directory preparation cancelled by user.", file=sys.stderr)
                return 1
        else:
            # Directory semantics (preferred)
            if not prepare_output_directory(out_root, args.force):
                print("ERROR: Output directory preparation cancelled by user.", file=sys.stderr)
                return 1
            out_path = out_root / qdec_filename
            print(f"[INFO] Output root: {out_root} (QDEC: {out_path})")
    except Exception as e:
        print(f"ERROR: Failed to prepare output directory: {e}", file=sys.stderr)
        return 1
        # Fallback: treat as file under current dir
        out_path = Path(qdec_filename)
        out_root = out_path.parent

    write_qdec(out_path, header, rows)
    print(f"Wrote Qdec file: {out_path}")
    # Detect headless environment early to reflect in effective config and downstream calls
    try:
        _disp = os.environ.get("DISPLAY", "").strip()
        _headless = (_disp == "")
    except Exception:
        _headless = True
    qc_surfaces_effective = bool(getattr(args, "qc_surfaces", False))
    if bool(getattr(args, "qc", False)) and qc_surfaces_effective and _headless:
        print("[INFO] No DISPLAY detected; fsqc surfaces will be disabled.", file=sys.stderr)
        qc_surfaces_effective = False
    # Save effective configuration for transparency/reproducibility
    try:
        eff_cfg = {
            "participants": str(args.participants),
            "subjects_dir": str(subj_dir),
            "output": str(out_path),
            "participant_column": args.participant_column,
            "session_column": args.session_column,
            "include_columns": args.include_columns,
            "strict": bool(args.strict),
            "force": bool(args.force),
            "bids": str(args.bids) if args.bids else None,
            "list_limit": int(args.list_limit),
            "verify_long": bool(args.verify_long),
            "link_long": bool(args.link_long),
            "link_dry_run": bool(args.link_dry_run),
            "link_force": bool(args.link_force),
            "aseg": bool(args.aseg),
            "aparc": bool(args.aparc),
            "aparc_parc": args.aparc_parc,
            "aparc_measures": args.aparc_measures,
            "aparc_hemis": args.aparc_hemis,
            "surf": bool(args.surf),
            "surf_target": args.surf_target,
            "surf_measures": args.surf_measures,
            "surf_hemis": args.surf_hemis,
            "smooth": _coerce_int_list(getattr(args, "smooth", None)) or [int(args.surf_fwhm)] if hasattr(args, "surf_fwhm") else None,
            "surf_outdir": str(args.surf_outdir) if args.surf_outdir else None,
            "qc": bool(args.qc),
            "qc_output": str(args.qc_output) if args.qc_output else None,
            "qc_from": args.qc_from,
            "qc_fastsurfer": bool(args.qc_fastsurfer),
            "qc_screenshots": bool(args.qc_screenshots),
            "qc_surfaces": bool(qc_surfaces_effective),
            "qc_skullstrip": bool(args.qc_skullstrip),
            "qc_outlier": bool(args.qc_outlier),
            "qc_html": bool(args.qc_html),
            "qc_skip_existing": bool(args.qc_skip_existing),
            "skip_sub": args.skip_sub,
            "skip_file": str(args.skip_file) if args.skip_file else None,
            "summary": {
                "bases": len(bases),
                "timepoints": len(timepoints),
            },
        }
        cfg_out = out_root / "prep_long.effective.json"
        with cfg_out.open("w") as fh:
            json.dump(eff_cfg, fh, indent=2, sort_keys=True)
        print(f"Wrote effective config: {cfg_out}")
    except Exception as e:
        print(f"[WARN] Failed to write effective config JSON: {e}", file=sys.stderr)
    # Optional consistency summary
    summarize_consistency(args.bids, subj_dir, participants_rows, participant_col, session_col, timepoints)
    # Optional FastSurfer .long symlink verification/creation for FreeSurfer tools compatibility
    if args.verify_long or args.link_long:
        verify_and_link_long(
            subj_dir,
            timepoints,
            link=args.link_long,
            dry_run=args.link_dry_run,
            force=args.link_force,
        )

    # Optional tables
    if args.aseg:
        if args.link_dry_run:
            print("[INFO] Skipping asegstats2table due to --link-dry-run (symlinks not actually created).")
        elif shutil.which("asegstats2table") is None:
            print("[WARN] asegstats2table not found in PATH; skipping --aseg. Ensure FreeSurfer is sourced.", file=sys.stderr)
        else:
            rc = run_asegstats2table(out_path, subj_dir)
            if rc != 0:
                return rc
    if args.aparc:
        if args.link_dry_run:
            print("[INFO] Skipping aparcstats2table due to --link-dry-run (symlinks not actually created).")
        elif shutil.which("aparcstats2table") is None:
            print("[WARN] aparcstats2table not found in PATH; skipping --aparc. Ensure FreeSurfer is sourced.", file=sys.stderr)
        else:
            rc = run_aparcstats2table(
                out_path,
                subj_dir,
                parc=args.aparc_parc,
                measures=args.aparc_measures,
                hemis=args.aparc_hemis,
            )
            if rc != 0:
                return rc
    # Optional mass-univariate surface data
    if args.surf:
        have_mris = shutil.which("mris_preproc") is not None
        have_surf2 = shutil.which("mri_surf2surf") is not None
        if not (have_mris and have_surf2):
            missing = [n for n, ok in (("mris_preproc", have_mris), ("mri_surf2surf", have_surf2)) if not ok]
            print(f"[WARN] Missing FreeSurfer binaries ({', '.join(missing)}); skipping --surf.", file=sys.stderr)
        else:
            # Determine smoothing kernels list
            smooth_list = _coerce_int_list(getattr(args, "smooth", None))
            if not smooth_list:
                # fallback to single kernel from --surf-fwhm
                try:
                    smooth_list = [int(args.surf_fwhm)]
                except Exception:
                    smooth_list = [10]
            rc = run_surf_mass_univariate(
                out_path,
                subj_dir,
                target=str(args.surf_target),
                measures=list(args.surf_measures),
                hemis=list(args.surf_hemis),
                smooth_kernels=smooth_list,
                outdir=args.surf_outdir,
                force=bool(args.force),
                dry_run=bool(args.link_dry_run),
            )
            if rc != 0:
                # do not fail the entire prep if surface prep tools missing; return code already logged
                pass
    # Optional fsqc QC
    if args.qc:
        _ = run_fsqc(
            out_path,
            subj_dir,
            outdir=args.qc_output,
            pick_from=args.qc_from,
            fastsurfer=bool(args.qc_fastsurfer),
            screenshots=bool(args.qc_screenshots),
            surfaces=bool(qc_surfaces_effective),
            skullstrip=bool(args.qc_skullstrip),
            outlier=bool(args.qc_outlier),
            html=bool(args.qc_html),
            skip_existing=bool(args.qc_skip_existing),
            force=bool(args.force),
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
