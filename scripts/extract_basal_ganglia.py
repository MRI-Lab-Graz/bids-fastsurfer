#!/usr/bin/env python3
"""
Extract longitudinal basal ganglia volumes (Caudate, Putamen, Pallidum,
Accumbens-area, both hemispheres) plus eTIV from standard FreeSurfer/
FastSurfer aseg.stats. Unlike the hippocampal subfield / amygdala nucleus
outputs, this comes from the BASE recon pipeline, not the specialized
segment_subregions run -- so coverage here is typically much more complete.

Emits the same two output formats as extract_hippo_subfields.py /
extract_amygdala_subfields.py:
  - a wide "long.table" compatible with scripts/fslmer_univariate.R
  - a tidy long TSV for the mixed-effects scripts under scripts/analysis/

Input layout expected:
  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/stats/aseg.stats
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SUBJECT_SES_RE = re.compile(r"^(?P<base>sub-[^/_]+)_(?P<ses>ses-[^/]+)$")

STRUCTURES = {
    "Left-Caudate": ("lh", "Caudate"),
    "Left-Putamen": ("lh", "Putamen"),
    "Left-Pallidum": ("lh", "Pallidum"),
    "Left-Accumbens-area": ("lh", "Accumbens"),
    "Right-Caudate": ("rh", "Caudate"),
    "Right-Putamen": ("rh", "Putamen"),
    "Right-Pallidum": ("rh", "Pallidum"),
    "Right-Accumbens-area": ("rh", "Accumbens"),
}


def parse_aseg_stats(stats_file: Path) -> Tuple[Dict[str, float], Optional[float]]:
    volumes: Dict[str, float] = {}
    etiv: Optional[float] = None
    with open(stats_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Measure EstimatedTotalIntraCranialVol") or line.startswith("# Measure eTIV"):
                parts = line.split(",")
                if len(parts) >= 4:
                    try:
                        etiv = float(parts[3].strip())
                    except ValueError:
                        pass
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            struct_name = parts[4]
            if struct_name in STRUCTURES:
                try:
                    volumes[struct_name] = float(parts[3])
                except ValueError:
                    continue
    return volumes, etiv


def find_timepoints(subjects_dir: Path) -> List[Tuple[str, str, Path]]:
    out = []
    for d in sorted(subjects_dir.iterdir()):
        if not d.is_dir():
            continue
        if ".long." in d.name:
            continue
        m = SUBJECT_SES_RE.match(d.name)
        if not m:
            continue
        out.append((m.group("base"), m.group("ses"), d))
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("subjects_dir", help="Path to the FastSurfer/FreeSurfer SUBJECTS_DIR")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--long-table-name", default="basal_ganglia.long.table")
    p.add_argument("--tidy-name", default="basal_ganglia_tidy.tsv")
    p.add_argument("--missing-log", default="missing_basal_ganglia_subjects.tsv")
    args = p.parse_args()

    subjects_dir = Path(args.subjects_dir)
    if not subjects_dir.is_dir():
        print(f"Error: {subjects_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    timepoints = find_timepoints(subjects_dir)
    if not timepoints:
        print(f"No sub-*_ses-* directories found under {subjects_dir}", file=sys.stderr)
        sys.exit(1)

    tidy_rows: List[Dict[str, object]] = []
    wide_rows: Dict[str, Dict[str, object]] = {}
    missing: List[Dict[str, str]] = []

    for subject_base, session, d in timepoints:
        stats_file = d / "stats" / "aseg.stats"
        if not stats_file.is_file():
            missing.append({
                "subject_id": subject_base, "session": session,
                "reason": "aseg.stats not found",
            })
            continue

        volumes, etiv = parse_aseg_stats(stats_file)
        if len(volumes) < len(STRUCTURES):
            found = set(volumes.keys())
            missing.append({
                "subject_id": subject_base, "session": session,
                "reason": f"missing structures: {', '.join(sorted(set(STRUCTURES) - found))}",
            })
            if not volumes:
                continue

        fsid = f"{subject_base}_{session}"
        fsid_long = f"{fsid}.long.{subject_base}"
        wide_row: Dict[str, object] = {
            "Measure.volume": fsid_long,
            "fsid": fsid,
            "fsid_base": subject_base,
            "eTIV": etiv,
        }

        for struct_name, vol in volumes.items():
            hemi, roi = STRUCTURES[struct_name]
            col = f"{hemi}_{roi}"
            wide_row[col] = vol
            tidy_rows.append({
                "subject_id": subject_base,
                "session": session,
                "fsid": fsid,
                "hemisphere": hemi,
                "roi": roi,
                "volume": vol,
                "etiv": etiv,
            })

        wide_rows[fsid_long] = wide_row

    if not wide_rows:
        print("No complete subject/session basal ganglia outputs found.", file=sys.stderr)
        sys.exit(1)

    long_table_path = outdir / args.long_table_name
    fixed_cols = ["Measure.volume", "fsid", "fsid_base", "eTIV"]
    data_cols = sorted({f"{hemi}_{roi}" for hemi, roi in STRUCTURES.values()})
    header = fixed_cols + data_cols
    with open(long_table_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for fsid_long in sorted(wide_rows):
            row = wide_rows[fsid_long]
            w.writerow([row.get(c, "NA") for c in header])
    print(f"Wrote {long_table_path} ({len(wide_rows)} rows, {len(data_cols)} ROI columns)")

    tidy_path = outdir / args.tidy_name
    tidy_cols = ["subject_id", "session", "fsid", "hemisphere", "roi", "volume", "etiv"]
    with open(tidy_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tidy_cols, delimiter="\t")
        w.writeheader()
        w.writerows(tidy_rows)
    print(f"Wrote {tidy_path} ({len(tidy_rows)} rows)")

    if missing:
        missing_path = outdir / args.missing_log
        with open(missing_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["subject_id", "session", "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(missing)
        print(f"Wrote {missing_path} ({len(missing)} subject/session timepoints missing basal ganglia output)", file=sys.stderr)
    else:
        print("No missing subject/session timepoints.", file=sys.stderr)


if __name__ == "__main__":
    main()
