#!/usr/bin/env python3
"""
Extract longitudinal amygdala nucleus volumes (FreeSurfer 8.2.0
segment_subregions hippo-amygdala -- the same command that produces the
hippocampal subfield output already handled by extract_hippo_subfields.py;
amygdala nuclei come out of that same run for free) plus eTIV.

Emits the same two output formats as extract_hippo_subfields.py:
  - a wide "long.table" compatible with scripts/fslmer_univariate.R
  - a tidy long TSV for the mixed-effects scripts under scripts/analysis/

Input layout expected:
  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/mri/lh.amygNucVolumes.long.txt
  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/mri/rh.amygNucVolumes.long.txt
  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/stats/aseg.stats   (for eTIV)
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SUBJECT_SES_RE = re.compile(r"^(?P<base>sub-[^/_]+)_(?P<ses>ses-[^/]+)$")


def sanitize_colname(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z_]", "_", name)


def parse_amyg_file(path: Path) -> Dict[str, float]:
    volumes: Dict[str, float] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            name, value = parts
            try:
                volumes[name] = float(value)
            except ValueError:
                continue
    return volumes


def parse_etiv(stats_file: Path) -> Optional[float]:
    if not stats_file.is_file():
        return None
    with open(stats_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Measure EstimatedTotalIntraCranialVol") or line.startswith(
                "# Measure eTIV"
            ):
                parts = line.split(",")
                if len(parts) >= 4:
                    try:
                        return float(parts[3].strip())
                    except ValueError:
                        return None
    return None


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
    p.add_argument("--long-table-name", default="amygdala.long.table")
    p.add_argument("--tidy-name", default="amygdala_tidy.tsv")
    p.add_argument("--missing-log", default="missing_amygdala_subjects.tsv")
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
    all_cols: set[str] = set()

    for subject_base, session, d in timepoints:
        mri_dir = d / "mri"
        lh_file = mri_dir / "lh.amygNucVolumes.long.txt"
        rh_file = mri_dir / "rh.amygNucVolumes.long.txt"
        etiv = parse_etiv(d / "stats" / "aseg.stats")

        missing_hemis = [h for h, fp in (("lh", lh_file), ("rh", rh_file)) if not fp.is_file()]
        if missing_hemis:
            missing.append(
                {
                    "subject_id": subject_base,
                    "session": session,
                    "missing_hemispheres": ",".join(missing_hemis),
                    "reason": "amygNucVolumes.long.txt not found (segmentation incomplete or failed)",
                }
            )
            continue

        fsid = f"{subject_base}_{session}"
        fsid_long = f"{fsid}.long.{subject_base}"
        wide_row: Dict[str, object] = {
            "Measure.volume": fsid_long,
            "fsid": fsid,
            "fsid_base": subject_base,
            "eTIV": etiv,
        }

        for hemi, fpath in (("lh", lh_file), ("rh", rh_file)):
            volumes = parse_amyg_file(fpath)
            for nucleus, vol in volumes.items():
                col = f"{hemi}_{sanitize_colname(nucleus)}"
                all_cols.add(col)
                wide_row[col] = vol
                tidy_rows.append(
                    {
                        "subject_id": subject_base,
                        "session": session,
                        "fsid": fsid,
                        "hemisphere": hemi,
                        "nucleus": nucleus,
                        "is_composite": nucleus == "Whole_amygdala",
                        "volume": vol,
                        "etiv": etiv,
                    }
                )

        wide_rows[fsid_long] = wide_row

    if not wide_rows:
        print("No complete subject/session amygdala outputs found.", file=sys.stderr)
        sys.exit(1)

    long_table_path = outdir / args.long_table_name
    fixed_cols = ["Measure.volume", "fsid", "fsid_base", "eTIV"]
    data_cols = sorted(all_cols)
    header = fixed_cols + data_cols
    with open(long_table_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for fsid_long in sorted(wide_rows):
            row = wide_rows[fsid_long]
            w.writerow([row.get(c, "NA") for c in header])
    print(f"Wrote {long_table_path} ({len(wide_rows)} rows, {len(data_cols)} nucleus columns)")

    tidy_path = outdir / args.tidy_name
    tidy_cols = ["subject_id", "session", "fsid", "hemisphere", "nucleus", "is_composite", "volume", "etiv"]
    with open(tidy_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tidy_cols, delimiter="\t")
        w.writeheader()
        w.writerows(tidy_rows)
    print(f"Wrote {tidy_path} ({len(tidy_rows)} rows)")

    if missing:
        missing_path = outdir / args.missing_log
        with open(missing_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["subject_id", "session", "missing_hemispheres", "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(missing)
        print(f"Wrote {missing_path} ({len(missing)} subject/session timepoints missing amygdala output)", file=sys.stderr)
    else:
        print("No missing subject/session timepoints.", file=sys.stderr)


if __name__ == "__main__":
    main()
