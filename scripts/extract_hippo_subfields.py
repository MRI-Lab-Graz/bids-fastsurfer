#!/usr/bin/env python3
"""
Extract longitudinal hippocampal subfield volumes (FreeSurfer 8.2.0
segment_subregions hippo-amygdala) plus eTIV, and emit two output formats:

  - a wide "long.table" compatible with scripts/fslmer_univariate.R
    (row id = sub-X_ses-Y.long.sub-X, one column per hemisphere/subfield)
  - a tidy long TSV (one row per subject/session/hemisphere/subfield)
    for the mixed-effects scripts under scripts/analysis/

Input layout expected (produced by scripts/run_subfield_segmentation.sh /
the hippoSF batch driver):

  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/mri/lh.hippoSfVolumes.long.txt
  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/mri/rh.hippoSfVolumes.long.txt
  <SUBJECTS_DIR>/sub-<ID>_ses-<N>/stats/aseg.stats   (for eTIV)

Subjects/sessions missing either hemisphere file are logged to
--missing-log rather than silently dropped.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SUBJECT_SES_RE = re.compile(r"^(?P<base>sub-[^/_]+)_(?P<ses>ses-[^/]+)$")

# Anatomical grouping for the FreeSurfer 8.x hippo-amygdala subfield atlas.
# "region" is the head/body/tail grouping used for hierarchical/nested models;
# composite measures (already summed by FreeSurfer) are marked as such so
# they can be excluded from subfield-level models by default.
REGION_MAP = {
    "Hippocampal_tail": ("tail", False),
    "hippocampal-fissure": ("other", False),  # CSF, not parenchyma
    "parasubiculum": ("other", False),
    "fimbria": ("other", False),
    "HATA": ("other", False),
    "Whole_hippocampus": ("whole", True),
    "Whole_hippocampal_body": ("body", True),
    "Whole_hippocampal_head": ("head", True),
}


def region_for(subfield: str) -> Tuple[str, bool]:
    """Return (region, is_composite) for a raw subfield label."""
    if subfield in REGION_MAP:
        return REGION_MAP[subfield]
    if subfield.endswith("-head"):
        return "head", False
    if subfield.endswith("-body"):
        return "body", False
    return "other", False


def base_subfield_name(subfield: str) -> str:
    """Subfield identity independent of head/body (e.g. CA1-head -> CA1)."""
    if subfield.endswith("-head") or subfield.endswith("-body"):
        return subfield.rsplit("-", 1)[0]
    return subfield


def sanitize_colname(name: str) -> str:
    """Make an R-safe (and shell-safe) column name."""
    return re.sub(r"[^0-9A-Za-z_]", "_", name)


def parse_hippo_sf_file(path: Path) -> Dict[str, float]:
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
    """Return (subject_base, session, dir) for every plain sub-X_ses-Y directory."""
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
    p.add_argument(
        "--long-table-name",
        default="hippo_subfields.long.table",
        help="Filename for the fslmer_univariate.R-compatible wide table [default: %(default)s]",
    )
    p.add_argument(
        "--tidy-name",
        default="hippo_subfields_tidy.tsv",
        help="Filename for the tidy long-format TSV [default: %(default)s]",
    )
    p.add_argument(
        "--missing-log",
        default="missing_hippo_subjects.tsv",
        help="Filename (in --outdir) logging subjects/sessions missing subfield output [default: %(default)s]",
    )
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
    all_subfield_cols: set[str] = set()

    for subject_base, session, d in timepoints:
        mri_dir = d / "mri"
        lh_file = mri_dir / "lh.hippoSfVolumes.long.txt"
        rh_file = mri_dir / "rh.hippoSfVolumes.long.txt"
        etiv = parse_etiv(d / "stats" / "aseg.stats")

        missing_hemis = [h for h, fp in (("lh", lh_file), ("rh", rh_file)) if not fp.is_file()]
        if missing_hemis:
            missing.append(
                {
                    "subject_id": subject_base,
                    "session": session,
                    "missing_hemispheres": ",".join(missing_hemis),
                    "reason": "hippoSfVolumes.long.txt not found (segmentation incomplete or failed)",
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
            volumes = parse_hippo_sf_file(fpath)
            for subfield, vol in volumes.items():
                region, is_composite = region_for(subfield)
                base_name = base_subfield_name(subfield)
                col = f"{hemi}_{sanitize_colname(subfield)}"
                all_subfield_cols.add(col)
                wide_row[col] = vol

                tidy_rows.append(
                    {
                        "subject_id": subject_base,
                        "session": session,
                        "fsid": fsid,
                        "hemisphere": hemi,
                        "subfield_raw": subfield,
                        "subfield": base_name,
                        "region": region,
                        "is_composite": is_composite,
                        "volume": vol,
                        "etiv": etiv,
                    }
                )

        wide_rows[fsid_long] = wide_row

    if not wide_rows:
        print("No complete subject/session subfield outputs found.", file=sys.stderr)
        sys.exit(1)

    # --- write wide long.table (fslmer_univariate.R compatible) ---
    long_table_path = outdir / args.long_table_name
    fixed_cols = ["Measure.volume", "fsid", "fsid_base", "eTIV"]
    subfield_cols = sorted(all_subfield_cols)
    header = fixed_cols + subfield_cols
    with open(long_table_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for fsid_long in sorted(wide_rows):
            row = wide_rows[fsid_long]
            w.writerow([row.get(c, "NA") for c in header])
    print(f"Wrote {long_table_path} ({len(wide_rows)} rows, {len(subfield_cols)} subfield columns)")

    # --- write tidy long TSV ---
    tidy_path = outdir / args.tidy_name
    tidy_cols = [
        "subject_id",
        "session",
        "fsid",
        "hemisphere",
        "subfield_raw",
        "subfield",
        "region",
        "is_composite",
        "volume",
        "etiv",
    ]
    with open(tidy_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tidy_cols, delimiter="\t")
        w.writeheader()
        w.writerows(tidy_rows)
    print(f"Wrote {tidy_path} ({len(tidy_rows)} rows)")

    # --- write missing log ---
    if missing:
        missing_path = outdir / args.missing_log
        with open(missing_path, "w", newline="") as f:
            w = csv.DictWriter(
                f, fieldnames=["subject_id", "session", "missing_hemispheres", "reason"], delimiter="\t"
            )
            w.writeheader()
            w.writerows(missing)
        print(
            f"Wrote {missing_path} ({len(missing)} subject/session timepoints missing subfield output)",
            file=sys.stderr,
        )
    else:
        print("No missing subject/session timepoints.", file=sys.stderr)


if __name__ == "__main__":
    main()
