#!/usr/bin/env python3
"""
Extract real acquisition dates per subject/session from BIDS T1w JSON
sidecars, and compute days elapsed since each subject's own baseline
(earliest) session.

Why this matters: nominal session labels (ses-1/ses-2/ses-3) imply fixed,
equal spacing, but actual inter-scan intervals vary substantially between
subjects (e.g. in this study, ses-1->ses-2 ranges 14-32 days, ses-2->ses-3
ranges 14-35 days). Modelling real elapsed time instead of categorical
session lets a longitudinal model estimate an actual rate of change
(volume per day) and correctly accounts for the fact that one subject's
"ses-2" may reflect a very different amount of elapsed training time than
another's.

Usage:
    python3 extract_scan_dates.py <bids_rawdata_dir> -o scan_dates.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from datetime import date
from pathlib import Path
from typing import Dict, Optional

DATE_RE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})")


def parse_date(dt_str: str) -> Optional[date]:
    m = DATE_RE.match(dt_str)
    if not m:
        return None
    return date(int(m.group(1)), int(m.group(2)), int(m.group(3)))


def find_t1w_json(anat_dir: Path) -> Optional[Path]:
    matches = sorted(anat_dir.glob("*T1w.json"))
    return matches[0] if matches else None


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("rawdata_dir", help="Path to the BIDS rawdata directory (contains sub-*/ses-*/anat/)")
    p.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = p.parse_args()

    rawdata_dir = Path(args.rawdata_dir)
    if not rawdata_dir.is_dir():
        print(f"Error: {rawdata_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    rows = []
    missing = []
    for sub_dir in sorted(rawdata_dir.glob("sub-*")):
        if not sub_dir.is_dir():
            continue
        subject_id = sub_dir.name
        for ses_dir in sorted(sub_dir.glob("ses-*")):
            session = ses_dir.name
            anat_dir = ses_dir / "anat"
            if not anat_dir.is_dir():
                missing.append((subject_id, session, "no anat/ directory"))
                continue
            json_path = find_t1w_json(anat_dir)
            if json_path is None:
                missing.append((subject_id, session, "no T1w.json found"))
                continue
            with open(json_path) as f:
                meta = json.load(f)
            dt_str = meta.get("AcquisitionDateTime")
            if not dt_str:
                missing.append((subject_id, session, "no AcquisitionDateTime in T1w.json"))
                continue
            d = parse_date(dt_str)
            if d is None:
                missing.append((subject_id, session, f"unparseable date: {dt_str}"))
                continue
            rows.append({"subject_id": subject_id, "session": session, "date": d})

    # Compute days since each subject's ses-1 scan specifically. Anchoring to
    # "whichever session is earliest for this subject" would silently
    # mislabel elapsed time for anyone missing ses-1 (their ses-2 would be
    # relabelled as day 0) -- instead, such subjects get days_since_baseline
    # left blank and are logged as missing, matching how they're already
    # excluded from baseline-referenced analyses elsewhere in this pipeline.
    by_subject: Dict[str, list] = {}
    for row in rows:
        by_subject.setdefault(row["subject_id"], []).append(row)

    out_rows = []
    for subject_id, sub_rows in by_subject.items():
        baseline_rows = [r for r in sub_rows if r["session"] == "ses-1"]
        if not baseline_rows:
            for r in sub_rows:
                out_rows.append({
                    "subject_id": subject_id,
                    "session": r["session"],
                    "acquisition_date": r["date"].isoformat(),
                    "days_since_baseline": "",
                })
            missing.append((subject_id, "ses-1", "no ses-1 scan; days_since_baseline left blank for all sessions"))
            continue
        baseline_date = baseline_rows[0]["date"]
        for r in sub_rows:
            out_rows.append({
                "subject_id": subject_id,
                "session": r["session"],
                "acquisition_date": r["date"].isoformat(),
                "days_since_baseline": (r["date"] - baseline_date).days,
            })

    out_rows.sort(key=lambda r: (r["subject_id"], r["session"]))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["subject_id", "session", "acquisition_date", "days_since_baseline"], delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {out_path} ({len(out_rows)} subject/session rows, {len(by_subject)} subjects)")

    if missing:
        missing_path = out_path.parent / "scan_dates_missing.tsv"
        with open(missing_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["subject_id", "session", "reason"])
            w.writerows(missing)
        print(f"Wrote {missing_path} ({len(missing)} missing subject/session dates)", file=sys.stderr)


if __name__ == "__main__":
    main()
