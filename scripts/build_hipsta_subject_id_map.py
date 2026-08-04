#!/usr/bin/env python3
"""
Build a short_id <-> long_id subject mapping between the hipsta batch's
source segmentations (derivatives/freesurfer_hpc, sequential pseudonymized
IDs sub-001..sub-279) and configs/participants.tsv (original study codes,
sub-1291003 etc: site/cohort prefix "129" + a 1-digit sub-cohort digit + a
3-digit subject number).

Mapping rule (confirmed against real data, 2026-08-04): the LAST THREE DIGITS
of a participants.tsv long-form ID match the sequential short-form ID's
number, e.g. sub-1291003 <-> sub-003. Checked against the actual 149 distinct
short IDs present in derivatives/freesurfer_hpc: 147/149 map to exactly one
participants.tsv row this way. The remaining two are NOT guessed:

  - sub-072: no participants.tsv row ends in "072" -- left unmapped.
  - sub-164: TWO participants.tsv rows end in "164" (sub-1292164 AND
    sub-1293164, different sub-cohort digit) -- ambiguous, left unmapped.

Both are written to the output with an explicit note rather than silently
dropped, so a human can resolve them (e.g. from an external randomization
log) before they're used in a merge.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

LONG_ID_RE = re.compile(r"^sub-1\d{6}$")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("participants", help="Path to participants.tsv (long-form IDs, e.g. configs/participants.tsv)")
    p.add_argument("hipsta_dir", help="Path to the hipsta output directory (or the source FreeSurfer derivatives dir) -- short-form IDs are read from its sub-*_ses-*.long.sub-* / sub-*_ses-* directory names")
    p.add_argument("-o", "--output", required=True, help="Output path for the short_id -> long_id mapping TSV")
    args = p.parse_args()

    participants_path = Path(args.participants)
    if not participants_path.is_file():
        raise SystemExit(f"participants file not found: {participants_path}")
    hipsta_dir = Path(args.hipsta_dir)
    if not hipsta_dir.is_dir():
        raise SystemExit(f"not a directory: {hipsta_dir}")

    by_last3: dict[str, list[str]] = defaultdict(list)
    with open(participants_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            long_id = row["subject_id"]
            if LONG_ID_RE.match(long_id):
                by_last3[long_id[-3:]].append(long_id)

    short_ids = set()
    for d in hipsta_dir.iterdir():
        if not d.is_dir():
            continue
        m = re.match(r"^(sub-\d+)_ses-", d.name)
        if m:
            short_ids.add(m.group(1))

    rows = []
    n_mapped = n_unmapped = n_ambiguous = 0
    for short_id in sorted(short_ids):
        num = short_id.replace("sub-", "")
        last3 = num[-3:].zfill(3)
        candidates = by_last3.get(last3, [])
        if len(candidates) == 1:
            rows.append({"short_id": short_id, "long_id": candidates[0], "match_type": "unique_last3"})
            n_mapped += 1
        elif len(candidates) == 0:
            rows.append({"short_id": short_id, "long_id": "", "match_type": "unmapped_no_candidate"})
            n_unmapped += 1
        else:
            rows.append(
                {
                    "short_id": short_id,
                    "long_id": "",
                    "match_type": "unmapped_ambiguous:" + "|".join(candidates),
                }
            )
            n_ambiguous += 1

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["short_id", "long_id", "match_type"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {outpath} ({len(rows)} short IDs: {n_mapped} mapped, {n_unmapped} unmapped, {n_ambiguous} ambiguous)")
    if n_unmapped or n_ambiguous:
        print("Unresolved entries (left blank in long_id, NOT guessed):")
        for r in rows:
            if not r["long_id"]:
                print(f"  {r['short_id']}: {r['match_type']}")


if __name__ == "__main__":
    main()
