#!/usr/bin/env python3
"""
Convert the study-129 FreeSurfer 8.2.0 longitudinal subregion tables pulled
from the remote DataLad repository (MRI-Lab_Repository/129/derivatives/
freesurfer/subregion_results/) into this project's standard subject-ID
convention and output shapes.

Why this exists
----------------
The remote tables use the short raw-data subject IDs (sub-003, sub-160, ...)
and one row per subject-timepoint ("sub-003_ses-1.long.sub-003"). Every other
table in this project (participants.tsv, the extract_*_subfields.py outputs,
scripts/analysis/*.R) keys on the long-form ID used throughout the local
derivatives tree instead: sub-129<G><NNN>, where G is a single digit
(1=Single, 2=Group, 3=Control) and NNN is the short numeric ID, zero-padded
to 3 digits. This script performs that one conversion, then emits both
output shapes used by the analysis stream (MRI-Lab-Graz/flex-analysis):

  - tidy long TSV   (one row per subject/session/hemisphere/region) --
    same column conventions as extract_hippo_subfields.py /
    extract_amygdala_subfields.py, for scripts/analysis/*.R there
  - wide TSV        (one row per subject-session, one column per region) --
    same shape as the *.long.table files, directly importable into SPSS

Subject ID mapping source
--------------------------
/data/local/129_PK01/rawdata/participants.tsv has columns
(participant_id, session_id, group, age, sex, size, weight) keyed on the
short ID, with `group` already the single digit (1/2/3) embedded in the
long-form ID. This is the authoritative source for the short->long mapping
(cross-checked against existing long-form derivatives directory names).

Input layout expected (as rsync'd from the remote):
  <indir>/lh.hippoSfVolumes_longitudinal_concat.csv
  <indir>/rh.hippoSfVolumes_longitudinal_concat.csv
  <indir>/lh.amygNucVolumes_longitudinal_concat.csv
  <indir>/rh.amygNucVolumes_longitudinal_concat.csv
  <indir>/ThalamicNuclei_longitudinal_concat.csv
  <indir>/brainstemSsLabels_longitudinal_concat.csv

eTIV
----
The remote concat tables don't carry eTIV themselves. Pass --etiv pointing at
a two-column TSV (timepoint, etiv) keyed on the same "sub-N_ses-N.long.sub-N"
identifier -- e.g. pulled from the matching FreeSurfer 8.2.0 aseg.stats files
on the same remote (same pipeline run as the subfields, so it's the correct
eTIV to pair with these volumes; don't mix in an eTIV from a different
FreeSurfer version). If omitted, output rows simply have no etiv column.

Two wide shapes are written per domain:
  - <name>.wide.tsv       one row per subject-SESSION (fsid/fsid_base kept
                           separate) -- matches this repo's existing
                           *.long.table convention, for scripts/fslmer_univariate.R
  - <name>.spss_wide.tsv  one row per SUBJECT, session folded into each
                           column name (e.g. lh_CA1_body_ses-1, ..._ses-2) --
                           the shape SPSS's repeated-measures GLM expects
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

TIMEPOINT_RE = re.compile(r"^(?P<short_id>sub-\d+)_(?P<ses>ses-\d+)\.long\.sub-\d+$")

# same grouping used by scripts/extract_hippo_subfields.py, kept in sync by hand
HIPPO_REGION_MAP = {
    "Hippocampal_tail": ("tail", False),
    "hippocampal-fissure": ("other", False),
    "parasubiculum": ("other", False),
    "fimbria": ("other", False),
    "HATA": ("other", False),
    "Whole_hippocampus": ("whole", True),
    "Whole_hippocampal_body": ("body", True),
    "Whole_hippocampal_head": ("head", True),
}


def hippo_region_for(subfield: str) -> Tuple[str, bool]:
    if subfield in HIPPO_REGION_MAP:
        return HIPPO_REGION_MAP[subfield]
    if subfield.endswith("-head"):
        return "head", False
    if subfield.endswith("-body"):
        return "body", False
    return "other", False


def hippo_base_name(subfield: str) -> str:
    if subfield.endswith("-head") or subfield.endswith("-body"):
        return subfield.rsplit("-", 1)[0]
    return subfield


def sanitize_colname(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z_]", "_", name)


def load_id_map(participants_tsv: Path) -> Dict[str, str]:
    """short_id (sub-003) -> long_id (sub-1291003), derived from group digit."""
    id_map: Dict[str, str] = {}
    with open(participants_tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            short_id = row["participant_id"].strip()
            group_digit = row["group"].strip()
            m = re.match(r"sub-(\d+)$", short_id)
            if not m:
                continue
            num = m.group(1).zfill(3)
            long_id = f"sub-129{group_digit}{num}"
            id_map[short_id] = long_id
    return id_map


def parse_timepoint(value: str) -> Optional[Tuple[str, str]]:
    m = TIMEPOINT_RE.match(value.strip())
    if not m:
        return None
    return m.group("short_id"), m.group("ses")


def load_concat_csv(path: Path) -> List[Dict[str, str]]:
    with open(path) as f:
        return list(csv.DictReader(f))


def load_etiv_map(path: Optional[Path]) -> Dict[str, float]:
    """timepoint (sub-N_ses-N.long.sub-N, short-ID form) -> eTIV."""
    if path is None:
        return {}
    etiv_map: Dict[str, float] = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                etiv_map[row["timepoint"].strip()] = float(row["etiv"])
            except (ValueError, KeyError):
                continue
    return etiv_map


def convert_hemi_pair(
    lh_path: Path, rh_path: Path, id_map: Dict[str, str], value_kind: str, etiv_map: Dict[str, float]
) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, object]], List[Dict[str, str]]]:
    """Shared logic for hippo (value_kind='hippo') and amygdala (value_kind='amyg')."""
    tidy_rows: List[Dict[str, object]] = []
    wide_rows: Dict[str, Dict[str, object]] = {}
    unmapped: List[Dict[str, str]] = []

    for hemi, path in (("lh", lh_path), ("rh", rh_path)):
        rows = load_concat_csv(path)
        data_cols = [c for c in rows[0].keys() if c != "timepoint"] if rows else []
        for row in rows:
            parsed = parse_timepoint(row["timepoint"])
            if parsed is None:
                unmapped.append({"timepoint": row["timepoint"], "reason": "did not match sub-N_ses-N.long.sub-N pattern"})
                continue
            short_id, session = parsed
            long_id = id_map.get(short_id)
            if long_id is None:
                unmapped.append({"timepoint": row["timepoint"], "reason": f"{short_id} not found in participants.tsv"})
                continue
            etiv = etiv_map.get(row["timepoint"].strip())

            fsid = f"{long_id}_{session}"
            fsid_long = f"{fsid}.long.{long_id}"
            wide_row = wide_rows.setdefault(
                fsid_long, {"Measure.volume": fsid_long, "fsid": fsid, "fsid_base": long_id, "eTIV": etiv}
            )

            for col in data_cols:
                try:
                    vol = float(row[col])
                except (ValueError, TypeError):
                    continue
                wide_row[f"{hemi}_{sanitize_colname(col)}"] = vol

                if value_kind == "hippo":
                    region, is_composite = hippo_region_for(col)
                    tidy_rows.append({
                        "subject_id": long_id, "session": session, "fsid": fsid, "hemisphere": hemi,
                        "subfield_raw": col, "subfield": hippo_base_name(col), "region": region,
                        "is_composite": is_composite, "volume": vol, "etiv": etiv,
                    })
                else:  # amyg
                    tidy_rows.append({
                        "subject_id": long_id, "session": session, "fsid": fsid, "hemisphere": hemi,
                        "nucleus": col, "is_composite": col == "Whole_amygdala", "volume": vol, "etiv": etiv,
                    })

    return tidy_rows, wide_rows, unmapped


def convert_thalamus(path: Path, id_map: Dict[str, str], etiv_map: Dict[str, float]) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, object]], List[Dict[str, str]]]:
    tidy_rows: List[Dict[str, object]] = []
    wide_rows: Dict[str, Dict[str, object]] = {}
    unmapped: List[Dict[str, str]] = []

    rows = load_concat_csv(path)
    data_cols = [c for c in rows[0].keys() if c != "timepoint"] if rows else []
    for row in rows:
        parsed = parse_timepoint(row["timepoint"])
        if parsed is None:
            unmapped.append({"timepoint": row["timepoint"], "reason": "did not match sub-N_ses-N.long.sub-N pattern"})
            continue
        short_id, session = parsed
        long_id = id_map.get(short_id)
        if long_id is None:
            unmapped.append({"timepoint": row["timepoint"], "reason": f"{short_id} not found in participants.tsv"})
            continue
        etiv = etiv_map.get(row["timepoint"].strip())

        fsid = f"{long_id}_{session}"
        fsid_long = f"{fsid}.long.{long_id}"
        wide_row = wide_rows.setdefault(
            fsid_long, {"Measure.volume": fsid_long, "fsid": fsid, "fsid_base": long_id, "eTIV": etiv}
        )

        for col in data_cols:
            try:
                vol = float(row[col])
            except (ValueError, TypeError):
                continue
            # columns are "Left-<nucleus>" / "Right-<nucleus>"
            if col.startswith("Left-"):
                hemi, nucleus = "lh", col[len("Left-"):]
            elif col.startswith("Right-"):
                hemi, nucleus = "rh", col[len("Right-"):]
            else:
                hemi, nucleus = "midline", col
            wide_row[f"{hemi}_{sanitize_colname(nucleus)}"] = vol
            tidy_rows.append({
                "subject_id": long_id, "session": session, "fsid": fsid, "hemisphere": hemi,
                "nucleus": nucleus, "is_composite": nucleus.startswith("Whole_"), "volume": vol, "etiv": etiv,
            })

    return tidy_rows, wide_rows, unmapped


def convert_brainstem(path: Path, id_map: Dict[str, str], etiv_map: Dict[str, float]) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, object]], List[Dict[str, str]]]:
    tidy_rows: List[Dict[str, object]] = []
    wide_rows: Dict[str, Dict[str, object]] = {}
    unmapped: List[Dict[str, str]] = []

    rows = load_concat_csv(path)
    data_cols = [c for c in rows[0].keys() if c != "timepoint"] if rows else []
    for row in rows:
        parsed = parse_timepoint(row["timepoint"])
        if parsed is None:
            unmapped.append({"timepoint": row["timepoint"], "reason": "did not match sub-N_ses-N.long.sub-N pattern"})
            continue
        short_id, session = parsed
        long_id = id_map.get(short_id)
        if long_id is None:
            unmapped.append({"timepoint": row["timepoint"], "reason": f"{short_id} not found in participants.tsv"})
            continue
        etiv = etiv_map.get(row["timepoint"].strip())

        fsid = f"{long_id}_{session}"
        fsid_long = f"{fsid}.long.{long_id}"
        wide_row = wide_rows.setdefault(
            fsid_long, {"Measure.volume": fsid_long, "fsid": fsid, "fsid_base": long_id, "eTIV": etiv}
        )

        for col in data_cols:
            try:
                vol = float(row[col])
            except (ValueError, TypeError):
                continue
            wide_row[sanitize_colname(col)] = vol
            tidy_rows.append({
                "subject_id": long_id, "session": session, "fsid": fsid, "hemisphere": "midline",
                "structure": col, "is_composite": col == "Whole_brainstem", "volume": vol, "etiv": etiv,
            })

    return tidy_rows, wide_rows, unmapped


def write_tidy(rows: List[Dict[str, object]], cols: List[str], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {path} ({len(rows)} rows)")


def write_wide(wide_rows: Dict[str, Dict[str, object]], path: Path) -> None:
    """One row per subject-SESSION -- matches this repo's *.long.table convention."""
    fixed_cols = ["Measure.volume", "fsid", "fsid_base", "eTIV"]
    data_cols = sorted({c for row in wide_rows.values() for c in row if c not in fixed_cols})
    header = fixed_cols + data_cols
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for fsid_long in sorted(wide_rows):
            row = wide_rows[fsid_long]
            w.writerow([row.get(c, "NA") for c in header])
    print(f"Wrote {path} ({len(wide_rows)} rows, {len(data_cols)} data columns)")


def pivot_to_subject_wide(wide_rows: Dict[str, Dict[str, object]]) -> Tuple[List[str], List[str], Dict[str, Dict[str, Dict[str, object]]]]:
    """Shared pivot for both the TSV and .sav SPSS-wide writers.
    Returns (sessions, data_cols, by_subject) where by_subject is
    fsid_base -> {session: row}."""
    fixed_cols = {"Measure.volume", "fsid", "fsid_base", "eTIV"}
    by_subject: Dict[str, Dict[str, Dict[str, object]]] = {}
    sessions_seen: set[str] = set()
    for row in wide_rows.values():
        subject = row["fsid_base"]
        session = row["fsid"].rsplit("_", 1)[-1]  # "sub-X_ses-N" -> "ses-N"
        by_subject.setdefault(subject, {})[session] = row
        sessions_seen.add(session)
    sessions = sorted(sessions_seen)
    data_cols = sorted({c for row in wide_rows.values() for c in row if c not in fixed_cols})
    return sessions, data_cols, by_subject


def write_spss_wide(wide_rows: Dict[str, Dict[str, object]], path: Path) -> None:
    """One row per SUBJECT, session folded into each column name -- what SPSS's
    repeated-measures GLM expects (fsid_base, then <col>_<session> per session)."""
    sessions, data_cols, by_subject = pivot_to_subject_wide(wide_rows)
    header = ["subject_id"] + [f"eTIV_{s}" for s in sessions] + [f"{c}_{s}" for s in sessions for c in data_cols]

    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for subject in sorted(by_subject):
            sess_rows = by_subject[subject]
            out_row = [subject]
            for s in sessions:
                out_row.append(sess_rows.get(s, {}).get("eTIV", "NA"))
            for s in sessions:
                row = sess_rows.get(s, {})
                for c in data_cols:
                    out_row.append(row.get(c, "NA"))
            w.writerow(out_row)
    print(f"Wrote {path} ({len(by_subject)} subjects, {len(sessions)} sessions x {len(data_cols)} data columns)")


def write_spss_sav(wide_rows: Dict[str, Dict[str, object]], path: Path) -> None:
    """Same one-row-per-subject shape as write_spss_wide, written as a native
    SPSS .sav file (via pyreadstat) instead of a TSV a colleague would have to
    import by hand."""
    import pandas as pd
    import pyreadstat

    sessions, data_cols, by_subject = pivot_to_subject_wide(wide_rows)

    def sav_name(name: str) -> str:
        # SPSS variable names: start with a letter, only letters/digits/underscore/dot,
        # no hyphens (our "ses-1" suffixes have one), max 64 bytes.
        n = re.sub(r"[^0-9A-Za-z_.]", "_", name)
        if not re.match(r"^[A-Za-z]", n):
            n = f"v_{n}"
        return n[:64]

    column_labels: Dict[str, str] = {}
    records = []
    for subject in sorted(by_subject):
        sess_rows = by_subject[subject]
        record: Dict[str, object] = {"subject_id": subject}
        for s in sessions:
            col = sav_name(f"eTIV_{s}")
            record[col] = sess_rows.get(s, {}).get("eTIV", None)
            column_labels[col] = f"Estimated total intracranial volume, {s}"
        for s in sessions:
            row = sess_rows.get(s, {})
            for c in data_cols:
                col = sav_name(f"{c}_{s}")
                record[col] = row.get(c, None)
                column_labels[col] = f"{c} ({s})"
        records.append(record)

    df = pd.DataFrame.from_records(records)
    # pyreadstat needs real NaN, not None/"NA" strings, for numeric columns
    for col in df.columns:
        if col != "subject_id":
            df[col] = pd.to_numeric(df[col], errors="coerce")

    pyreadstat.write_sav(df, str(path), column_labels=column_labels)
    print(f"Wrote {path} ({len(by_subject)} subjects, {len(sessions)} sessions x {len(data_cols)} data columns)")


def write_unmapped(unmapped: List[Dict[str, str]], path: Path) -> None:
    if not unmapped:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["timepoint", "reason"], delimiter="\t")
        w.writeheader()
        w.writerows(unmapped)
    print(f"Wrote {path} ({len(unmapped)} unmapped rows)", file=sys.stderr)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("indir", help="Directory containing the rsync'd *_longitudinal_concat.csv files")
    p.add_argument("-o", "--outdir", required=True, help="Output directory for tidy/wide TSVs")
    p.add_argument(
        "--participants",
        default="/data/local/129_PK01/rawdata/participants.tsv",
        help="Short-ID participants TSV used to derive the short->long subject ID map [default: %(default)s]",
    )
    p.add_argument(
        "--etiv",
        default=None,
        help="Optional TSV (timepoint, etiv) keyed on the short-ID timepoint string, e.g. pulled from the "
             "matching FreeSurfer run's aseg.stats. Omit to leave eTIV blank.",
    )
    args = p.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    id_map = load_id_map(Path(args.participants))
    print(f"Loaded {len(id_map)} short->long subject ID mappings from {args.participants}")

    etiv_map = load_etiv_map(Path(args.etiv) if args.etiv else None)
    print(f"Loaded {len(etiv_map)} eTIV values" + (f" from {args.etiv}" if args.etiv else " (none provided)"))

    all_unmapped: List[Dict[str, str]] = []

    # --- hippocampus ---
    tidy, wide, unmapped = convert_hemi_pair(
        indir / "lh.hippoSfVolumes_longitudinal_concat.csv",
        indir / "rh.hippoSfVolumes_longitudinal_concat.csv",
        id_map, "hippo", etiv_map,
    )
    write_tidy(tidy, ["subject_id", "session", "fsid", "hemisphere", "subfield_raw", "subfield", "region", "is_composite", "volume", "etiv"],
               outdir / "hippo_subfields_tidy_remote.tsv")
    write_wide(wide, outdir / "hippo_subfields_remote.wide.tsv")
    write_spss_wide(wide, outdir / "hippo_subfields_remote.spss_wide.tsv")
    write_spss_sav(wide, outdir / "hippo_subfields_remote.sav")
    all_unmapped += unmapped

    # --- amygdala ---
    tidy, wide, unmapped = convert_hemi_pair(
        indir / "lh.amygNucVolumes_longitudinal_concat.csv",
        indir / "rh.amygNucVolumes_longitudinal_concat.csv",
        id_map, "amyg", etiv_map,
    )
    write_tidy(tidy, ["subject_id", "session", "fsid", "hemisphere", "nucleus", "is_composite", "volume", "etiv"],
               outdir / "amygdala_tidy_remote.tsv")
    write_wide(wide, outdir / "amygdala_remote.wide.tsv")
    write_spss_wide(wide, outdir / "amygdala_remote.spss_wide.tsv")
    write_spss_sav(wide, outdir / "amygdala_remote.sav")
    all_unmapped += unmapped

    # --- thalamic nuclei ---
    tidy, wide, unmapped = convert_thalamus(indir / "ThalamicNuclei_longitudinal_concat.csv", id_map, etiv_map)
    write_tidy(tidy, ["subject_id", "session", "fsid", "hemisphere", "nucleus", "is_composite", "volume", "etiv"],
               outdir / "thalamic_nuclei_tidy_remote.tsv")
    write_wide(wide, outdir / "thalamic_nuclei_remote.wide.tsv")
    write_spss_wide(wide, outdir / "thalamic_nuclei_remote.spss_wide.tsv")
    write_spss_sav(wide, outdir / "thalamic_nuclei_remote.sav")
    all_unmapped += unmapped

    # --- brainstem substructures ---
    tidy, wide, unmapped = convert_brainstem(indir / "brainstemSsLabels_longitudinal_concat.csv", id_map, etiv_map)
    write_tidy(tidy, ["subject_id", "session", "fsid", "hemisphere", "structure", "is_composite", "volume", "etiv"],
               outdir / "brainstem_tidy_remote.tsv")
    write_wide(wide, outdir / "brainstem_remote.wide.tsv")
    write_spss_wide(wide, outdir / "brainstem_remote.spss_wide.tsv")
    write_spss_sav(wide, outdir / "brainstem_remote.sav")
    all_unmapped += unmapped

    write_unmapped(all_unmapped, outdir / "unmapped_timepoints.tsv")
    if not all_unmapped:
        print("No unmapped timepoints across any table.")


if __name__ == "__main__":
    main()
