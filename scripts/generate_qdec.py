#!/usr/bin/env python3
"""
Generate a FreeSurfer longitudinal Qdec file from a BIDS participants.tsv and a
FastSurfer/FreeSurfer subjects directory.

This script scans the subjects_dir for subject base templates and longitudinal
timepoints, then merges covariates from participants.tsv to produce a Qdec file
with the required columns:

  - fsid:       the longitudinal timepoint subject id (e.g., sub-001_ses-1)
  - fsid-base:  the within-subject template id (e.g., sub-001)

Additional columns are carried over from participants.tsv (e.g., age, sex, group)
and a numeric 'tp' column is derived from the session label (ses-<number> when
available). Rows are sorted by fsid-base and tp.

References:
  - FreeSurfer Longitudinal Statistics:
    https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalStatistics

Example:
  python scripts/generate_qdec.py \
    --participants /path/to/participants.tsv \
    --subjects-dir /path/to/subjects_dir \
    --output qdec.table.dat
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set


SUBJECT_DIR_PATTERN = re.compile(r"^(?P<base>sub-[^/]+?)(?:_(?P<ses>ses-[^/]+))?$")
SES_NUM_PATTERN = re.compile(r"^ses-(?P<num>\d+)$")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate FreeSurfer longitudinal Qdec file")
    p.add_argument("--participants", required=True, type=Path, help="Path to BIDS participants.tsv")
    p.add_argument("--subjects-dir", required=True, type=Path, help="Path to FastSurfer/FreeSurfer subjects directory")
    p.add_argument("--output", type=Path, default=Path("qdec.table.dat"), help="Output Qdec TSV file (default: qdec.table.dat)")
    p.add_argument("--participant-column", default="participant_id", help="Column name for participant id (default: participant_id)")
    p.add_argument("--session-column", default="session_id", help="Column name for session id if present (default: session_id)")
    p.add_argument(
        "--include-columns",
        nargs="*",
        default=None,
        help="Optional explicit list of covariate columns to include from participants.tsv. "
             "If omitted, include all columns except participant and session columns.",
    )
    p.add_argument("--strict", action="store_true", help="Fail if a subjects_dir timepoint has no matching participants row")
    p.add_argument("--inspect", action="store_true", help="Print participants.tsv columns and exit")
    p.add_argument("--bids", type=Path, default=None, help="Optional BIDS root to cross-check subjects/sessions consistency")
    return p.parse_args(argv)


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
    for fsid, base, ses in timepoints:
        r = find_row(base, ses)
        if r is None:
            if strict:
                raise ValueError(
                    f"No participants.tsv row found for subject {base} session {ses!r}"
                )
            # fill NA values
            values = ["n/a" for _ in cols_to_include]
        else:
            values = [r.get(c, "n/a") for c in cols_to_include]

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
        print(", ".join(only_in_participants[:20]) + (" ..." if len(only_in_participants) > 20 else ""))
    if only_in_subjects_dir:
        print(f"Subjects in subjects_dir but missing in participants.tsv: {len(only_in_subjects_dir)}")
        print(", ".join(only_in_subjects_dir[:20]) + (" ..." if len(only_in_subjects_dir) > 20 else ""))

    if bids_root:
        bids_subjects, bids_pairs = scan_bids_subjects(bids_root)
        print(f"BIDS subjects: {len(bids_subjects)}")
        missing_in_sd = sorted(bids_subjects - sd_subjects)
        missing_in_parts = sorted(bids_subjects - parts_subjects)
        if missing_in_sd:
            print(f"BIDS subjects missing in subjects_dir: {len(missing_in_sd)}")
            print(", ".join(missing_in_sd[:20]) + (" ..." if len(missing_in_sd) > 20 else ""))
        if missing_in_parts:
            print(f"BIDS subjects missing in participants.tsv: {len(missing_in_parts)}")
            print(", ".join(missing_in_parts[:20]) + (" ..." if len(missing_in_parts) > 20 else ""))



def write_qdec(output_path: Path, header: List[str], rows: List[List[str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, dialect=csv.excel_tab)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    # Existence checks
    if not args.participants.exists():
        print(f"ERROR: participants.tsv not found: {args.participants}", file=sys.stderr)
        return 2
    if not args.subjects_dir.exists() or not args.subjects_dir.is_dir():
        print(f"ERROR: subjects_dir not found or not a directory: {args.subjects_dir}", file=sys.stderr)
        return 2

    fieldnames, participants_rows, participant_col, session_col = read_participants(
        args.participants, args.participant_column, args.session_column
    )
    if args.inspect:
        print("participants.tsv columns:")
        for fn in fieldnames:
            print(f"- {fn}")
        return 0
    timepoints = scan_subjects_dir(args.subjects_dir)
    header, rows = build_qdec_rows(
        timepoints,
        participants_rows,
        participant_col,
        session_col,
        args.include_columns,
        args.strict,
    )
    write_qdec(args.output, header, rows)
    print(f"Wrote Qdec file: {args.output}")
    # Optional consistency summary
    summarize_consistency(args.bids, args.subjects_dir, participants_rows, participant_col, session_col, timepoints)
    return 0


if __name__ == "__main__":
    sys.exit(main())
