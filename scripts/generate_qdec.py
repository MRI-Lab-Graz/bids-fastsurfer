#!/usr/bin/env python3
"""
Generate a FreeSurfer Qdec file from a BIDS participants.tsv and a
FastSurfer/FreeSurfer subjects directory.

Supports both cross-sectional and longitudinal studies:
  - Cross-sectional: Single timepoint per subject (fsid = fsid-base)
  - Longitudinal: Multiple timepoints per subject (fsid contains session, fsid-base is subject)

This script scans the subjects_dir for subjects and timepoints, then merges
covariates from participants.tsv to produce a Qdec file with the required columns:

  - fsid:       the subject/timepoint id (e.g., sub-001 or sub-001_ses-1)
  - fsid-base:  the within-subject template id (e.g., sub-001)
  - tp:         numeric timepoint (derived from ses-<number> for longitudinal)

Additional columns are carried over from participants.tsv (e.g., age, sex, group).
Rows are sorted by fsid-base and tp.

Note: For advanced analysis (aseg/aparc tables, surface preprocessing, QC),
      use analyse_qdec.py after generating the Qdec file.

References:
  - FreeSurfer Longitudinal Statistics:
    https://surfer.nmr.mgh.harvard.edu/fswiki/LongitudinalStatistics

Examples:
  # Longitudinal study
  python scripts/generate_qdec.py \
    --participants /path/to/participants.tsv \
    --subjects-dir /path/to/subjects_dir \
    --output qdec.table.dat \
    --verify-long --link-long

  # Cross-sectional study
  python scripts/generate_qdec.py \
    --participants /path/to/participants.tsv \
    --subjects-dir /path/to/subjects_dir \
    --output qdec.table.dat

  # Verbose output
  python scripts/generate_qdec.py \
    --participants /path/to/participants.tsv \
    --subjects-dir /path/to/subjects_dir \
    --verbose
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

# Constants
SUBJECT_DIR_PATTERN = re.compile(r"^(?P<base>sub-[^/]+?)(?:_(?P<ses>ses-[^/]+))?$")
SES_NUM_PATTERN = re.compile(r"^ses-(?P<num>\d+)$")
MISSING_TOKENS = {"", "na", "n/a", "nan", "null"}
DEFAULT_LIST_LIMIT = 20
DEFAULT_OUTPUT_FILENAME = "qdec.table.dat"

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate FreeSurfer Qdec file (cross-sectional or longitudinal)"
    )
    p.add_argument("--participants", required=True, type=Path, help="Path to BIDS participants.tsv")
    p.add_argument(
        "--subjects-dir",
        required=True,
        type=Path,
        help="Path to FastSurfer/FreeSurfer subjects directory",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path(DEFAULT_OUTPUT_FILENAME),
        help=(
            f"Output Qdec TSV file (default: {DEFAULT_OUTPUT_FILENAME}). "
            "If a directory path is provided, the file will be created inside it."
        ),
    )
    p.add_argument(
        "--participant-column",
        default="participant_id",
        help="Column name for participant id (default: participant_id)",
    )
    p.add_argument(
        "--session-column",
        default="session_id",
        help="Column name for session id if present (default: session_id)",
    )
    p.add_argument(
        "--include-columns",
        nargs="*",
        default=None,
        help="Optional explicit list of covariate columns to include from participants.tsv. "
        "If omitted, include all columns except participant and session columns.",
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="Fail if a subjects_dir timepoint has no matching participants row",
    )
    p.add_argument("--inspect", action="store_true", help="Print participants.tsv columns and exit")
    p.add_argument(
        "--bids",
        type=Path,
        default=None,
        help="Optional BIDS root to cross-check subjects/sessions consistency",
    )
    p.add_argument(
        "--list-limit",
        type=int,
        default=DEFAULT_LIST_LIMIT,
        help=f"Max number of IDs to show when listing missing subjects (default: {DEFAULT_LIST_LIMIT})",
    )
    # FastSurfer compatibility with FreeSurfer .long directories
    p.add_argument(
        "--verify-long",
        action="store_true",
        help="Verify presence of <fsid>.long.<base>/stats/aseg.stats for each timepoint",
    )
    p.add_argument(
        "--link-long",
        action="store_true",
        help="Create <fsid>.long.<base> symlinks pointing to the timepoint directories when missing. "
        "If not specified, automatically enabled when multiple sessions per subject are detected.",
    )
    p.add_argument(
        "--no-link-long",
        action="store_true",
        help="Disable automatic linking even when multiple sessions are detected",
    )
    p.add_argument(
        "--link-dry-run",
        action="store_true",
        help="Print the symlink actions without making changes",
    )
    p.add_argument(
        "--link-force",
        action="store_true",
        help="If an existing symlink points elsewhere, replace it (does not delete real directories)",
    )
    p.add_argument(
        "--skip-sub", default=None, help="Comma-separated fsid_base IDs to skip (exclude) from QDEC"
    )
    p.add_argument(
        "--skip-file",
        type=Path,
        default=None,
        help="File with fsid_base IDs (one per line) to skip from QDEC",
    )
    p.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
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


def has_multiple_sessions(timepoints: List[Tuple[str, str, Optional[str]]]) -> bool:
    """Check if any subject has multiple sessions (timepoints)."""
    subject_sessions = {}
    for fsid, base, ses in timepoints:
        if ses:  # Only count actual sessions
            if base not in subject_sessions:
                subject_sessions[base] = []
            subject_sessions[base].append(ses)

    # Check if any subject has more than one session
    return any(len(sessions) > 1 for sessions in subject_sessions.values())


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
            if norm_sex is None or norm_sex in MISSING_TOKENS:
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
        logger.warning(
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
    parts_subjects: Set[str] = set(
        r.get(participant_col, "") for r in participants_rows if r.get(participant_col)
    )
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

    logger.info("=== Qdec/Subjects summary ===")
    logger.info(f"subjects_dir: {subjects_dir}")
    logger.info(f"participants.tsv subjects: {len(parts_subjects)}")
    logger.info(f"subjects_dir subjects (with any timepoints): {len(sd_subjects)}")
    logger.info(f"subjects_dir timepoints: {len(timepoints)}")

    only_in_participants = sorted(parts_subjects - sd_subjects)
    only_in_subjects_dir = sorted(sd_subjects - parts_subjects)
    if only_in_participants:
        logger.info(
            f"Subjects in participants.tsv but missing in subjects_dir: {len(only_in_participants)}"
        )
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", DEFAULT_LIST_LIMIT)
        logger.info(
            ", ".join(only_in_participants[:limit])
            + (" ..." if len(only_in_participants) > limit else "")
        )
    if only_in_subjects_dir:
        logger.info(
            f"Subjects in subjects_dir but missing in participants.tsv: {len(only_in_subjects_dir)}"
        )
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", DEFAULT_LIST_LIMIT)
        logger.info(
            ", ".join(only_in_subjects_dir[:limit])
            + (" ..." if len(only_in_subjects_dir) > limit else "")
        )

    if bids_root:
        bids_subjects, bids_pairs = scan_bids_subjects(bids_root)
        logger.info(f"BIDS subjects: {len(bids_subjects)}")
        missing_in_sd = sorted(bids_subjects - sd_subjects)
        missing_in_parts = sorted(bids_subjects - parts_subjects)
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", DEFAULT_LIST_LIMIT)
        # Avoid repeating the exact same list twice. If participants-missing and BIDS-missing-in-sd are identical,
        # suppress the second detailed list and only print counts.
        if missing_in_sd:
            logger.info(f"BIDS subjects missing in subjects_dir: {len(missing_in_sd)}")
            if missing_in_sd != only_in_participants:
                logger.info(
                    ", ".join(missing_in_sd[:limit])
                    + (" ..." if len(missing_in_sd) > limit else "")
                )
        if missing_in_parts:
            logger.info(f"BIDS subjects missing in participants.tsv: {len(missing_in_parts)}")
            # This is generally different; still avoid printing if equal to only_in_subjects_dir
            if missing_in_parts != only_in_subjects_dir:
                logger.info(
                    ", ".join(missing_in_parts[:limit])
                    + (" ..." if len(missing_in_parts) > limit else "")
                )


def write_qdec(output_path: Path, header: List[str], rows: List[List[str]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f, dialect=csv.excel_tab)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def _ensure_symlink(
    link_path: Path, target_path: Path, dry_run: bool = True, force: bool = False
) -> Tuple[bool, str]:
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
            return (
                False,
                f"exists (symlink to different target, use --link-force to update): {link_path}",
            )
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

    for fsid, base, ses in timepoints:
        if ".long." in fsid:
            skipped += 1
            logger.debug(f"skipping: {fsid} (already a .long entry)")
            continue
        tp_dir = subjects_dir / fsid
        long_dir = subjects_dir / f"{fsid}.long.{base}"
        stats_path = tp_dir / "stats" / "aseg.stats"

        if long_dir.exists() and long_dir.is_dir():
            present += 1
            continue

        # If long_dir is missing, check whether stats file exists in tp_dir
        if not stats_path.exists():
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
            logger.info(msg)
        else:
            logger.info(
                f"would link: {long_dir} -> {tp_dir} (use --link-long to create){' [MISSING]' if not stats_path.exists() else ''}"
            )
            skipped += 1

    logger.info("=== Long symlink verification ===")
    logger.info(f"Existing long dirs: {present}")
    logger.info(f"Created: {created}, Updated: {updated}, Skipped: {skipped}")
    if missing_stats:
        logger.warning(
            f"Timepoints missing stats/aseg.stats in {subjects_dir}: {len(missing_stats)}"
        )
        limit = getattr(sys.modules[__name__], "_LIST_LIMIT", DEFAULT_LIST_LIMIT)
        sample = ", ".join(sorted(missing_stats)[:limit])
        logger.warning(sample + (" ..." if len(missing_stats) > limit else ""))


def build_skip_set(args) -> Set[str]:
    """Build the set of subjects to skip from CLI arguments."""
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
            logger.warning(f"Failed reading --skip-file {args.skip_file}: {e}")
    return skip_set


def resolve_output_path(args) -> Path:
    """Resolve the output path, handling directory inputs intelligently.

    If the path looks like a file (has extension like .tsv, .dat, .txt), treat as file.
    If the path looks like a directory (no extension, ends with /, or is existing dir),
    create directory and put default file inside it.
    If path exists as a file but no extension, append .dat extension.
    """
    out_path = args.output

    # Check if it looks like a file (has extension)
    path_str = str(out_path)
    if "." in path_str and not path_str.endswith("."):
        # Has extension, treat as file
        return out_path

    # No extension - check if it exists
    if out_path.exists():
        if out_path.is_file():
            # Exists as file, append .dat extension
            out_path = out_path.with_suffix(".dat")
            logger.info(f"--output exists as file; writing to: {out_path}")
            return out_path
        elif out_path.is_dir():
            # Exists as directory, add default filename
            out_path = out_path / DEFAULT_OUTPUT_FILENAME
            logger.info(f"--output is existing directory; writing file to: {out_path}")
            return out_path

    # Doesn't exist or looks like directory path
    # Create directory and put default file inside
    out_path = out_path / DEFAULT_OUTPUT_FILENAME
    logger.info(f"--output looks like a directory; writing file to: {out_path}")

    return out_path


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    # Configure logging level based on verbose flag
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    # Existence checks
    if not args.participants.exists():
        logger.error(f"participants.tsv not found: {args.participants}")
        return 2
    if not args.subjects_dir.exists() or not args.subjects_dir.is_dir():
        logger.error(f"subjects_dir not found or not a directory: {args.subjects_dir}")
        return 2

    fieldnames, participants_rows, participant_col, session_col = read_participants(
        args.participants, args.participant_column, args.session_column
    )

    if args.inspect:
        logger.info("participants.tsv columns:")
        for fn in fieldnames:
            logger.info(f"- {fn}")

        # Analyze values in each column
        logger.info("\nparticipants.tsv value analysis:")
        for col in fieldnames:
            all_values = [row.get(col, "") for row in participants_rows]
            non_empty_values = [v for v in all_values if v.strip()]
            missing_count = len(all_values) - len(non_empty_values)

            if not non_empty_values:
                logger.info(f"- {col}: (empty column)")
                continue

            # Try to determine if numeric
            numeric_values = []
            for v in non_empty_values:
                try:
                    numeric_values.append(float(v))
                except (ValueError, TypeError):
                    pass

            unique_vals = set(non_empty_values)

            if len(numeric_values) > len(non_empty_values) * 0.8:  # Mostly numeric
                min_val = min(numeric_values)
                max_val = max(numeric_values)
                unique_numeric = sorted(set(numeric_values))
                if len(unique_numeric) <= 10:
                    logger.info(
                        f"- {col}: numeric, values {unique_numeric}, range [{min_val:.2f}, {max_val:.2f}]"
                    )
                else:
                    logger.info(
                        f"- {col}: numeric, range [{min_val:.2f}, {max_val:.2f}], {len(unique_vals)} unique values"
                    )
            else:
                # Categorical
                if len(unique_vals) <= 10:
                    sorted_unique = sorted(unique_vals)
                    logger.info(
                        f"- {col}: categorical, {len(unique_vals)} unique values: {', '.join(sorted_unique)}"
                    )
                else:
                    logger.info(f"- {col}: categorical, {len(unique_vals)} unique values")
                    # Show most common values
                    counts = Counter(non_empty_values)
                    most_common = counts.most_common(5)
                    logger.info(
                        f"    top values: {', '.join(f'{val}({count})' for val, count in most_common)}"
                    )

            if missing_count > 0:
                logger.info(f"    missing: {missing_count} rows")

        logger.info(f"\nTotal rows: {len(participants_rows)}")
        return 0

    timepoints = scan_subjects_dir(args.subjects_dir)
    skip_set = build_skip_set(args)

    # Determine if linking should be enabled
    multiple_sessions = has_multiple_sessions(timepoints)
    enable_linking = False

    if args.link_long:
        # Explicitly requested
        enable_linking = True
        logger.info("Longitudinal linking enabled (--link-long specified)")
    elif args.no_link_long:
        # Explicitly disabled
        enable_linking = False
        if multiple_sessions:
            logger.warning(
                "Multiple sessions detected but linking disabled (--no-link-long specified)"
            )
    elif multiple_sessions:
        # Auto-enable when multiple sessions found
        enable_linking = True
        logger.info(
            f"Multiple sessions detected ({multiple_sessions}), automatically enabling longitudinal linking"
        )
    else:
        enable_linking = False

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
        logger.info(f"Skipped subjects (fsid-base) provided: {len(skip_set)}")

    # Set list limit globally for summary printing
    setattr(sys.modules[__name__], "_LIST_LIMIT", max(0, int(args.list_limit)))

    out_path = resolve_output_path(args)
    write_qdec(out_path, header, rows)
    logger.info(f"Wrote Qdec file: {out_path}")

    # Optional consistency summary
    summarize_consistency(
        args.bids, args.subjects_dir, participants_rows, participant_col, session_col, timepoints
    )

    # Optional FastSurfer .long symlink verification/creation for FreeSurfer tools compatibility
    if args.verify_long or enable_linking:
        verify_and_link_long(
            args.subjects_dir,
            timepoints,
            link=enable_linking,
            dry_run=args.link_dry_run,
            force=args.link_force,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
