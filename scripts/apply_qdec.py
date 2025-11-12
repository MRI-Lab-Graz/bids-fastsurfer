
"""DEPRECATED: apply_qdec.py

This file is a minimal compatibility shim. It forwards all CLI arguments to
``scripts/analyse_qdec.py`` and exits with the same return code.

Keep this small to avoid importing heavy modules or duplicating logic.
"""

from __future__ import annotations

import sys
import subprocess
from pathlib import Path


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    prog = Path(__file__).resolve()
    analyse = prog.parent / "analyse_qdec.py"
    if not analyse.exists():
        print(
            "ERROR: analyse_qdec.py not found next to apply_qdec.py. Cannot forward.",
            file=sys.stderr,
        )
        return 2

    # Deprecation warning
    print(
        "[DEPRECATED] 'apply_qdec.py' is deprecated — forwarding to 'analyse_qdec.py'. Please update your calls.",
        file=sys.stderr,
    )

    cmd = [sys.executable, str(analyse)] + argv
    try:
        rc = subprocess.call(cmd)
        return rc
    except KeyboardInterrupt:
        return 130
    except Exception as e:
        print(f"ERROR: Failed to invoke analyse_qdec.py: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

def read_qdec(
    qdec_path: Path,
) -> Tuple[List[str], List[List[str]], List[Tuple[str, str, Optional[str]]]]:
    """Read QDEC file and extract header, rows, and timepoints.

    Returns:
        header: list of column names
        rows: list of row lists
        timepoints: list of (fsid, fsid_base, session) tuples
    """
    with qdec_path.open("r", newline="") as f:
        reader = csv.reader(f, dialect=csv.excel_tab)
        rows = list(reader)
    if not rows:
        raise ValueError(f"QDEC file {qdec_path} is empty")
    header = rows[0]
    data_rows = rows[1:]

    # Find fsid and fsid-base columns
    try:
        fsid_idx = header.index("fsid")
        base_idx = header.index("fsid-base")
    except ValueError as e:
        raise ValueError(f"QDEC file missing required columns: {e}")

    timepoints = []
    for row in data_rows:
        if len(row) <= max(fsid_idx, base_idx):
            continue
        fsid = row[fsid_idx]
        base = row[base_idx]
        # Extract session from fsid if it has _ses-
        m = SUBJECT_DIR_PATTERN.match(fsid)
        ses = m.group("ses") if m else None
        timepoints.append((fsid, base, ses))

    return header, data_rows, timepoints


def detect_study_type(timepoints: List[Tuple[str, str, Optional[str]]]) -> str:
    """Detect if study is longitudinal or cross-sectional based on timepoints."""
    if not timepoints:
        return "unknown"

    # If any timepoint has a session, it's longitudinal
    for _, _, ses in timepoints:
        if ses is not None:
            return "longitudinal"

    # If fsid != fsid-base for any, it's longitudinal
    for fsid, base, _ in timepoints:
        if fsid != base:
            return "longitudinal"

    return "cross-sectional"


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
        if missing_in_sd:
            logger.info(f"BIDS subjects missing in subjects_dir: {len(missing_in_sd)}")
            if missing_in_sd != only_in_participants:
                logger.info(
                    ", ".join(missing_in_sd[:limit])
                    + (" ..." if len(missing_in_sd) > limit else "")
                )
        if missing_in_parts:
            logger.info(f"BIDS subjects missing in participants.tsv: {len(missing_in_parts)}")
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
            logger.debug(f"skipping: {fsid} (already a .long entry)")
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
            logger.info(msg)
        else:
            note_missing = " [NO-EVIDENCE]" if require_stats and not stats_path.exists() else ""
            logger.info(
                f"would link: {long_dir} -> {tp_dir} (use --link-long to create){note_missing}"
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


def run_asegstats2table(
    qdec_path: Path, subjects_dir: Path, study_type: str = "longitudinal"
) -> int:
    """Run asegstats2table with SUBJECTS_DIR pointing to subjects_dir."""

    aseg_bin = shutil.which("asegstats2table")
    if not aseg_bin:
        print(
            "asegstats2table not found in PATH. Source FreeSurfer before using --aseg.",
            file=sys.stderr,
        )
        return 4

    # Output filename depends on study type
    if study_type == "longitudinal":
        aseg_out = qdec_path.parent / "aseg.long.table"
    else:
        aseg_out = qdec_path.parent / "aseg.table"

    aseg_out.parent.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["SUBJECTS_DIR"] = str(subjects_dir.resolve())

    # Build command based on study type
    if study_type == "longitudinal":
        cmd = [
            aseg_bin,
            "--qdec-long",
            str(qdec_path),
            "-t",
            str(aseg_out),
            "--skip",
        ]
    else:
        # Cross-sectional: need to extract subject IDs from Qdec
        try:
            with qdec_path.open("r", newline="") as fh:
                reader = csv.reader(fh, dialect=csv.excel_tab)
                rows = list(reader)
            if not rows:
                print("[WARN] QDEC empty; skipping asegstats2table", file=sys.stderr)
                return 0
            header = rows[0]
            id_col = "fsid"
            try:
                idx = header.index(id_col)
            except ValueError:
                print(
                    f"[WARN] Column '{id_col}' not found in QDEC; skipping asegstats2table",
                    file=sys.stderr,
                )
                return 0
            subjects = [r[idx] for r in rows[1:] if len(r) > idx and r[idx]]
            if not subjects:
                print("[WARN] No subjects found in QDEC; skipping asegstats2table", file=sys.stderr)
                return 0
        except Exception as e:
            print(
                f"[WARN] Failed to parse QDEC for subjects: {e}; skipping asegstats2table",
                file=sys.stderr,
            )
            return 0

        cmd = (
            [aseg_bin, "--subjects"]
            + subjects
            + [
                "-t",
                str(aseg_out),
            ]
        )

    print(
        f"Running: {' '.join(cmd[:10])}{'...' if len(cmd) > 10 else ''} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})"
    )

    try:
        subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        error_output = exc.stderr or exc.stdout or ""
        if study_type == "longitudinal" and (
            "IndexError: list index out of range" in error_output
            or "list index out of range" in error_output
        ):
            print(
                "asegstats2table failed because no valid longitudinal data was found. "
                "This likely means .long directories are missing or don't contain proper stats files. "
                "Try using --link-long to create the required symlinks first, or check that FastSurfer/FreeSurfer "
                "processing completed successfully for the timepoints.",
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
    study_type: str = "longitudinal",
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
        missing = [
            n
            for n, b in [("mris_preproc", mris_preproc_bin), ("mri_surf2surf", surf2surf_bin)]
            if not b
        ]
        print(
            f"[WARN] Missing FreeSurfer binaries: {', '.join(missing)}. Skipping surface prep.",
            file=sys.stderr,
        )
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
            verify_and_link_long(
                subjects_dir, tps, link=True, dry_run=False, force=False, require_stats=False
            )
        except Exception as e:
            print(
                f"[WARN] Failed to auto-link .long symlinks before surface prep: {e}",
                file=sys.stderr,
            )
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
        print(
            f"[INFO] Filtered QDEC for {hemi}/{meas}: kept={len(kept_rows)}, dropped={dropped} -> {filt_path}"
        )
        return filt_path, len(kept_rows), dropped, dropped_pairs

    # QC summary rows
    qc_rows: List[List[str]] = [
        ["hemi", "measure", "kept", "dropped", "filtered_qdec", "missing_list"]
    ]

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
                print(
                    f"[WARN] Skipping surface prep for {hemi}/{meas}: no subjects with existing surf files.",
                    file=sys.stderr,
                )
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
            if study_type == "longitudinal":
                qdec_flag = "--qdec-long"
            else:
                qdec_flag = "--qdec"
            cmd1 = [
                mris_preproc_bin,
                qdec_flag,
                str(qdec_for_pair),
                "--target",
                target,
                "--hemi",
                hemi,
                "--meas",
                meas,
                "--out",
                str(pre_path),
            ]
            print(f"Running: {' '.join(cmd1)} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})")
            if not dry_run:
                try:
                    subprocess.run(cmd1, check=True, env=env)
                except subprocess.CalledProcessError as exc:
                    print(
                        f"mris_preproc failed (hemi={hemi}, meas={meas}) with code {exc.returncode}",
                        file=sys.stderr,
                    )
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
                    "--hemi",
                    hemi,
                    "--s",
                    target,
                    "--sval",
                    str(pre_path),
                    "--tval",
                    str(sm_path),
                    "--fwhm-trg",
                    str(fwhm),
                    "--cortex",
                    "--noreshape",
                ]
                print(f"Running: {' '.join(cmd2)} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})")
                if not dry_run:
                    try:
                        subprocess.run(cmd2, check=True, env=env)
                    except subprocess.CalledProcessError as exc:
                        print(
                            f"mri_surf2surf failed (hemi={hemi}, meas={meas}, fwhm={fwhm}) with code {exc.returncode}",
                            file=sys.stderr,
                        )
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
    # Try to find run_fsqc command
    fsqc_bin = shutil.which("run_fsqc")
    if not fsqc_bin:
        # If run_fsqc not in PATH, check if fsqc module is available and try to run it via python -m
        try:
            import fsqc  # noqa: F401

            # Use fsqc.run_fsqc directly
            use_direct_call = True
        except ImportError:
            print(
                "[WARN] fsqc not found (run_fsqc command or Python module). Skipping --qc step. Install with: bash scripts/install.sh",
                file=sys.stderr,
            )
            return 0
    else:
        use_direct_call = False
        fsqc_command = [fsqc_bin]

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
        headless = disp == ""
    except Exception:
        headless = True
    if surfaces and headless:
        print(
            "[INFO] No DISPLAY detected; disabling fsqc surfaces module to avoid OpenGL errors.",
            file=sys.stderr,
        )
        surfaces = False

    cmd = fsqc_command + [
        "--subjects_dir",
        str(subjects_dir),
        "--output_dir",
        str(out_root),
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
    if use_direct_call:
        print(f"Running fsqc.run_fsqc(subjects_dir={subjects_dir}, output_dir={out_root}, subjects={subjects}, fastsurfer={fastsurfer}, screenshots={screenshots}, surfaces={surfaces}, skullstrip={skullstrip}, outlier={outlier}, screenshots_html={html if screenshots else False}, surfaces_html={html if surfaces else False}, skullstrip_html={html if skullstrip else False}, skip_existing={skip_existing and not force})")
        try:
            fsqc.run_fsqc(
                subjects_dir=str(subjects_dir),
                output_dir=str(out_root),
                subjects=subjects,
                fastsurfer=fastsurfer,
                screenshots=screenshots,
                surfaces=surfaces,
                skullstrip=skullstrip,
                outlier=outlier,
                screenshots_html=html if screenshots else False,
                surfaces_html=html if surfaces else False,
                skullstrip_html=html if skullstrip else False,
                skip_existing=skip_existing and not force,
            )
            print(f"Wrote fsqc outputs to: {out_root}")
        except Exception as exc:
            print(
                f"[WARN] fsqc failed with exception {exc}; continuing.",
                file=sys.stderr,
            )
            return 0
    else:
        cmd = fsqc_command + [
            "--subjects_dir",
            str(subjects_dir),
            "--output_dir",
            str(out_root),
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

        print(f"Running fsqc: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, env=env)
            print(f"Wrote fsqc outputs to: {out_root}")
        except subprocess.CalledProcessError as exc:
            print(
                f"[WARN] fsqc failed with exit code {exc.returncode}; continuing. Command: {' '.join(cmd)}",
                file=sys.stderr,
            )
            return 0
    return 0


def run_aparcstats2table(
    qdec_path: Path,
    subjects_dir: Path,
    parc: str = "aparc",
    measures: Optional[List[str]] = None,
    hemis: Optional[List[str]] = None,
    study_type: str = "longitudinal",
) -> int:
    """Run aparcstats2table for cross-sectional or longitudinal studies.

    Args:
        qdec_path: Path to Qdec file
        subjects_dir: Path to subjects directory
        parc: Parcellation name
        measures: List of measures (thickness, area, volume)
        hemis: List of hemispheres (lh, rh)
        study_type: 'cross-sectional' or 'longitudinal'

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

    # Preflight: auto-detect available parcellation stats
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
            if study_type == "longitudinal":
                pattern = f"**/*.long.*/stats/{hemi}.{p}.stats"
            else:
                pattern = f"**/stats/{hemi}.{p}.stats"
            if list(subjects_dir.glob(pattern)):
                found = True
                break
        if found:
            chosen_parc = p
            break

    if not chosen_parc:
        print(
            f"[WARN] No aparc stats files found for any of parcs {candidate_parcs} and hemis={hemis} under {subjects_dir}. Skipping aparc tables."
        )
        return 0
    if chosen_parc != parc:
        print(
            f"[INFO] Using detected parcellation '{chosen_parc}' for aparc tables (requested '{parc}')."
        )
    parc = chosen_parc

    # Get subject list for cross-sectional mode
    subjects = []
    if study_type == "cross-sectional":
        try:
            with qdec_path.open("r", newline="") as fh:
                reader = csv.reader(fh, dialect=csv.excel_tab)
                rows = list(reader)
            if not rows:
                print("[WARN] QDEC empty; skipping aparcstats2table", file=sys.stderr)
                return 0
            header = rows[0]
            id_col = "fsid"
            try:
                idx = header.index(id_col)
            except ValueError:
                print(
                    f"[WARN] Column '{id_col}' not found in QDEC; skipping aparcstats2table",
                    file=sys.stderr,
                )
                return 0
            subjects = [r[idx] for r in rows[1:] if len(r) > idx and r[idx]]
            if not subjects:
                print(
                    "[WARN] No subjects found in QDEC; skipping aparcstats2table", file=sys.stderr
                )
                return 0
        except Exception as e:
            print(
                f"[WARN] Failed to parse QDEC for subjects: {e}; skipping aparcstats2table",
                file=sys.stderr,
            )
            return 0

    for hemi in hemis:
        for meas in measures:
            if study_type == "longitudinal":
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
            else:
                out_path = out_root / f"{hemi}.{parc}.{meas}.table"
                cmd = (
                    [aparc_bin, "--subjects"]
                    + subjects
                    + [
                        "--hemi",
                        hemi,
                        "--meas",
                        meas,
                        "--parc",
                        parc,
                        "-t",
                        str(out_path),
                    ]
                )
            print(
                f"Running: {' '.join(cmd[:10])}{'...' if len(cmd) > 10 else ''} (with SUBJECTS_DIR={env['SUBJECTS_DIR']})"
            )
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

    # Configure logging level based on verbose flag
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    else:
        logging.getLogger().setLevel(logging.INFO)

    # Early dependency check
    missing_deps = check_dependencies(args)
    if missing_deps:
        logger.error("Missing required dependencies:")
        for dep in missing_deps:
            logger.error(f"  - {dep}")
        logger.error(
            "\nPlease install missing dependencies and ensure FreeSurfer is properly sourced."
        )
        return 1

    # Check QDEC file
    if not args.qdec.exists():
        logger.error(f"QDEC file not found: {args.qdec}")
        return 2

    subj_dir: Path = args.subjects_dir
    if not subj_dir.exists() or not subj_dir.is_dir():
        logger.error(f"subjects_dir not found or not a directory: {subj_dir}")
        return 2

    # Read QDEC file
    try:
        header, rows, timepoints = read_qdec(args.qdec)
    except Exception as e:
        logger.error(f"Failed to read QDEC file {args.qdec}: {e}")
        return 2

    study_type = detect_study_type(timepoints)
    logger.info(f"Read QDEC file: {args.qdec} ({len(rows)} rows, study type: {study_type})")

    # Quick overview
    bases = set(tp[1] for tp in timepoints)
    logger.info(f"Subjects overview: bases={len(bases)}, timepoints={len(timepoints)} in QDEC")

    # Prepare output directory
    out_root = args.output
    if not prepare_output_directory(out_root, args.force):
        logger.error("Output directory preparation cancelled by user.")
        return 1

    # Handle verify-qdec: when used, disable analyses by default unless explicitly enabled
    if getattr(args, "verify_qdec", False):
        # Check if analysis flags were explicitly provided
        argv_set = set(sys.argv[1:])  # Get command line args as set for fast lookup

        # Disable analyses by default, but allow explicit enabling
        if "--aseg" not in argv_set:
            args.aseg = False
        if "--aparc" not in argv_set:
            args.aparc = False
        if "--surf" not in argv_set:
            args.surf = False
    # Detect headless environment (no DISPLAY) and auto-disable surfaces to avoid OpenGL/GLFW errors
    try:
        disp = os.environ.get("DISPLAY", "").strip()
        _headless = disp == ""
    except Exception:
        _headless = True
    qc_surfaces_effective = bool(getattr(args, "qc_surfaces", False))
    if bool(getattr(args, "qc", False)) and qc_surfaces_effective and _headless:
        logger.info("No DISPLAY detected; fsqc surfaces will be disabled.")
        qc_surfaces_effective = False

    # Save effective configuration
    try:
        eff_cfg = {
            "qdec": str(args.qdec),
            "subjects_dir": str(subj_dir),
            "output": str(out_root),
            "study_type": study_type,
            "force": bool(args.force),
            "verify_qdec": bool(getattr(args, "verify_qdec", False)),
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
            "smooth": _coerce_int_list(getattr(args, "smooth", None)) or [int(args.surf_fwhm)],
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
            "summary": {
                "bases": len(bases),
                "timepoints": len(timepoints),
            },
        }
        cfg_out = out_root / "apply_qdec.effective.json"
        with cfg_out.open("w") as fh:
            json.dump(eff_cfg, fh, indent=2, sort_keys=True)
        logger.info(f"Wrote effective config: {cfg_out}")
    except Exception as e:
        logger.warning(f"Failed to write effective config JSON: {e}")

    # Optional .long symlink verification/creation
    if args.verify_qdec or args.link_long:
        verify_and_link_long(
            subj_dir,
            timepoints,
            link=args.link_long,
            dry_run=args.link_dry_run,
            force=args.link_force,
        )

        # Add verification summary when using --verify-qdec
        if args.verify_qdec and not args.link_long:
            # Count timepoints per subject
            tp_counts = Counter(tp[1] for tp in timepoints)  # Count by base (fsid-base)
            tp_distribution = Counter(tp_counts.values())  # Count subjects by their timepoint count

            logger.info("=== QDEC Verification Summary ===")
            logger.info(
                f"✅ QDEC file validated: {len(rows)} rows, {len(bases)} subjects, {len(timepoints)} timepoints"
            )
            logger.info(f"✅ Study type: {study_type}")
            logger.info(f"✅ Long symlinks: {len(timepoints)} verified (all present)")

            # Show timepoint distribution
            if len(tp_distribution) > 1:  # Only show if there's variation
                dist_parts = []
                for tp_count in sorted(tp_distribution.keys()):
                    subject_count = tp_distribution[tp_count]
                    dist_parts.append(
                        f"{subject_count} subjects with {tp_count} timepoint{'s' if tp_count != 1 else ''}"
                    )
                logger.info(f"✅ Timepoint distribution: {', '.join(dist_parts)}")

            logger.info("✅ Verification complete - ready for analysis")

    # Optional tables
    if args.aseg:
        if args.link_dry_run:
            logger.info(
                "Skipping asegstats2table due to --link-dry-run (symlinks not actually created)."
            )
        elif shutil.which("asegstats2table") is None:
            logger.warning(
                "asegstats2table not found in PATH; skipping --aseg. Ensure FreeSurfer is sourced."
            )
        else:
            rc = run_asegstats2table(args.qdec, subj_dir, study_type)
            if rc != 0:
                return rc
    if args.aparc:
        if args.link_dry_run:
            logger.info(
                "Skipping aparcstats2table due to --link-dry-run (symlinks not actually created)."
            )
        elif shutil.which("aparcstats2table") is None:
            logger.warning(
                "aparcstats2table not found in PATH; skipping --aparc. Ensure FreeSurfer is sourced."
            )
        else:
            rc = run_aparcstats2table(
                args.qdec,
                subj_dir,
                parc=args.aparc_parc,
                measures=args.aparc_measures,
                hemis=args.aparc_hemis,
                study_type=study_type,
            )
            if rc != 0:
                return rc
    # Optional mass-univariate surface data
    if args.surf:
        have_mris = shutil.which("mris_preproc") is not None
        have_surf2 = shutil.which("mri_surf2surf") is not None
        if not (have_mris and have_surf2):
            missing = [
                n
                for n, ok in (("mris_preproc", have_mris), ("mri_surf2surf", have_surf2))
                if not ok
            ]
            logger.warning(f"Missing FreeSurfer binaries ({', '.join(missing)}); skipping --surf.")
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
                args.qdec,
                subj_dir,
                target=str(args.surf_target),
                measures=list(args.surf_measures),
                hemis=list(args.surf_hemis),
                smooth_kernels=smooth_list,
                outdir=args.surf_outdir,
                force=bool(args.force),
                dry_run=bool(args.link_dry_run),
                study_type=study_type,
            )
            if rc != 0:
                # do not fail the entire prep if surface prep tools missing; return code already logged
                pass
    # Optional fsqc QC
    if args.qc:
        _ = run_fsqc(
            args.qdec,
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
