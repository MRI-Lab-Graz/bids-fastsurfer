#!/usr/bin/env python3
"""
DataLad ingest stage for the flex pipeline: pulls only the small per-subject
text tables that scripts/flex/extract_freesurfer.py needs (stats files and
subregion volumes.txt tables) -- never the .mgz volumes -- for every
datalad-backed dataset referenced by a study config's --config measures.

Reuses the directory discovery and dir_pattern regex from
extract_freesurfer.py (same directory layout, same override mechanism) so
the two stages can never disagree about which subject/session dirs exist.

Mirrors the warn-don't-fail, single-batched-get pattern already proven in
scripts/run_hipsta_segmentation.sh (lines ~116-134): compute the exact file
list up front, skip what's already present, one `datalad get -J <jobs>`
call, and log anything that still can't be fetched rather than aborting the
whole run over one missing subject.

Usage:
  python scripts/flex/ingest.py --config configs/flex/study.pk01.json
  python scripts/flex/ingest.py --config ... --dry-run
  python scripts/flex/ingest.py --config ... --jobs 8
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ConfigError, dataset_path, load_config  # noqa: E402
from extract_freesurfer import DEFAULT_DIR_PATTERN, find_timepoints  # noqa: E402

# Per-source, per-timepoint relative file paths needed by extract_freesurfer.py.
# stats/aseg.stats is always included separately (every source's eTIV lookup
# reads it, regardless of which source is being extracted -- see
# extract_freesurfer.py's main loop).
SOURCE_FILES = {
    "aseg": ["stats/aseg.stats"],
    "aparc": ["stats/lh.aparc.stats", "stats/rh.aparc.stats"],
    "aparc.a2009s": ["stats/lh.aparc.a2009s.stats", "stats/rh.aparc.a2009s.stats"],
    "hippo-amygdala": [
        "mri/lh.hippoSfVolumes.long.txt", "mri/rh.hippoSfVolumes.long.txt",
        "mri/lh.amygNucVolumes.long.txt", "mri/rh.amygNucVolumes.long.txt",
    ],
}
# thalamus/brainstem use a glob in extract_freesurfer.py (filename varies by
# FreeSurfer version); resolved against the actual directory listing below
# rather than hardcoded here.
SOURCE_GLOBS = {
    "thalamus": ["mri/ThalamicNuclei*.txt"],
    "brainstem": ["mri/brainstemSsLabels*.txt"],
}
HIPSTA_GLOBS = [
    "thickness/lh.grid-segments-z.csv", "thickness/rh.grid-segments-z.csv",
    "thickness/lh.mid-surface.hsf.csv", "thickness/rh.mid-surface.hsf.csv",
    "hipsta-status_lh.txt", "hipsta-status_rh.txt",
]


def files_needed_for(d: Path, source: str) -> List[str]:
    if source == "hipsta":
        return HIPSTA_GLOBS
    rels = list(SOURCE_FILES.get(source, []))
    for pattern in SOURCE_GLOBS.get(source, []):
        rels.extend(str(p.relative_to(d)) for p in d.glob(pattern))
    return rels


def is_git_annex_repo(d: Path) -> bool:
    r = subprocess.run(["git", "rev-parse", "--git-dir"], cwd=d, capture_output=True, text=True)
    return r.returncode == 0


def annex_missing_content(ds_path: Path, rel_paths: List[str]) -> Set[str]:
    """Ask git-annex which of these tracked paths lack local content.

    Deliberately NOT a symlink/size heuristic: this repo's annex is in v7
    "unlocked"/adjusted-branch mode, where an un-fetched annexed file is a
    small ordinary-looking text file (a pointer to /annex/objects/<key>),
    not a broken symlink -- a symlink-based check silently finds nothing to
    fetch. `git annex find --batch --not --in=here` is mode-agnostic and
    authoritative regardless of how the annex happens to be checked out.
    """
    if not rel_paths:
        return set()
    proc = subprocess.run(
        ["git", "annex", "find", "--not", "--in=here", "--batch"],
        cwd=ds_path, input="\n".join(rel_paths) + "\n",
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        print(f"WARNING: git annex find failed in {ds_path}: {proc.stderr.strip()}", file=sys.stderr)
    return {line for line in proc.stdout.splitlines() if line}


def plan_dataset(ds_path: Path, dir_pattern: str, sources: Set[str], jobs: int) -> Tuple[List[Path], List[str]]:
    """Return (files_to_fetch, warnings) for one dataset directory."""
    timepoints = find_timepoints(ds_path, dir_pattern)
    if not timepoints:
        return [], [f"no directories matching dir_pattern under {ds_path}"]

    candidates: List[Path] = []
    seen: Set[Path] = set()
    for _base, _ses, d in timepoints:
        rels = set(SOURCE_FILES["aseg"])  # eTIV is always read
        for source in sources:
            rels.update(files_needed_for(d, source))
        for rel in rels:
            fpath = d / rel
            if fpath in seen:
                continue
            seen.add(fpath)
            # Present on disk at all (pointer file, symlink, or real
            # content) -- git-annex is the authority on whether the
            # *content* is here; a path missing entirely is a genuine
            # missing input, not something to fetch.
            if fpath.is_file() or fpath.is_symlink():
                candidates.append(fpath)

    rel_candidates = [str(f.relative_to(ds_path)) for f in candidates]
    missing_rel = annex_missing_content(ds_path, rel_candidates)
    to_fetch = [f for f, rel in zip(candidates, rel_candidates) if rel in missing_rel]
    return to_fetch, []


def datalad_get(ds_path: Path, files: List[Path], jobs: int) -> Tuple[int, int]:
    rels = [str(f.relative_to(ds_path)) for f in files]
    cmd = ["datalad", "get", "-J", str(jobs), *rels]
    proc = subprocess.run(cmd, cwd=ds_path)
    if proc.returncode != 0:
        return len(files), 0
    return 0, len(files)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="Path to a flex study.json")
    p.add_argument("--jobs", type=int, default=4, help="Parallel datalad get jobs [default: %(default)s]")
    p.add_argument("--dry-run", action="store_true", help="Print the file list/estimate without fetching")
    p.add_argument("--missing-log", default=None,
                    help="Where to log unfetchable files [default: <output.root>/logs/ingest_missing.tsv]")
    args = p.parse_args()

    try:
        cfg = load_config(args.config)
    except ConfigError as e:
        print(f"INVALID CONFIG: {e}", file=sys.stderr)
        raise SystemExit(1)

    # Group requested sources by dataset, restricted to datalad-backed ones.
    by_dataset: Dict[str, Set[str]] = {}
    for m in cfg["measures"]:
        ds_name = m["dataset"]
        ds_cfg = cfg["datasets"][ds_name]
        if not ds_cfg.get("datalad", False):
            continue
        by_dataset.setdefault(ds_name, set()).add(m["source"])

    if not by_dataset:
        print("No datalad-backed datasets referenced by this config's measures -- nothing to ingest.")
        return

    missing_log_path = Path(args.missing_log) if args.missing_log else \
        Path(cfg["output"]["root"]) / "logs" / "ingest_missing.tsv"

    all_unfetchable: List[Dict[str, str]] = []
    total_planned = 0
    total_fetched = 0

    for ds_name, sources in by_dataset.items():
        ds_path = dataset_path(cfg, ds_name)
        ds_cfg = cfg["datasets"][ds_name]
        dir_pattern = ds_cfg.get("dir_pattern", DEFAULT_DIR_PATTERN)

        if not ds_path.is_dir():
            print(f"WARNING: dataset '{ds_name}' path does not exist: {ds_path}", file=sys.stderr)
            continue
        if not is_git_annex_repo(ds_path):
            print(f"NOTE: '{ds_name}' is marked datalad=true but {ds_path} is not a git repo -- skipping ingest for it")
            continue

        to_fetch, warnings = plan_dataset(ds_path, dir_pattern, sources, args.jobs)
        for w in warnings:
            print(f"WARNING [{ds_name}]: {w}", file=sys.stderr)

        print(f"[{ds_name}] sources={sorted(sources)}: {len(to_fetch)} file(s) need fetching")
        total_planned += len(to_fetch)

        if args.dry_run:
            for f in to_fetch[:20]:
                print(f"  would fetch: {f.relative_to(ds_path)}")
            if len(to_fetch) > 20:
                print(f"  ... and {len(to_fetch) - 20} more")
            continue

        if not to_fetch:
            continue

        n_failed, n_ok = datalad_get(ds_path, to_fetch, args.jobs)
        if n_failed:
            print(f"WARNING [{ds_name}]: datalad get reported a failure; re-checking with git-annex "
                  "to isolate which inputs are actually still missing", file=sys.stderr)
            rel_to_fetch = [str(f.relative_to(ds_path)) for f in to_fetch]
            still_missing = annex_missing_content(ds_path, rel_to_fetch)
            for f, rel in zip(to_fetch, rel_to_fetch):
                if rel in still_missing:
                    all_unfetchable.append({"dataset": ds_name, "path": rel,
                                             "reason": "datalad get failed / content unavailable"})
            total_fetched += len(to_fetch) - len(still_missing)
        else:
            total_fetched += n_ok

    if args.dry_run:
        print(f"\nDRY RUN: {total_planned} file(s) would be fetched across {len(by_dataset)} dataset(s).")
        return

    if all_unfetchable:
        missing_log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(missing_log_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["dataset", "path", "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(all_unfetchable)
        print(f"WARNING: {len(all_unfetchable)} file(s) could not be fetched; logged to {missing_log_path}",
              file=sys.stderr)
        print("Affected subjects/sessions will be logged as missing by the extract stage rather than failing the run.")

    print(f"\nIngest done: {total_fetched}/{total_planned} file(s) fetched.")


if __name__ == "__main__":
    main()
