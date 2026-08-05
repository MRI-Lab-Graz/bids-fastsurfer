#!/usr/bin/env python3
"""
Top-level runner for the flex pipeline: validates a study.json, then runs
ingest (datalad get) -> extract (tidy tables) -> analyse (R test battery)
in order, or a single stage on its own. This is what a new study actually
invokes -- see docs/FLEX_PIPELINE.md for the "porting a new study" walkthrough.

Usage:
  python scripts/flex/run_study.py --config configs/flex/study.pk01.json
  python scripts/flex/run_study.py --config ... --validate-only
  python scripts/flex/run_study.py --config ... --stage ingest --dry-run
  python scripts/flex/run_study.py --config ... --stage extract
  python scripts/flex/run_study.py --config ... --stage analyse --only lmm,moderator
  python scripts/flex/run_study.py --config ... --stage analyse --only lmm --measures hippo,cortex

Each stage is independently re-runnable: ingest/extract always re-check
their own inputs (ingest via git-annex, extract by re-parsing the source
files) so a repeat run is safe and just redoes cheap work; analyse always
overwrites its own output directory. There is no cross-stage skip-if-done
tracking -- re-run only the stage you need via --stage.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ConfigError, load_config, resolve_path  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
FLEX_DIR = Path(__file__).resolve().parent


def run_ingest(config_path: str, jobs: int, dry_run: bool) -> None:
    cmd = [sys.executable, str(FLEX_DIR / "ingest.py"), "--config", config_path, "--jobs", str(jobs)]
    if dry_run:
        cmd.append("--dry-run")
    print(f"\n=== ingest ===\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def run_extract(config_path: str) -> None:
    cmd = [sys.executable, str(FLEX_DIR / "extract_freesurfer.py"), "--config", config_path]
    print(f"\n=== extract ===\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def run_analyse(config_path: str, only: str | None, measures: str | None) -> None:
    cmd = ["Rscript", str(FLEX_DIR / "R" / "run_analysis.R"), "--config", config_path]
    if only:
        cmd += ["--only", only]
    if measures:
        cmd += ["--measures", measures]
    print(f"\n=== analyse ===\n$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def write_manifest(cfg: dict, stages_run: list[str]) -> None:
    out_root = resolve_path(cfg, cfg["output"]["root"])
    out_root.mkdir(parents=True, exist_ok=True)
    manifest = {
        "study_id": cfg["study"]["id"],
        "config_path": cfg["_config_path"],
        "stages_run": stages_run,
        "run_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    with open(out_root / "run_manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nWrote {out_root / 'run_manifest.json'}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, help="Path to a flex study.json")
    p.add_argument("--stage", choices=["ingest", "extract", "analyse"], default=None,
                    help="Run only this stage [default: ingest, extract, analyse in order]")
    p.add_argument("--validate-only", action="store_true", help="Validate the config and exit")
    p.add_argument("--dry-run", action="store_true", help="Ingest stage only: plan without fetching")
    p.add_argument("--jobs", type=int, default=4, help="Ingest stage: parallel datalad get jobs [default: %(default)s]")
    p.add_argument("--only", default=None, help="Analyse stage: comma-separated test names to run")
    p.add_argument("--measures", default=None, help="Analyse stage: comma-separated measure names to restrict to")
    args = p.parse_args()

    try:
        cfg = load_config(args.config)
    except ConfigError as e:
        print(f"INVALID CONFIG: {e}", file=sys.stderr)
        raise SystemExit(1)
    print(f"Config OK: {args.config} ({cfg['study']['id']}) -- "
          f"{len(cfg['measures'])} measures, {len(cfg.get('analyses', []))} analyses")

    if args.validate_only:
        return

    stages = [args.stage] if args.stage else ["ingest", "extract", "analyse"]
    stages_run = []

    try:
        if "ingest" in stages:
            run_ingest(args.config, args.jobs, args.dry_run)
            stages_run.append("ingest")
        if "extract" in stages and not args.dry_run:
            run_extract(args.config)
            stages_run.append("extract")
        if "analyse" in stages and not args.dry_run:
            run_analyse(args.config, args.only, args.measures)
            stages_run.append("analyse")
    except subprocess.CalledProcessError as e:
        print(f"\nSTAGE FAILED: {e}", file=sys.stderr)
        raise SystemExit(1)

    if not args.dry_run and stages_run:
        write_manifest(cfg, stages_run)
    print(f"\nDone: stage(s) run = {stages_run or '(none -- dry-run)'}")


if __name__ == "__main__":
    main()
