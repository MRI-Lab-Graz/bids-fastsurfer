#!/usr/bin/env python3
"""
Generic FreeSurfer longitudinal extractor: pulls region-level measures from
any of several FreeSurfer output types into ONE unified tidy long-format TSV,
for use with the flexible analysis battery in scripts/flex/R/.

Supported --source values:
  aseg            aseg.stats (subcortical volumes) -- single file per timepoint
  aparc           lh/rh.aparc.stats (Desikan-Killiany cortical parcellation)
  aparc.a2009s    lh/rh.aparc.a2009s.stats (Destrieux cortical parcellation)
  hippo-amygdala  lh/rh.hippoSfVolumes.long.txt + lh/rh.amygNucVolumes.long.txt
                  (segment_subregions hippo-amygdala; both come out together)
  thalamus        ThalamicNuclei.long.<base>.T1.v13.txt (segment_subregions thalamus)
  brainstem       brainstemSsLabels.long.<base>.v13.txt (segment_subregions brainstem)

For aseg/aparc/aparc.a2009s, --metric selects which column to extract
(default: Volume_mm3 for aseg, ThickAvg for aparc/aparc.a2009s -- also
available: SurfArea, GrayVol, NumVert, etc, exactly as named in the stats
file's ColHeaders). Cortical thickness does NOT need an eTIV correction the
way volume/area measures do (thickness doesn't scale with head size) --
eTIV is still extracted and included for reference, but --etiv-covariate is
left as a downstream analysis choice, not decided here.

Output tidy schema (superset of this project's earlier per-domain tables):
  subject_id, session, fsid, hemisphere, source, region_raw, region,
  is_composite, value, etiv

Input layout expected: <SUBJECTS_DIR>/sub-<ID>_ses-<N>.long.sub-<ID>/{mri,stats}/...
(longitudinal FreeSurfer output; override the subject/session directory
pattern with --dir-pattern for other naming conventions).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DEFAULT_DIR_PATTERN = r"^(?P<base>sub-[^/_]+)_(?P<ses>ses-[^/]+)\.long\.(?P=base)$"

DEFAULT_METRIC = {
    "aseg": "Volume_mm3",
    "aparc": "ThickAvg",
    "aparc.a2009s": "ThickAvg",
}

# Anatomical grouping for the FreeSurfer 8.x hippo-amygdala subfield atlas
# (same convention as the earlier per-study extractors).
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


def sanitize_colname(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z_]", "_", name)


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


def find_timepoints(subjects_dir: Path, dir_pattern: str) -> List[Tuple[str, str, Path]]:
    """Return (subject_base, session, dir) for every directory matching dir_pattern."""
    rx = re.compile(dir_pattern)
    out = []
    for d in sorted(subjects_dir.iterdir()):
        if not d.is_dir():
            continue
        m = rx.match(d.name)
        if not m:
            continue
        out.append((m.group("base"), m.group("ses"), d))
    return out


def parse_etiv(stats_file: Path) -> Optional[float]:
    if not stats_file.is_file():
        return None
    with open(stats_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("# Measure EstimatedTotalIntraCranialVol") or line.startswith("# Measure eTIV"):
                parts = line.split(",")
                if len(parts) >= 4:
                    try:
                        return float(parts[3].strip())
                    except ValueError:
                        return None
    return None


# ---------------------------------------------------------------------------
# aseg.stats: single file, `Index SegId NVoxels Volume_mm3 StructName ...`
# ---------------------------------------------------------------------------
def parse_aseg_stats(path: Path, metric: str) -> Dict[str, float]:
    values: Dict[str, float] = {}
    if not path.is_file():
        return values
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            struct_name = parts[4]
            try:
                vol = float(parts[3])
            except ValueError:
                continue
            if metric == "Volume_mm3":
                values[struct_name] = vol
    return values


# ---------------------------------------------------------------------------
# aparc/aparc.a2009s stats: per-hemisphere, `# TableCol N ColHeader X` header
# block followed by a StructName + N numeric-column table.
# ---------------------------------------------------------------------------
def parse_aparc_stats(path: Path, metric: str) -> Dict[str, float]:
    values: Dict[str, float] = {}
    if not path.is_file():
        return values
    col_headers: Dict[int, str] = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        m = re.match(r"#\s*TableCol\s+(\d+)\s+ColHeader\s+(\S+)", line)
        if m:
            col_headers[int(m.group(1))] = m.group(2)
    if not col_headers:
        raise ValueError(f"No '# TableCol N ColHeader X' lines found in {path} -- unexpected stats format")
    metric_col = None
    for idx, name in col_headers.items():
        if name == metric:
            metric_col = idx
            break
    if metric_col is None:
        available = ", ".join(sorted(col_headers.values()))
        raise ValueError(f"--metric '{metric}' not found in {path}. Available columns: {available}")
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        # column 1 is StructName (text); TableCol indices are 1-based and
        # include StructName as column 1, so the Nth numeric field is at
        # parts[metric_col - 1]
        if len(parts) < metric_col:
            continue
        struct_name = parts[0]
        try:
            values[struct_name] = float(parts[metric_col - 1])
        except ValueError:
            continue
    return values


# ---------------------------------------------------------------------------
# segment_subregions single-column-per-line outputs (hippo-amygdala, and by
# the same documented format, thalamus/brainstem)
# ---------------------------------------------------------------------------
def parse_subregion_txt(path: Path) -> Dict[str, float]:
    values: Dict[str, float] = {}
    if not path.is_file():
        return values
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            name, val = parts
            try:
                values[name] = float(val)
            except ValueError:
                continue
    return values


def extract_source(source: str, d: Path, metric: str) -> List[Tuple[str, str, float, bool]]:
    """Return list of (hemisphere, region_raw, value, is_composite) for one subject/session dir."""
    rows: List[Tuple[str, str, float, bool]] = []

    if source == "aseg":
        vals = parse_aseg_stats(d / "stats" / "aseg.stats", metric)
        for name, v in vals.items():
            hemi = "lh" if name.startswith("Left-") else ("rh" if name.startswith("Right-") else "midline")
            rows.append((hemi, name, v, False))

    elif source in ("aparc", "aparc.a2009s"):
        for hemi in ("lh", "rh"):
            stats_path = d / "stats" / f"{hemi}.{source}.stats"
            vals = parse_aparc_stats(stats_path, metric)
            for name, v in vals.items():
                rows.append((hemi, name, v, False))

    elif source == "hippo-amygdala":
        for hemi in ("lh", "rh"):
            hippo_vals = parse_subregion_txt(d / "mri" / f"{hemi}.hippoSfVolumes.long.txt")
            for name, v in hippo_vals.items():
                rows.append((hemi, f"hippo:{name}", v, hippo_region_for(name)[1]))
            amyg_vals = parse_subregion_txt(d / "mri" / f"{hemi}.amygNucVolumes.long.txt")
            for name, v in amyg_vals.items():
                rows.append((hemi, f"amygdala:{name}", v, name == "Whole_amygdala"))

    elif source == "thalamus":
        # documented FreeSurfer naming; not locally verified against this
        # project's data (never run locally -- only available via the
        # remote-processed concat tables used elsewhere in this project)
        for cand in d.glob("mri/ThalamicNuclei*.txt"):
            vals = parse_subregion_txt(cand)
            for name, v in vals.items():
                hemi = "lh" if name.startswith("Left-") else ("rh" if name.startswith("Right-") else "midline")
                rows.append((hemi, name, v, name in ("Left-Whole_thalamus", "Right-Whole_thalamus")))

    elif source == "brainstem":
        # same caveat as thalamus above
        for cand in d.glob("mri/brainstemSsLabels*.txt"):
            vals = parse_subregion_txt(cand)
            for name, v in vals.items():
                rows.append(("midline", name, v, name == "Whole_brainstem"))

    else:
        raise ValueError(f"Unknown --source '{source}'")

    return rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("subjects_dir", help="Path to the FreeSurfer SUBJECTS_DIR")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--source", required=True, action="append", choices=sorted(DEFAULT_METRIC.keys() | {"hippo-amygdala", "thalamus", "brainstem"}),
                   help="FreeSurfer output type to extract (repeatable, e.g. --source aparc.a2009s --source aseg)")
    p.add_argument("--metric", default=None,
                   help="Column to extract for aseg/aparc/aparc.a2009s sources [default: Volume_mm3 for aseg, ThickAvg for aparc*]")
    p.add_argument("--dir-pattern", default=DEFAULT_DIR_PATTERN,
                   help=f"Regex (named groups 'base','ses') matching per-timepoint directories [default: {DEFAULT_DIR_PATTERN}]")
    p.add_argument("--tidy-name", default="tidy.tsv", help="Output tidy TSV filename [default: %(default)s]")
    p.add_argument("--missing-log", default="missing.tsv", help="Filename (in --outdir) logging subjects/sessions with no data for a source")
    args = p.parse_args()

    subjects_dir = Path(args.subjects_dir)
    if not subjects_dir.is_dir():
        print(f"Error: {subjects_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    timepoints = find_timepoints(subjects_dir, args.dir_pattern)
    if not timepoints:
        print(f"No directories matching --dir-pattern found under {subjects_dir}", file=sys.stderr)
        sys.exit(1)

    tidy_rows: List[Dict[str, object]] = []
    missing: List[Dict[str, str]] = []

    for source in args.source:
        metric = args.metric or DEFAULT_METRIC.get(source, "")
        for subject_base, session, d in timepoints:
            fsid = f"{subject_base}_{session}"
            etiv = parse_etiv(d / "stats" / "aseg.stats")
            rows = extract_source(source, d, metric)
            if not rows:
                missing.append({"subject_id": subject_base, "session": session, "source": source,
                                 "reason": "no data found for this source"})
                continue
            for hemi, region_raw, value, is_composite in rows:
                region = hippo_base_name(region_raw.split(":", 1)[-1]) if source == "hippo-amygdala" else region_raw
                tidy_rows.append({
                    "subject_id": subject_base, "session": session, "fsid": fsid, "hemisphere": hemi,
                    "source": source, "region_raw": region_raw, "region": region,
                    "is_composite": is_composite, "value": value, "etiv": etiv,
                })

    if not tidy_rows:
        print("No data extracted for any requested source.", file=sys.stderr)
        sys.exit(1)

    tidy_path = outdir / args.tidy_name
    tidy_cols = ["subject_id", "session", "fsid", "hemisphere", "source", "region_raw", "region", "is_composite", "value", "etiv"]
    with open(tidy_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tidy_cols, delimiter="\t")
        w.writeheader()
        w.writerows(tidy_rows)
    print(f"Wrote {tidy_path} ({len(tidy_rows)} rows, sources={args.source})")

    if missing:
        missing_path = outdir / args.missing_log
        with open(missing_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["subject_id", "session", "source", "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(missing)
        print(f"Wrote {missing_path} ({len(missing)} missing subject/session/source combinations)", file=sys.stderr)


if __name__ == "__main__":
    main()
