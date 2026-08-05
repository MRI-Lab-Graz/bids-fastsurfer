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
  thalamus        ThalamicNuclei.long.volumes.txt (segment_subregions thalamus)
  brainstem       brainstemSsLabels.long.volumes.txt (segment_subregions brainstem)
  hipsta          grid-segments-z.csv + mid-surface.hsf.csv (hipsta shape/
                  thickness, https://github.com/Deep-MI/Hipsta) -- ported from
                  extract_hipsta_thickness.py; only reachable via --config mode
                  (its directory layout/status-stamp semantics differ from the
                  other sources, see extract_hipsta_source() below)

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

--config mode: pass --config <study.json> (a flex study definition, see
configs/flex/study.schema.json) instead of the positional/--source CLI to
extract every measure the config declares in one pass, across all of its
datasets, writing <output.root>/tidy/<measure_name>.tsv per measure. This is
the mode scripts/flex/run_study.py uses; the plain CLI above remains for
ad hoc single-source extraction and is what the legacy per-domain
scripts/extract_*.py wrappers were replaced by.
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


def strip_hemi_prefix(name: str) -> str:
    """Strip a leading 'Left-'/'Right-' so the same structure's two
    hemispheres share one ROI name (paired with the separate `hemisphere`
    column) -- e.g. aseg's 'Left-Hippocampus'/'Right-Hippocampus' or
    thalamus's 'Left-Whole_thalamus'/'Right-Whole_thalamus' both become
    'Hippocampus'/'Whole_thalamus'. Without this, a --roi-set/roi_set entry
    written once (as for every other source) would silently match neither
    hemisphere's row."""
    for prefix in ("Left-", "Right-"):
        if name.startswith(prefix):
            return name[len(prefix):]
    return name


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
        for cand in d.glob("mri/ThalamicNuclei*.txt"):
            vals = parse_subregion_txt(cand)
            for name, v in vals.items():
                hemi = "lh" if name.startswith("Left-") else ("rh" if name.startswith("Right-") else "midline")
                rows.append((hemi, name, v, name in ("Left-Whole_thalamus", "Right-Whole_thalamus")))

    elif source == "brainstem":
        for cand in d.glob("mri/brainstemSsLabels*.txt"):
            vals = parse_subregion_txt(cand)
            for name, v in vals.items():
                rows.append(("midline", name, v, name == "Whole_brainstem"))

    else:
        raise ValueError(f"Unknown --source '{source}'")

    return rows


TIDY_COLS = ["subject_id", "session", "fsid", "hemisphere", "source", "region_raw", "region", "is_composite", "value", "etiv"]


def load_subject_id_map(path: Path) -> Dict[str, str]:
    """Load a short_id/long_id/match_type TSV (scripts/build_hipsta_subject_id_map.py).
    Rows with an empty long_id (unresolved/ambiguous, deliberately not guessed
    by that script) are omitted -- their short_id is then treated as unmapped."""
    out: Dict[str, str] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("long_id"):
                out[row["short_id"]] = row["long_id"]
    return out


def remap_subject(subject_base: str, id_map: Optional[Dict[str, str]]) -> Optional[str]:
    """Apply an optional short_id->long_id map; returns None if id_map is
    given but has no entry for this subject (caller should log and skip)."""
    if id_map is None:
        return subject_base
    return id_map.get(subject_base)


def extract_single_source(
    subjects_dir: Path, source: str, metric: str, dir_pattern: str,
    sessions_include: Optional[List[str]] = None,
    id_map: Optional[Dict[str, str]] = None,
) -> Tuple[List[Dict[str, object]], List[Dict[str, str]]]:
    """Extract one source across every timepoint under subjects_dir into the
    unified tidy-row schema. Shared by the legacy CLI (main(), possibly
    combining several sources into one file) and --config mode (one call
    per (dataset, source, metric) combination, cached across measures that
    share it -- see run_config_extraction())."""
    timepoints = find_timepoints(subjects_dir, dir_pattern)
    if sessions_include is not None:
        timepoints = [(b, s, d) for b, s, d in timepoints if s in sessions_include]

    tidy_rows: List[Dict[str, object]] = []
    missing: List[Dict[str, str]] = []
    for raw_subject_base, session, d in timepoints:
        subject_base = remap_subject(raw_subject_base, id_map)
        if subject_base is None:
            missing.append({"subject_id": raw_subject_base, "session": session, "source": source,
                             "reason": "no entry for this short ID in subject_id_map"})
            continue
        fsid = f"{subject_base}_{session}"
        etiv = parse_etiv(d / "stats" / "aseg.stats")
        rows = extract_source(source, d, metric)
        if not rows:
            missing.append({"subject_id": subject_base, "session": session, "source": source,
                             "reason": "no data found for this source"})
            continue
        for hemi, region_raw, value, is_composite in rows:
            if source == "hippo-amygdala":
                region = hippo_base_name(region_raw.split(":", 1)[-1])
            elif source in ("aseg", "thalamus"):
                region = strip_hemi_prefix(region_raw)
            else:
                region = region_raw
            tidy_rows.append({
                "subject_id": subject_base, "session": session, "fsid": fsid, "hemisphere": hemi,
                "source": source, "region_raw": region_raw, "region": region,
                "is_composite": is_composite, "value": value, "etiv": etiv,
            })
    return tidy_rows, missing


def write_tsv(path: Path, cols: List[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# hipsta shape/thickness source (ported from scripts/extract_hipsta_thickness.py
# -- see that module's docstring for the input layout and grid/hsf join this
# mirrors exactly; only reachable from --config mode, since its status-stamp
# gating and grid-point aggregation don't fit the per-timepoint stats/mri
# layout the other sources share).
# ---------------------------------------------------------------------------

# Hipsta's internal, reduced FreeSurfer LUT (hipsta/cfg/atlases.py, lut ==
# "freesurfer"). Only presubiculum/subiculum/CA1/CA2+CA3(merged)/molecular_layer
# actually appear in mid-surface.hsf.csv in practice (HSFLIST in atlases.py),
# but CA4/DG are included here too in case a future hipsta version or a
# different --lut run emits them.
HSF_LABELS = {
    0: "unlabeled",
    234: "presubiculum",
    236: "subiculum",
    238: "CA1",
    240: "CA2_CA3",
    242: "CA4",
    244: "DG",
    246: "molecular_layer",
}


def hsf_label_name(value: int) -> str:
    return HSF_LABELS.get(value, f"unknown_label_{value}")


def parse_grid_csv(path: Path) -> Dict[Tuple[int, int], float]:
    """Parse a <hemi>.grid-segments-*.csv into {(x_idx, y_idx): value}."""
    out: Dict[Tuple[int, int], float] = {}
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        y_idx = [int(c[1:]) for c in header[1:]]
        for row in reader:
            if not row or not row[0]:
                continue
            x = int(row[0][1:])
            for y, val in zip(y_idx, row[1:], strict=True):
                if val == "" or val.lower() == "nan":
                    continue
                out[(x, y)] = float(val)
    return out


def parse_hsf_csv(path: Path) -> Dict[Tuple[int, int], int]:
    """Parse a <hemi>.mid-surface.hsf.csv into {(x_idx, y_idx): hsf_label}.

    Grid points at the medial/lateral edge of the parametrization sometimes
    have no label (hipsta writes an empty field there) -- those points are
    omitted, and end up excluded by the (x,y) intersection with
    grid-segments-z.csv.
    """
    out: Dict[Tuple[int, int], int] = {}
    with open(path, newline="") as f:
        for row in csv.reader(f):
            if not row:
                continue
            x, y, label = row
            if label == "" or label.lower() == "nan":
                continue
            out[(int(round(float(x))), int(round(float(y))))] = int(round(float(label)))
    return out


def hemi_status_ok(session_dir: Path, hemi: str) -> Tuple[bool, str]:
    stamp = session_dir / f"hipsta-status_{hemi}.txt"
    if not stamp.is_file():
        return False, "hipsta-status file not found (segmentation not run or in progress)"
    with open(stamp) as f:
        first_line = f.readline().strip()
    if first_line == "status: OK":
        return True, ""
    return False, f"hipsta-status reports failure ({first_line or 'empty status file'})"


def extract_hipsta(
    hipsta_dir: Path, dir_pattern: str, sessions_include: Optional[List[str]] = None,
    id_map: Optional[Dict[str, str]] = None,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]], List[Dict[str, str]]]:
    """Returns (subfield_rows, grid_rows, missing_rows).

    subfield_rows follow the unified TIDY_COLS schema plus median_value/
    sd_value/n_vertices (grid points aggregated by subfield label per
    subject/session/hemisphere; etiv is left blank -- hipsta thickness
    doesn't use an eTIV covariate, see the "thickness" profile). grid_rows
    are one row per subject/session/hemisphere/grid-vertex, for
    test_pointwise.R.
    """
    timepoints = find_timepoints(hipsta_dir, dir_pattern)
    if sessions_include is not None:
        timepoints = [(b, s, d) for b, s, d in timepoints if s in sessions_include]

    grid_rows: List[Dict[str, object]] = []
    subfield_agg: Dict[Tuple[str, str, str], Dict[str, List[float]]] = {}
    missing: List[Dict[str, str]] = []

    for raw_subject_base, session, session_dir in timepoints:
        subject_base = remap_subject(raw_subject_base, id_map)
        if subject_base is None:
            missing.append({"subject_id": raw_subject_base, "session": session, "hemisphere": "both",
                             "reason": "no entry for this short ID in subject_id_map"})
            continue
        fsid = f"{subject_base}_{session}"
        for hemi in ("lh", "rh"):
            ok, reason = hemi_status_ok(session_dir, hemi)
            if not ok:
                missing.append({"subject_id": subject_base, "session": session, "hemisphere": hemi, "reason": reason})
                continue

            grid_file = session_dir / "thickness" / f"{hemi}.grid-segments-z.csv"
            hsf_file = session_dir / "thickness" / f"{hemi}.mid-surface.hsf.csv"
            if not grid_file.is_file() or not hsf_file.is_file():
                missing.append({"subject_id": subject_base, "session": session, "hemisphere": hemi,
                                 "reason": "status OK but thickness/grid-segments-z.csv or mid-surface.hsf.csv missing"})
                continue

            try:
                thickness = parse_grid_csv(grid_file)
                hsf = parse_hsf_csv(hsf_file)
            except (OSError, ValueError) as e:
                missing.append({"subject_id": subject_base, "session": session, "hemisphere": hemi,
                                 "reason": f"failed to parse thickness/hsf files: {e}"})
                continue

            common_idx = sorted(set(thickness) & set(hsf))
            if not common_idx:
                missing.append({"subject_id": subject_base, "session": session, "hemisphere": hemi,
                                 "reason": "grid-segments-z.csv and mid-surface.hsf.csv share no (x,y) index"})
                continue

            for x, y in common_idx:
                subfield = hsf_label_name(hsf[(x, y)])
                val = thickness[(x, y)]
                grid_rows.append({
                    "subject_id": subject_base, "session": session, "fsid": fsid, "hemisphere": hemi,
                    "x": x, "y": y, "thickness": val, "hsf_label": hsf[(x, y)], "subfield": subfield,
                })
                key = (subject_base, session, hemi)
                subfield_agg.setdefault(key, {}).setdefault(subfield, []).append(val)

    def mean(xs: List[float]) -> float:
        return sum(xs) / len(xs)

    def median(xs: List[float]) -> float:
        s = sorted(xs)
        n = len(s)
        mid = n // 2
        return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2

    def stdev(xs: List[float], m: float) -> Optional[float]:
        if len(xs) < 2:
            return None
        var = sum((x - m) ** 2 for x in xs) / (len(xs) - 1)
        return var ** 0.5

    subfield_rows: List[Dict[str, object]] = []
    for (subject_base, session, hemi), by_subfield in subfield_agg.items():
        fsid = f"{subject_base}_{session}"
        for subfield, vals in by_subfield.items():
            m = mean(vals)
            subfield_rows.append({
                "subject_id": subject_base, "session": session, "fsid": fsid, "hemisphere": hemi,
                "source": "hipsta", "region_raw": subfield, "region": subfield, "is_composite": False,
                "value": m, "etiv": "",
                "median_value": median(vals), "sd_value": stdev(vals, m), "n_vertices": len(vals),
            })

    return subfield_rows, grid_rows, missing


# ---------------------------------------------------------------------------
# --config mode: extract every measure a flex study.json declares
# ---------------------------------------------------------------------------
def _region_matches(row_region: str, roi_set) -> bool:
    if roi_set == "all" or roi_set is None:
        return True
    return row_region in roi_set


def filter_measure_rows(raw_rows: List[Dict[str, object]], region_filter: Optional[Dict], roi_set) -> List[Dict[str, object]]:
    out = raw_rows
    if region_filter:
        prefix = region_filter.get("prefix")
        include = region_filter.get("include")
        if prefix is not None:
            out = [r for r in out if str(r["region_raw"]).startswith(prefix)]
        if include is not None:
            out = [r for r in out if str(r["region_raw"]) in set(include)]
    out = [r for r in out if _region_matches(str(r["region"]), roi_set)]
    return out


def run_config_extraction(config_path: str) -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from config import ConfigError, dataset_path, load_config, resolve_path  # noqa: E402

    try:
        cfg = load_config(config_path)
    except ConfigError as e:
        print(f"INVALID CONFIG: {e}", file=sys.stderr)
        sys.exit(1)

    out_root = resolve_path(cfg, cfg["output"]["root"])
    tidy_dir = out_root / "tidy"
    logs_dir = out_root / "logs"

    sessions_include = cfg.get("sessions", {}).get("include")

    raw_cache: Dict[Tuple[str, str, str], List[Dict[str, object]]] = {}
    id_map_cache: Dict[str, Optional[Dict[str, str]]] = {}

    def id_map_for(ds_name: str, ds_cfg: Dict) -> Optional[Dict[str, str]]:
        if ds_name not in id_map_cache:
            map_path = ds_cfg.get("subject_id_map")
            if map_path:
                resolved = resolve_path(cfg, map_path)
                id_map_cache[ds_name] = load_subject_id_map(resolved)
                print(f"[{ds_name}] loaded subject_id_map: {len(id_map_cache[ds_name])} entries from {resolved}")
            else:
                id_map_cache[ds_name] = None
        return id_map_cache[ds_name]

    for m in cfg["measures"]:
        ds_name = m["dataset"]
        ds_cfg = cfg["datasets"][ds_name]
        ds_path = dataset_path(cfg, ds_name)
        dir_pattern = ds_cfg.get("dir_pattern", DEFAULT_DIR_PATTERN)
        source = m["source"]
        id_map = id_map_for(ds_name, ds_cfg)

        if not ds_path.is_dir():
            print(f"WARNING: dataset '{ds_name}' path does not exist: {ds_path} -- skipping measure '{m['name']}'",
                  file=sys.stderr)
            continue

        if source == "hipsta":
            subfield_rows, grid_rows, missing = extract_hipsta(ds_path, dir_pattern, sessions_include, id_map)
            subfield_cols = TIDY_COLS + ["median_value", "sd_value", "n_vertices"]
            write_tsv(tidy_dir / f"{m['name']}.tsv", subfield_cols, subfield_rows)
            msg = f"[{m['name']}] wrote {len(subfield_rows)} subfield rows"
            if m.get("pointwise"):
                grid_cols = ["subject_id", "session", "fsid", "hemisphere", "x", "y", "thickness", "hsf_label", "subfield"]
                write_tsv(tidy_dir / f"{m['name']}_grid.tsv", grid_cols, grid_rows)
                msg += f", {len(grid_rows)} grid rows"
            if missing:
                write_tsv(logs_dir / f"{m['name']}_missing.tsv",
                          ["subject_id", "session", "hemisphere", "reason"], missing)
            print(msg)
            continue

        metric = m.get("metric") or DEFAULT_METRIC.get(source, "")
        cache_key = (ds_name, source, metric)
        if cache_key not in raw_cache:
            raw_rows, missing = extract_single_source(ds_path, source, metric, dir_pattern, sessions_include, id_map)
            raw_cache[cache_key] = raw_rows
            if missing:
                write_tsv(logs_dir / f"{ds_name}_{source}_missing.tsv",
                          ["subject_id", "session", "source", "reason"], missing)
        raw_rows = raw_cache[cache_key]

        filtered = filter_measure_rows(raw_rows, m.get("region_filter"), m.get("roi_set", "all"))
        write_tsv(tidy_dir / f"{m['name']}.tsv", TIDY_COLS, filtered)
        print(f"[{m['name']}] wrote {len(filtered)} rows (source={source}, dataset={ds_name})")

    print(f"\nExtraction done. Tidy tables in {tidy_dir}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("subjects_dir", nargs="?", help="Path to the FreeSurfer SUBJECTS_DIR (omit when using --config)")
    p.add_argument("-o", "--outdir", help="Output directory (required unless --config)")
    p.add_argument("--config", default=None,
                   help="Path to a flex study.json -- extracts every configured measure and ignores "
                        "subjects_dir/-o/--source/--metric/--dir-pattern/--tidy-name/--missing-log")
    p.add_argument("--source", action="append",
                    choices=sorted(DEFAULT_METRIC.keys() | {"hippo-amygdala", "thalamus", "brainstem"}),
                    help="FreeSurfer output type to extract (repeatable, e.g. --source aparc.a2009s --source aseg)")
    p.add_argument("--metric", default=None,
                   help="Column to extract for aseg/aparc/aparc.a2009s sources [default: Volume_mm3 for aseg, ThickAvg for aparc*]")
    p.add_argument("--dir-pattern", default=DEFAULT_DIR_PATTERN,
                   help=f"Regex (named groups 'base','ses') matching per-timepoint directories [default: {DEFAULT_DIR_PATTERN}]")
    p.add_argument("--tidy-name", default="tidy.tsv", help="Output tidy TSV filename [default: %(default)s]")
    p.add_argument("--missing-log", default="missing.tsv", help="Filename (in --outdir) logging subjects/sessions with no data for a source")
    args = p.parse_args()

    if args.config:
        run_config_extraction(args.config)
        return

    if not args.subjects_dir or not args.outdir or not args.source:
        p.error("subjects_dir, -o/--outdir and --source are required unless --config is given")

    subjects_dir = Path(args.subjects_dir)
    if not subjects_dir.is_dir():
        print(f"Error: {subjects_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tidy_rows: List[Dict[str, object]] = []
    missing: List[Dict[str, str]] = []

    for source in args.source:
        metric = args.metric or DEFAULT_METRIC.get(source, "")
        rows, miss = extract_single_source(subjects_dir, source, metric, args.dir_pattern)
        tidy_rows.extend(rows)
        missing.extend(miss)

    if not tidy_rows:
        print("No data extracted for any requested source.", file=sys.stderr)
        sys.exit(1)

    write_tsv(outdir / args.tidy_name, TIDY_COLS, tidy_rows)
    print(f"Wrote {outdir / args.tidy_name} ({len(tidy_rows)} rows, sources={args.source})")

    if missing:
        write_tsv(outdir / args.missing_log, ["subject_id", "session", "source", "reason"], missing)
        print(f"Wrote {outdir / args.missing_log} ({len(missing)} missing subject/session/source combinations)", file=sys.stderr)


if __name__ == "__main__":
    main()
