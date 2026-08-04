#!/usr/bin/env python3
"""
Extract longitudinal hippocampal subfield THICKNESS (hipsta shape/thickness
analysis, https://github.com/Deep-MI/Hipsta) and emit two tidy long TSVs:

  - a grid-point-level tidy TSV (one row per subject/session/hemisphere/grid
    vertex) for the point-wise scripts under scripts/analysis/ (31, 32)
  - a subfield-level tidy TSV (grid points aggregated by subfield label, one
    row per subject/session/hemisphere/subfield) for the ROI-style scripts
    (29, 30) -- column names ("region", "value") deliberately match
    aparc_a2009s_tidy_2tp_fs82.tsv's contract (extract_freesurfer.py) rather
    than hippo_subfields_tidy.tsv's ("subfield", "volume"), so that
    27_cortical_thickness_lmm.R's formula/FDR logic could, in principle, be
    pointed at either file unmodified.

Input layout expected (produced by scripts/run_hipsta_segmentation.sh):

  <HIPSTA_DIR>/sub-<ID>_ses-<N>.long.sub-<ID>/hipsta-status_lh.txt
  <HIPSTA_DIR>/sub-<ID>_ses-<N>.long.sub-<ID>/hipsta-status_rh.txt
  <HIPSTA_DIR>/sub-<ID>_ses-<N>.long.sub-<ID>/thickness/lh.grid-segments-z.csv
  <HIPSTA_DIR>/sub-<ID>_ses-<N>.long.sub-<ID>/thickness/lh.mid-surface.hsf.csv
  (and the rh.* equivalents)

grid-segments-z.csv is a 41 (x, medial->lateral) x 21 (y, posterior->anterior)
table of thickness values (mm), with a header row/column of "x<i>"/"y<j>"
labels. mid-surface.hsf.csv is an unheaded (x_idx, y_idx, hsf_label) table on
the SAME grid, giving each grid point's FreeSurfer subfield label. The two are
joined on (x_idx, y_idx).

Sessions/hemispheres without an OK status stamp, or with unreadable/missing
thickness files, are logged to --missing-log rather than silently dropped.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SESSION_DIR_RE = re.compile(r"^(?P<base>sub-[^/_]+)_(?P<ses>ses-[^/.]+)\.long\.sub-[^/]+$")

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
    have no label (hipsta writes an empty field there, seen e.g. for the
    x=0 column on some hemispheres) -- those points are simply omitted, and
    end up excluded by the (x,y) intersection with grid-segments-z.csv.
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


def find_sessions(hipsta_dir: Path) -> List[Tuple[str, str, Path]]:
    """Return (subject_base, session, dir) for every sub-X_ses-Y.long.sub-X directory."""
    out = []
    for d in sorted(hipsta_dir.iterdir()):
        if not d.is_dir():
            continue
        m = SESSION_DIR_RE.match(d.name)
        if not m:
            continue
        out.append((m.group("base"), m.group("ses"), d))
    return out


def hemi_status_ok(session_dir: Path, hemi: str) -> Tuple[bool, str]:
    """Check hipsta-status_<hemi>.txt; returns (ok, reason_if_not)."""
    stamp = session_dir / f"hipsta-status_{hemi}.txt"
    if not stamp.is_file():
        return False, "hipsta-status file not found (segmentation not run or in progress)"
    with open(stamp) as f:
        first_line = f.readline().strip()
    if first_line == "status: OK":
        return True, ""
    return False, f"hipsta-status reports failure ({first_line or 'empty status file'})"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("hipsta_dir", help="Path to the hipsta output directory (run_hipsta_segmentation.sh's OUTPUT_BASE)")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument(
        "--grid-tidy-name",
        default="hippo_thickness_grid_tidy.tsv",
        help="Filename for the grid-point-level tidy TSV [default: %(default)s]",
    )
    p.add_argument(
        "--subfield-tidy-name",
        default="hippo_thickness_tidy.tsv",
        help="Filename for the subfield-aggregated tidy TSV [default: %(default)s]",
    )
    p.add_argument(
        "--missing-log",
        default="missing_hipsta_subjects.tsv",
        help="Filename (in --outdir) logging subjects/sessions/hemispheres missing thickness output [default: %(default)s]",
    )
    args = p.parse_args()

    hipsta_dir = Path(args.hipsta_dir)
    if not hipsta_dir.is_dir():
        print(f"Error: {hipsta_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sessions = find_sessions(hipsta_dir)
    if not sessions:
        print(f"No sub-*_ses-*.long.sub-* directories found under {hipsta_dir}", file=sys.stderr)
        sys.exit(1)

    grid_rows: List[Dict[str, object]] = []
    subfield_agg: Dict[Tuple[str, str, str], Dict[str, List[float]]] = {}
    missing: List[Dict[str, str]] = []

    for subject_base, session, session_dir in sessions:
        fsid = f"{subject_base}_{session}"

        for hemi in ("lh", "rh"):
            ok, reason = hemi_status_ok(session_dir, hemi)
            if not ok:
                missing.append(
                    {"subject_id": subject_base, "session": session, "hemisphere": hemi, "reason": reason}
                )
                continue

            grid_file = session_dir / "thickness" / f"{hemi}.grid-segments-z.csv"
            hsf_file = session_dir / "thickness" / f"{hemi}.mid-surface.hsf.csv"
            if not grid_file.is_file() or not hsf_file.is_file():
                missing.append(
                    {
                        "subject_id": subject_base,
                        "session": session,
                        "hemisphere": hemi,
                        "reason": "status OK but thickness/grid-segments-z.csv or mid-surface.hsf.csv missing",
                    }
                )
                continue

            try:
                thickness = parse_grid_csv(grid_file)
                hsf = parse_hsf_csv(hsf_file)
            except (OSError, ValueError) as e:
                missing.append(
                    {
                        "subject_id": subject_base,
                        "session": session,
                        "hemisphere": hemi,
                        "reason": f"failed to parse thickness/hsf files: {e}",
                    }
                )
                continue

            common_idx = sorted(set(thickness) & set(hsf))
            if not common_idx:
                missing.append(
                    {
                        "subject_id": subject_base,
                        "session": session,
                        "hemisphere": hemi,
                        "reason": "grid-segments-z.csv and mid-surface.hsf.csv share no (x,y) index",
                    }
                )
                continue

            for x, y in common_idx:
                subfield = hsf_label_name(hsf[(x, y)])
                val = thickness[(x, y)]
                grid_rows.append(
                    {
                        "subject_id": subject_base,
                        "session": session,
                        "fsid": fsid,
                        "hemisphere": hemi,
                        "x": x,
                        "y": y,
                        "thickness": val,
                        "hsf_label": hsf[(x, y)],
                        "subfield": subfield,
                    }
                )
                key = (subject_base, session, hemi)
                subfield_agg.setdefault(key, {}).setdefault(subfield, []).append(val)

    if not grid_rows:
        print("No complete subject/session/hemisphere thickness outputs found.", file=sys.stderr)
        sys.exit(1)

    # --- write grid-point-level tidy TSV ---
    grid_path = outdir / args.grid_tidy_name
    grid_cols = ["subject_id", "session", "fsid", "hemisphere", "x", "y", "thickness", "hsf_label", "subfield"]
    with open(grid_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=grid_cols, delimiter="\t")
        w.writeheader()
        w.writerows(grid_rows)
    print(f"Wrote {grid_path} ({len(grid_rows)} rows)")

    # --- write subfield-aggregated tidy TSV ---
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
        return var**0.5

    subfield_rows: List[Dict[str, object]] = []
    for (subject_base, session, hemi), by_subfield in subfield_agg.items():
        fsid = f"{subject_base}_{session}"
        for subfield, vals in by_subfield.items():
            m = mean(vals)
            subfield_rows.append(
                {
                    "subject_id": subject_base,
                    "session": session,
                    "fsid": fsid,
                    "hemisphere": hemi,
                    "region": subfield,
                    "value": m,
                    "median_value": median(vals),
                    "sd_value": stdev(vals, m),
                    "n_vertices": len(vals),
                }
            )

    subfield_path = outdir / args.subfield_tidy_name
    subfield_cols = [
        "subject_id",
        "session",
        "fsid",
        "hemisphere",
        "region",
        "value",
        "median_value",
        "sd_value",
        "n_vertices",
    ]
    with open(subfield_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=subfield_cols, delimiter="\t")
        w.writeheader()
        w.writerows(subfield_rows)
    print(f"Wrote {subfield_path} ({len(subfield_rows)} rows)")

    # --- write missing log ---
    if missing:
        missing_path = outdir / args.missing_log
        with open(missing_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["subject_id", "session", "hemisphere", "reason"], delimiter="\t")
            w.writeheader()
            w.writerows(missing)
        print(
            f"Wrote {missing_path} ({len(missing)} subject/session/hemisphere combinations missing thickness output)",
            file=sys.stderr,
        )
    else:
        print("No missing subject/session/hemisphere combinations.", file=sys.stderr)


if __name__ == "__main__":
    main()
