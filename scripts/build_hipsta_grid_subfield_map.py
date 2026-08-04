#!/usr/bin/env python3
"""
Build a static grid-point -> subfield mapping for hipsta's fixed 41x21 thickness
grid, for use as the --clusters argument of
scripts/analysis/31_hippo_thickness_cluster_multivariate.R.

Why this needs to exist at all: hipsta's cube parametrization puts grid index
(x, y) at the same relative anatomical position for every subject (that's the
entire point of the parametrization -- no separate surface registration step
is needed, unlike e.g. Destrieux/aparc.a2009s regions on native cortical
surfaces). But the exact subfield BOUNDARY within that shared grid is still
estimated per subject (from each subject's own segmentation), so it wobbles
by a grid cell or two from one subject to the next. For a per-subfield
multivariate analysis (31, mirroring 28_cortical_cluster_multivariate.R's
per-cluster design) the grid columns that make up "CA1" must be the SAME
861-point-indexed columns for every subject -- so this script takes a
majority vote, per (hemisphere, x, y), across the whole cohort's
hippo_thickness_grid_tidy.tsv (from extract_hipsta_thickness.py), and writes
that as a fixed lookup table, analogous to how configs/destrieux_functional_clusters.tsv
is a fixed region->cluster table rather than something recomputed per subject.

Grid points where the winning label's agreement is below --min-agreement are
dropped (--min-agreement default 0.6) rather than assigned to a subfield they
don't reliably belong to -- these sit near a subfield boundary that itself
moves across subjects, so no single subfield assignment would be meaningful
there. This means 31's per-subfield analysis is deliberately conservative:
it excludes boundary-adjacent grid points rather than mislabel them.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("grid_tidy", help="Path to hippo_thickness_grid_tidy.tsv (from extract_hipsta_thickness.py)")
    p.add_argument("-o", "--output", required=True, help="Output path for the grid_id -> subfield TSV")
    p.add_argument(
        "--min-agreement",
        type=float,
        default=0.6,
        help="Minimum fraction of subjects agreeing on a grid point's subfield label to keep it [default: %(default)s]",
    )
    p.add_argument(
        "--exclude-unlabeled",
        action="store_true",
        default=True,
        help="Drop grid points whose majority label is 'unlabeled' (grid-edge points outside any subfield) [default: on]",
    )
    args = p.parse_args()

    grid_path = Path(args.grid_tidy)
    if not grid_path.is_file():
        raise SystemExit(f"grid tidy file not found: {grid_path}")

    # votes[(hemisphere, x, y)][subfield] = count of (subject, session) observations
    votes: dict[tuple[str, int, int], Counter] = defaultdict(Counter)
    n_obs: dict[tuple[str, int, int], int] = defaultdict(int)

    with open(grid_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["hemisphere"], int(row["x"]), int(row["y"]))
            votes[key][row["subfield"]] += 1
            n_obs[key] += 1

    out_rows = []
    n_dropped_agreement = 0
    n_dropped_unlabeled = 0
    for (hemi, x, y), counter in votes.items():
        winner, count = counter.most_common(1)[0]
        agreement = count / n_obs[(hemi, x, y)]
        if args.exclude_unlabeled and winner == "unlabeled":
            n_dropped_unlabeled += 1
            continue
        if agreement < args.min_agreement:
            n_dropped_agreement += 1
            continue
        out_rows.append(
            {
                "grid_id": f"{hemi}_x{x}_y{y}",
                "hemisphere": hemi,
                "x": x,
                "y": y,
                "subfield": winner,
                "n_observations": n_obs[(hemi, x, y)],
                "agreement_pct": round(100 * agreement, 1),
            }
        )

    out_rows.sort(key=lambda r: (r["hemisphere"], r["x"], r["y"]))

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["grid_id", "hemisphere", "x", "y", "subfield", "n_observations", "agreement_pct"], delimiter="\t"
        )
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {outpath} ({len(out_rows)} grid points)")
    print(f"  dropped {n_dropped_unlabeled} points whose majority label was 'unlabeled'")
    print(f"  dropped {n_dropped_agreement} points below --min-agreement={args.min_agreement} (boundary-adjacent)")
    subfield_counts = Counter(r["subfield"] for r in out_rows)
    for hemi in sorted({r["hemisphere"] for r in out_rows}):
        counts = Counter(r["subfield"] for r in out_rows if r["hemisphere"] == hemi)
        print(f"  {hemi}: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
