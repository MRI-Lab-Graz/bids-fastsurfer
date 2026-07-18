#!/usr/bin/env python3
"""
Extract longitudinal psychological moderator variables from
data/survey/long_en/combined_survey.csv (study 129 -- "Brainhearthlon:
Running away from depression") into a tidy, subject x session TSV.

Handles three real data-quality issues found in the source file:

1. ID/subject_id swap: for a subset of rows, the `subject_id` column holds
   a stray small integer while the properly-formatted "sub-129XXXX" id is
   in the `ID` column instead (and vice versa for other rows). Fixed by
   taking whichever of the two columns matches the expected pattern.

2. Corrupted time fields: PSQI01 (bed time), PSQI03 (wake time), and PSQI04
   (sleep duration) are stored as HH:MM:SS strings with a literal "NA"
   erroneously inserted into the minutes/seconds fields (e.g.
   "22:NA30:NA00" instead of "22:30:00"). Verified fully decodable by
   stripping the literal "NA" substring; every resulting value across the
   dataset parses to a valid, plausible time/duration.

3. Missing PSQI sub-items: PSQI05-01 through PSQI05-10 (the sleep
   disturbance checklist) are entirely empty for every subject/session --
   not corrupted, never collected. This makes the standard 7-component
   PSQI Global Score impossible to compute validly (missing the Sleep
   Disturbances component and part of Sleep Latency). This script instead
   computes a clearly-labelled PARTIAL score from the 5 available
   components (Quality, Duration, Efficiency, Medication, Daytime
   Dysfunction) -- this is NOT the validated PSQI Global Score and must
   not be interpreted against the standard >5 "poor sleep" cutoff.

Scoring assumptions (ASSUMED STANDARD instrument versions -- the source
codebooks contain no item-level labels or scoring documentation, so these
follow the canonical published item order/reverse-scoring for each
instrument; verify against the actual questionnaire used before treating
results as final):
  - STAI (20 items, 1-4 scale): trait-anxiety form, reverse-scored items
    1, 6, 7, 10, 13, 16, 19 (Spielberger/Laux et al. German STAI-Trait,
    Form X-2), summed to a 20-80 total.
  - FSozU (54 items, 1-5 scale): standard F-SozU-54 has no reverse items;
    total = mean of all valid items.
  - ADS / PANAS-trait / TSDZ: already have precomputed scores in the
    source file (ADS_score, panastrait_PA, panastrait_NA,
    tsdz_total_score) -- passed through unchanged, no scoring needed.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Optional

SUBJECT_ID_RE = re.compile(r"^sub-129\d+$")

STAI_REVERSE_ITEMS = {1, 6, 7, 10, 13, 16, 19}


def resolve_subject_id(row: dict) -> Optional[str]:
    sid = (row.get("subject_id") or "").strip()
    idcol = (row.get("ID") or "").strip()
    if SUBJECT_ID_RE.match(sid):
        return sid
    if SUBJECT_ID_RE.match(idcol):
        return idcol
    return None


def decode_time_field(raw: str) -> Optional[float]:
    """Strip the literal 'NA' corruption and parse HH:MM:SS -> hours (float)."""
    if not raw:
        return None
    cleaned = raw.replace("NA", "")
    parts = cleaned.split(":")
    if len(parts) != 3 or not all(p.isdigit() for p in parts):
        return None
    h, m, s = (int(p) for p in parts)
    if not (0 <= h < 24 and 0 <= m < 60 and 0 <= s < 60):
        return None
    return h + m / 60 + s / 3600


def to_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def score_stai(row: dict) -> Optional[float]:
    total = 0
    n = 0
    for i in range(1, 21):
        v = to_float(row.get(f"STAI-{i:02d}", ""))
        if v is None:
            continue
        total += (5 - v) if i in STAI_REVERSE_ITEMS else v
        n += 1
    if n < 20:  # require complete data -- partial STAI sums aren't comparable
        return None
    return total


def score_fsozu(row: dict) -> Optional[float]:
    vals = []
    for i in range(1, 55):
        v = to_float(row.get(f"FSozU-{i:02d}", ""))
        if v is not None:
            vals.append(v)
    if len(vals) < 54:
        return None
    return sum(vals) / len(vals)


def recode_c1_quality(psqi06: Optional[float]) -> Optional[float]:
    return psqi06  # PSQI06 is already the 0-3 "overall sleep quality" rating


def recode_c3_duration(hours: Optional[float]) -> Optional[float]:
    if hours is None:
        return None
    if hours > 7:
        return 0
    if hours >= 6:
        return 1
    if hours >= 5:
        return 2
    return 3


def recode_c4_efficiency(bed_h: Optional[float], wake_h: Optional[float], sleep_h: Optional[float]) -> Optional[float]:
    if bed_h is None or wake_h is None or sleep_h is None:
        return None
    time_in_bed = (wake_h - bed_h) % 24
    if time_in_bed <= 0:
        return None
    efficiency = 100 * sleep_h / time_in_bed
    if efficiency >= 85:
        return 0
    if efficiency >= 75:
        return 1
    if efficiency >= 65:
        return 2
    return 3


def recode_c6_medication(psqi07: Optional[float]) -> Optional[float]:
    return psqi07  # already 0-3 frequency scale


def recode_c7_daytime(psqi08: Optional[float], psqi09: Optional[float]) -> Optional[float]:
    if psqi08 is None or psqi09 is None:
        return None
    combined = psqi08 + psqi09
    if combined == 0:
        return 0
    if combined <= 2:
        return 1
    if combined <= 4:
        return 2
    return 3


def score_psqi_partial(row: dict) -> tuple[Optional[float], int]:
    bed_h = decode_time_field(row.get("PSQI01", ""))
    wake_h = decode_time_field(row.get("PSQI03", ""))
    sleep_h = decode_time_field(row.get("PSQI04", ""))

    c1 = recode_c1_quality(to_float(row.get("PSQI06", "")))
    c3 = recode_c3_duration(sleep_h)
    c4 = recode_c4_efficiency(bed_h, wake_h, sleep_h)
    c6 = recode_c6_medication(to_float(row.get("PSQI07", "")))
    c7 = recode_c7_daytime(to_float(row.get("PSQI08", "")), to_float(row.get("PSQI09", "")))

    components = [c1, c3, c4, c6, c7]
    n_available = sum(c is not None for c in components)
    if n_available == 0:
        return None, 0
    return sum(c for c in components if c is not None), n_available


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("survey_csv", help="Path to combined_survey.csv under data/survey/long_en/")
    p.add_argument("-o", "--output", required=True, help="Output TSV path")
    args = p.parse_args()

    survey_path = Path(args.survey_csv)
    if not survey_path.is_file():
        print(f"Error: {survey_path} not found.", file=sys.stderr)
        sys.exit(1)

    with open(survey_path) as f:
        rows = list(csv.DictReader(f))

    out_rows = []
    unresolved = []
    for row in rows:
        subject_id = resolve_subject_id(row)
        session = row.get("session", "").strip()
        if not subject_id or not session:
            unresolved.append(row.get("participant_id", "?"))
            continue

        psqi_partial, psqi_n_components = score_psqi_partial(row)
        fsozu_mean = score_fsozu(row)
        stai_total = score_stai(row)

        out_rows.append({
            "subject_id": subject_id,
            "session": session,
            "ads_score": row.get("ADS_score", "") or "",
            "ads_psy_symptom": row.get("ADS_PsySymptom", "") or "",
            "ads_soma_symptom": row.get("ADS_SomaSymptom", "") or "",
            "panas_pa": row.get("panastrait_PA", "") or "",
            "panas_na": row.get("panastrait_NA", "") or "",
            "tsdz_total": row.get("tsdz_total_score", "") or "",
            "fsozu_mean": fsozu_mean if fsozu_mean is not None else "",
            "stai_total": stai_total if stai_total is not None else "",
            "psqi_partial_score": psqi_partial if psqi_partial is not None else "",
            "psqi_n_components_of_5": psqi_n_components,
        })

    out_rows.sort(key=lambda r: (r["subject_id"], r["session"]))
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["subject_id", "session", "ads_score", "ads_psy_symptom", "ads_soma_symptom",
                  "panas_pa", "panas_na", "tsdz_total", "fsozu_mean", "stai_total",
                  "psqi_partial_score", "psqi_n_components_of_5"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {out_path} ({len(out_rows)} subject/session rows)")
    if unresolved:
        print(f"WARNING: {len(unresolved)} rows had no resolvable subject_id or session and were dropped "
              f"(participant_ids: {unresolved[:10]}{'...' if len(unresolved) > 10 else ''})", file=sys.stderr)
    print("NOTE: psqi_partial_score is NOT the validated PSQI Global Score -- it covers only "
          "5 of 7 standard components (Sleep Disturbances and part of Sleep Latency are unavailable, "
          "the underlying sub-items were never collected). Do not compare against the standard >5 cutoff.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
