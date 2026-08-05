# The flex pipeline: one JSON file per study

`scripts/analysis/*.R` (33 scripts) is the reviewed, frozen statistical battery
for **this study**. It stays exactly as-is — every one of its scripts, defaults,
and results is what got reviewed and reported.

The `flex` pipeline (`scripts/flex/`, `configs/flex/`) is a **second, general
implementation** of the same battery, driven entirely by one JSON definition
file per study, so a new study takes a new config file instead of a fork of 33
scripts. It was built by generalizing the existing scripts and verifying, ROI
by ROI, that every module reproduces the original scripts' numbers exactly
(see "Parity" below) — it is not a rewrite-and-hope.

```
DataLad remote (HPC)  →  ingest.py  →  extract_freesurfer.py  →  tidy/*.tsv  →  R test battery  →  results/
                         (datalad get)   (--config mode)                        (run_analysis.R)
```

For a hands-on guide with real commands, real output, and a full worked
example of porting a new study, see
[`FLEX_WALKTHROUGH.md`](FLEX_WALKTHROUGH.md) — this document covers the
concepts; that one covers doing it.

## Quick start

```bash
# 1. Validate a config without touching any data
python scripts/flex/run_study.py --config configs/flex/study.pk01.json --validate-only

# 2. See what would be fetched from the HPC, without fetching it
python scripts/flex/run_study.py --config configs/flex/study.pk01.json --stage ingest --dry-run

# 3. Run everything: ingest -> extract -> analyse
python scripts/flex/run_study.py --config configs/flex/study.pk01.json

# 4. Re-run just a subset of the battery
python scripts/flex/run_study.py --config configs/flex/study.pk01.json \
  --stage analyse --only lmm,moderator --measures hippo,cortex
```

Outputs land under `<output.root>/` from the config (e.g.
`results/flex/pk01/`): `tidy/` (the extracted long-format tables),
`results/<test>_<measure>/` (per-test-per-measure summaries and model dumps),
`logs/` (missing-data logs from ingest/extract), and a top-level
`run_manifest.json`.

## The study definition file

One JSON file (`configs/flex/study.<id>.json`, validated against
`configs/flex/study.schema.json`) describes a study end to end:

| Section | What it does |
|---|---|
| `study` | id/label |
| `datasets` | Where the raw FreeSurfer/hipsta output lives, and whether it's a DataLad sibling (`datalad: true` triggers the ingest stage) |
| `sessions` | Which sessions to include, and which one is baseline |
| `measures` | Named extractions: a FreeSurfer `source` + a `profile` + an ROI set. One measure = one tidy table. |
| `profiles` | Optional custom measure profiles (rarely needed — three built-ins cover almost everything, see below) |
| `design` | Participant file, column-name roles, group-value contrasts, factorial axes, covariates, moderators |
| `analyses` | Which test modules to run against which measures |
| `output` | Where results go, alpha, FDR scope |

See `configs/flex/study.TEMPLATE.json` for a minimal starting point and
`configs/flex/study.pk01.json` for a fully worked example (the current
study, reproduced from the actual `run_pipeline.sh`). Field-by-field
reference: `configs/flex/README.md`.

## The `profile`: the single biggest generalization

Seven of the original 33 scripts (`01_primary_lmm.R`, `10_amygdala_lmm.R`,
`11_basal_ganglia_lmm.R`, `24_thalamic_lmm.R`, `25_brainstem_lmm.R`,
`27_cortical_thickness_lmm.R`, `29_hippo_thickness_lmm.R`) are the *same*
mixed-effects model. They differ in exactly three things:

| | `transform` | `etiv_covariate` | `hemisphere` |
|---|---|---|---|
| **`volume`** (hippocampus, amygdala, thalamus, basal ganglia) | `log` | `true` | `pooled` |
| **`thickness`** (cortex, hipsta) | `identity` | `false` | `pooled` |
| **`volume_midline`** (brainstem) | `log` | `true` | `none` |

A `measure` declares which profile it uses; `test_lmm.R` assembles the
formula from the profile + the study's declared covariates. This is why
seven scripts collapse into one module (`scripts/flex/R/test_lmm.R`) — no
per-study code, just three profile flags per measure. The rationale for
each choice (why volume needs eTIV, why thickness doesn't, why raw scale
for thickness) is carried over verbatim from `01_primary_lmm.R`'s and
`29_hippo_thickness_lmm.R`'s own headers — read those if you need the
underlying reasoning, not just the recipe.

Custom profiles can be declared per-study under `profiles` if a new
study's measure genuinely doesn't fit any of the three built-ins.

## Test modules implemented

| `analyses[].test` | Replaces | Notes |
|---|---|---|
| `lmm` | 01, 10, 11, 24, 25, 27, 29 | Driven entirely by `profile` |
| `change_ancova` | 05 | Baseline-corrected ANCOVA on change scores |
| `multivariate` | 06, 28 | MANOVA + PCA; set `region_groups` on the measure to run once per cluster (28's design) |
| `moderator` | 08 | Age/sex + every `design.moderators.baseline_cols` column, threshold-gated LOSO |
| `factorial` | 13, 30 | Crosses every axis in `design.factorial` (or `spec.factorial_axes`) |
| `loso` | 07 | LOSO for the change-score model and the MANOVA test |
| `rm_anova` | 03, 22, 23 | `engine: "afex"` (classical rm-ANOVA, configurable `between`) or `"lmm_interaction"` (N-way crossing + optional continuous moderator) |

**Not yet ported** (use the original `scripts/analysis/*.R` script directly
for these): `02` (hierarchical brms), `04` (exploratory RF), `09`
(continuous-time LMM — a `time: "continuous"` option on `lmm` is the
natural home for this), `12` (responder heterogeneity), `14` (coupled
change), `15`/`18` (pooled-studies / proportion-of-whole variants of
`change_ancova`), `16` (general time pattern), `17` (trajectory
clustering), `19`/`20`/`33` (moderator variants), `21` (latent growth
curve), `26` (jmv-based full-covariate GLM — the `jmv` package isn't in
this project's `renv.lock`), `31`/`32` (grid-point/pointwise hipsta
analyses — a genuinely different feature space from every other module
here). Running an unported test via `run_analysis.R` prints a `SKIP`
message naming the module and pointing at the frozen script instead of
failing the whole run.

## Parity: how "generalized" was verified, not assumed

Every module above was checked against its original script on the actual
study-129 data, not just read and trusted:

- **Extraction**: `hippo_thick` matches `extract_hipsta_thickness.py`'s
  output byte-for-byte (2537/2537 rows, max diff 0.0). `thalamic` and
  `brainstem` match the "final" merged reference tables exactly
  (13000/13000 and 1250/1250 rows). `hippo`/`amygdala`/`cortex` match
  exactly except for 2 subjects whose upstream FreeSurfer segmentation was
  reprocessed *after* the reference tables were built — a real data
  update, not an extraction bug (confirmed via file mtimes).
- **Statistics**: every module's summary output was diffed against its
  original script's output on the same extracted data. F-statistics,
  p-values, and estimates agree to ~1e-8 or exactly (0.0 diff) in every
  case checked — `test_lmm.R` (01/27/29, three different profiles),
  `test_change_ancova.R` (05), `test_multivariate.R` (06), `test_moderator.R`
  (08, 56/56 rows), `test_factorial.R` (13), `test_loso.R` (07, 889 + 127
  LOSO combinations), `test_rm_anova.R` (03, 22, 23 — three different
  engine/factor configurations).
- **One deliberate, documented divergence**: `05`/`06`/`13` hardcode
  raw-scale `age_z`, while `01`/`27`/`29`/`33` use rank-transformed
  `age_z` — a pre-existing inconsistency across the original 33 scripts,
  not a considered design choice. This engine applies `design.covariates`
  uniformly across every test (one declared transform per study), which is
  the more defensible default. Parity checks against the raw-z scripts
  used a locally rank-z-patched copy for a fair, isolated comparison of
  everything *except* that one difference; every such check still matched
  to ~1e-8.

If you port a new module, hold it to the same bar: diff its output against
the original script on real data, not just against your own expectations
of what the formula should produce.

## Porting a new study

1. Copy `configs/flex/study.TEMPLATE.json` to `configs/flex/study.<id>.json`.
2. Point `datasets` at wherever the new study's FreeSurfer/hipsta output
   lives. Set `datalad: true` if it's a DataLad sibling you need to `datalad
   get` from (e.g. an HPC-processed dataset) — the ingest stage handles the
   rest.
3. If the raw derivatives directory uses different subject IDs than your
   participants file (e.g. a pseudonymized short ID vs. the study's own
   coding), build a mapping with `scripts/build_hipsta_subject_id_map.py`
   and point `datasets.<name>.subject_id_map` at it. Don't guess mappings —
   that script only emits unambiguous matches and logs the rest for a human
   to resolve.
4. List the `measures` you need, each with a `source` and a `profile`
   (pick one of the three built-ins unless your measure genuinely needs a
   different transform/eTIV/hemisphere combination).
5. Fill in `design`: your participants file's actual column names under
   `roles`, the group values that make up your intervention/control
   contrast (or omit `contrasts` entirely if the study isn't an
   intervention design — `lmm`/`change_ancova`/etc. that need a contrast
   will say so explicitly rather than silently misbehaving), any factorial
   axes, and covariates.
6. Pick `analyses` from the table above.
7. `python scripts/flex/run_study.py --config configs/flex/study.<id>.json
   --validate-only`, fix whatever it flags, then run for real.

No file under `scripts/flex/` should ever need editing to onboard a new
study that fits the existing measures/tests — if you find yourself editing
`test_*.R` or `extract_freesurfer.py` to accommodate a new study's
specifics, that specific thing probably belongs in the JSON instead.
