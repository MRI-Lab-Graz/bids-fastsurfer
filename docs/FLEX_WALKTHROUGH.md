# flex pipeline: hands-on walkthrough

This is a worked-example companion to [`FLEX_PIPELINE.md`](FLEX_PIPELINE.md)
(which explains the concepts — profiles, test modules, parity verification).
Everything below was actually run against real data; commands and output
are copied verbatim, not illustrative.

---

## 1. Running an existing study end to end

`configs/flex/study.pk01.json` is the current study, fully configured.

### 1.1 Validate first, always

```bash
python scripts/flex/run_study.py --config configs/flex/study.pk01.json --validate-only
```
```
Config OK: configs/flex/study.pk01.json (pk01) -- 6 measures, 9 analyses
```

This catches typos before touching any data. Two examples of what it catches
(both tested against this exact config):

**Missing required field** (deleted `design.roles.subject`):
```
INVALID CONFIG: study config failed schema validation (1 error(s)):
  /design/roles: 'subject' is a required property
```

**Semantic error** (added `"Control"` to both `contrasts.intervention` and `contrasts.control`):
```
INVALID CONFIG: /design/contrasts: group level(s) ['Control'] appear in
both 'intervention' and 'control'
```

The error always names the JSON pointer path (`/design/roles`,
`/design/contrasts`) so you know exactly where to look.

### 1.2 See what would be fetched, without fetching

```bash
python scripts/flex/run_study.py --config configs/flex/study.pk01.json --stage ingest --dry-run
```
```
[freesurfer] sources=['aparc.a2009s', 'brainstem', 'hippo-amygdala', 'thalamus']: 2511 file(s) need fetching
  would fetch: sub-001_ses-1.long.sub-001/mri/brainstemSsLabels.long.volumes.txt
  would fetch: sub-001_ses-1.long.sub-001/mri/lh.amygNucVolumes.long.txt
  ...
DRY RUN: 2511 file(s) would be fetched across 1 dataset(s).
```

Only `stats/*.stats` and `mri/*.txt` tables are ever fetched — never the
`.mgz` volumes. This ran in under 2 seconds against a 711-subject/session
DataLad dataset.

### 1.3 Run for real

```bash
python scripts/flex/run_study.py --config configs/flex/study.pk01.json
```

This runs `ingest` → `extract` → `analyse` in sequence. On the full pk01
config (6 measures, 9 analyses) this took **4m20s** and produced:

```
Done: 16 test/measure run(s) completed, 0 skipped/failed.
```

You'll see `datalad get` progress during ingest, `[measure] wrote N rows`
lines during extract, and `=== test / measure -> outdir ===` lines during
analyse. `boundary (singular) fit` messages from `lme4` are normal —
they're the trigger for the automatic random-slope → random-intercept
fallback, not an error.

### 1.4 Run just a slice of it

```bash
# Only the LMM battery, only two measures
python scripts/flex/run_study.py --config configs/flex/study.pk01.json \
  --stage analyse --only lmm --measures hippo,cortex

# Everything except the slow LOSO checks
python scripts/flex/run_study.py --config configs/flex/study.pk01.json \
  --stage analyse --only lmm,change_ancova,multivariate,moderator,factorial

# Just re-extract (e.g. after editing an ROI set), skip ingest/analyse
python scripts/flex/run_study.py --config configs/flex/study.pk01.json --stage extract
```

`--stage` picks one of `ingest`/`extract`/`analyse`; `--only` and
`--measures` (analyse-stage only) further restrict which
`analyses[]` entries run. Each stage is independently safe to re-run.

---

## 2. Reading the output

Everything lands under `output.root` from the config (`results/flex/pk01/`
here):

```
results/flex/pk01/
├── tidy/                          # one TSV per measure (extract stage)
│   ├── hippo.tsv
│   ├── amygdala.tsv
│   ├── cortex.tsv
│   ├── hippo_thick.tsv
│   ├── hippo_thick_grid.tsv       # companion grid-point file (pointwise: true)
│   └── ...
├── logs/                          # missing-data logs from ingest/extract
├── results/
│   ├── lmm_hippo/
│   │   ├── summary.csv            # one row per ROI, FDR-corrected
│   │   ├── models/
│   │   │   ├── CA1_model.txt                          # full lmer() summary
│   │   │   ├── CA1_intervention_vs_control_by_time.csv # emmeans contrasts
│   │   │   └── CA1_style_contrasts_by_time.csv         # secondary "style" contrast
│   │   ├── all_roi_models.rds     # every fitted model object, for follow-up in R
│   │   └── run_manifest.json      # what produced this
│   ├── change_ancova_hippo/
│   ├── multivariate_hippo/  multivariate_cortex/
│   ├── moderator_hippo/            # + loso_moderator_hits.csv for screened hits
│   ├── factorial_hippo/  factorial_hippo_thick/
│   ├── loso_hippo/                 # loso_change_score_by_roi.csv, loso_manova.csv, ...
│   ├── rm_anova_afex_hippo/
│   ├── rm_anova_colleague_replication_hippo/
│   └── rm_anova_colleague_cesd_hippo/
└── run_manifest.json               # top-level: which stages ran, when
```

Every `results/<test>_<measure>/run_manifest.json` records the exact config
file, its md5 hash, the R version, and the row/subject counts that
produced that output — e.g.:

```json
{
  "test": "lmm",
  "config_path": "/data/local/software/fast-surfer-workshop/configs/flex/study.pk01.json",
  "config_sha1": "86473126f47393e8dfd618aea4437cfe",
  "run_time": "2026-08-05T07:39:58+0200",
  "r_version": "R version 4.6.1 (2026-06-24)",
  "measure": "hippo",
  "n_subjects": 127,
  "n_rois": 7
}
```

`summary.csv` for `lmm` looks like this (real output, `lmm_hippo`):

```
roi,re_structure,n_obs,interaction_F,interaction_p,intervention_effect_last_timepoint,intervention_effect_last_timepoint_p,interaction_p_fdr,significant_interaction_p,...
CA1,random_intercept_only,508,0.0049,0.9442,-0.0036,0.8152,0.9839,FALSE,...
```

Every module's `summary.csv` follows the same shape: one row per ROI, raw
p-values plus `<col>_fdr` and `significant_<col>` columns added
automatically by `flex_write_summary()` (`core_report.R`).

---

## 3. Porting a brand-new study, from scratch

This is the actual sequence used to onboard study 134 (`configs/flex/study.af19.json`)
as a portability test — a genuinely different design from pk01: 3 sessions
instead of 2, a 3-arm `ballet`/`contemporary`/`control` design instead of
5-arm, different subject IDs, a different moderators file. **No code was
changed** to support it.

### 3.1 Look at what you actually have

```bash
head -3 /data/local/134_AF19/participants.tsv
```
```
subject_id      group           age  sex
sub-134001      control         27   F
sub-134002      control         22   M
```
```bash
ls /data/local/134_AF19/derivatives/fastsurfer/ | grep -E "^sub-[0-9]+_ses-[0-9]+\.long\."| head -3
```
```
sub-134001_ses-1.long.sub-134001
sub-134001_ses-2.long.sub-134001
sub-134001_ses-3.long.sub-134001
```

Subject IDs in `participants.tsv` (`sub-134001`) match the FreeSurfer
directory names directly here — no `subject_id_map` needed. (pk01 needed
one, because its DataLad-sourced directories use short pseudonymized IDs
like `sub-001` while `participants.tsv` uses `sub-1291003` — see
`configs/hipsta_subject_id_map.tsv` and
`scripts/build_hipsta_subject_id_map.py`. If your new study has the same
mismatch, build a map the same way rather than guessing at a rename rule.)

### 3.2 Write the config

```json
{
  "study": { "id": "af19", "label": "134 AF19 dance intervention (ballet+contemporary vs control), 3-timepoint" },

  "datasets": {
    "freesurfer": {
      "path": "/data/local/134_AF19/derivatives/fastsurfer",
      "datalad": false,
      "dir_pattern": "^(?P<base>sub-[^/_]+)_(?P<ses>ses-[^/]+)\\.long\\.(?P=base)$"
    }
  },

  "sessions": { "include": ["ses-1", "ses-2", "ses-3"], "baseline": "ses-1" },

  "measures": [
    {
      "name": "hippo", "dataset": "freesurfer", "source": "hippo-amygdala", "profile": "volume",
      "region_filter": { "prefix": "hippo:" },
      "roi_set": ["Whole_hippocampus", "GC-ML-DG", "CA1", "CA3", "CA4", "subiculum", "molecular_layer_HP"]
    }
  ],

  "design": {
    "participants": "/data/local/134_AF19/participants.tsv",
    "roles": { "subject": "subject_id", "group": "group", "age": "age", "sex": "sex" },
    "contrasts": { "intervention": ["ballet", "contemporary"], "control": ["control"] },
    "covariates": [
      { "col": "age", "type": "numeric", "transform": "rank_z", "as": "age_z" },
      { "col": "sex", "type": "factor" }
    ],
    "moderators": {
      "file": "/data/local/134_AF19/moderators.tsv",
      "baseline_cols": ["psqi_global", "pss_total", "ads_score"]
    }
  },

  "analyses": [
    { "test": "lmm", "measures": ["hippo"] },
    { "test": "change_ancova", "measures": ["hippo"] }
  ],

  "output": { "root": "results/flex/af19", "alpha": 0.05, "fdr_within": "measure" }
}
```

Every value here came directly from step 3.1 — no invented conventions.
`datalad: false` because this data was already materialized locally
(check with `git annex find --not --in=here` from inside the dataset dir;
see `docs/FLEX_PIPELINE.md`'s ingest section for what that's checking).

### 3.3 Validate, extract, analyse

```bash
python scripts/flex/run_study.py --config configs/flex/study.af19.json --validate-only
```
```
Config OK: configs/flex/study.af19.json (af19) -- 1 measures, 2 analyses
```
```bash
python scripts/flex/run_study.py --config configs/flex/study.af19.json --stage extract
```
```
[hippo] wrote 7566 rows (source=hippo-amygdala, dataset=freesurfer)
```
```bash
python scripts/flex/run_study.py --config configs/flex/study.af19.json --stage analyse
```
```
[lmm/hippo] 4074 rows, 102 subjects, 7 ROIs
Wrote .../lmm_hippo/summary.csv (7 rows)
[change_ancova/hippo] baseline=ses-1, follow-up=ses-2, ses-3
Wrote .../change_ancova_hippo/summary.csv (7 rows)
```

Note `change_ancova` printed `follow-up=ses-2, ses-3` — it correctly used
the *multi-followup* branch (`intervention * followup_f` interaction term)
because this study has 3 sessions, versus pk01's 2-session
`intervention`-only branch. Same code, different design, correct behavior
in both — that branch selection is exactly what `test_change_ancova.R`'s
`length(followup_sessions) > 1` check is for.

---

## 4. Extending an existing study's config

These are small, common edits — none require touching `scripts/flex/`.

### 4.1 Add a new measure (e.g. basal ganglia via `aseg`)

Append to `measures[]`:
```json
{
  "name": "basal_ganglia",
  "dataset": "freesurfer",
  "source": "aseg",
  "profile": "volume",
  "roi_set": ["Caudate", "Putamen", "Pallidum", "Accumbens-area"]
}
```
Then add it to whichever `analyses[]` entries should test it, e.g.:
```json
{ "test": "lmm", "measures": ["hippo", "amygdala", "basal_ganglia"] }
```
`aseg` region names come back with hemisphere stripped already (e.g.
`Caudate`, not `Left-Caudate`/`Right-Caudate` — `extract_freesurfer.py`'s
`strip_hemi_prefix()` handles that), so the `roi_set` above matches both
hemispheres' rows in one go.

### 4.2 Add a new factorial axis

```json
"design": {
  "factorial": {
    "responder": { "high": ["arm_a"], "low": ["arm_b", "arm_c"] }
  }
}
```
```json
{ "test": "factorial", "measures": ["hippo"], "factorial_axes": ["responder"] }
```
Any test that calls `flex_derive_factorial()` (currently `factorial` and
`rm_anova`'s `between`/`interaction_factors`) can reference this axis by
name — no new R code needed for a new grouping scheme, only a new JSON
block.

### 4.3 Run a second `rm_anova` variant against the same measure

Multiple analyses can target the same `(test, measure)` pair — give each
one an explicit `id` so their output directories don't collide:

```json
{ "test": "rm_anova", "id": "rm_anova_afex", "measures": ["hippo"], "engine": "afex" },
{ "test": "rm_anova", "id": "rm_anova_by_arm", "measures": ["hippo"], "engine": "afex",
  "between": "group3", "interaction_factors": ["group3", "time_f", "hemisphere"] }
```
Without `id`, both would try to write to `results/rm_anova_hippo/` and the
second run would silently overwrite the first's output — `id` is exactly
the escape hatch for that.

### 4.4 Restrict which sessions are included

```json
"sessions": { "include": ["ses-1", "ses-2"], "baseline": "ses-1" }
```
Useful for e.g. excluding a 3rd wave that isn't complete yet, without
touching the underlying data. Drop `include` entirely to use every session
found.

---

## 5. One test, two measure types: volume vs. thickness

This is the concrete payoff of the `profile` mechanism, not just a claim
about it. The original battery has near-duplicate scripts for hippocampal
subfield *volume* vs. *thickness* — `08_moderator_analysis.R` /
`33_hippo_thickness_moderator_analysis.R`, `23_hem_time_group_cesd.R` /
`36_hippo_thickness_hem_time_group_cesd.R`, `19_age_x_social_duration.R` /
`34_hippo_thickness_age_x_social_duration.R`, and so on — that differ
*only* in the volume-vs-thickness details the `profile` already encodes
(log-transform, eTIV covariate, aggregation FUN). In `flex`, one module
runs both:

```json
{ "test": "moderator", "measures": ["hippo", "hippo_thick"], "moderators": "all" }
```

Real output, side by side (same command, same run, just a different
`measure` picking up the `hippo` (`volume` profile) vs. `hippo_thick`
(`thickness` profile) tidy table):

```bash
Rscript scripts/flex/R/run_analysis.R --config configs/flex/study.pk01.json \
  --only moderator --measures hippo,hippo_thick
```
```
[moderator/hippo] testing 8 candidate moderator(s) across 7 ROI(s)
[moderator/hippo_thick] testing 8 candidate moderator(s) across 4 ROI(s)
```

`results/moderator_hippo/summary.csv` and
`results/moderator_hippo_thick/summary.csv` are produced by the exact same
R code (`test_moderator.R`), diffed against `08_moderator_analysis.R` and
`33_hippo_thickness_moderator_analysis.R` respectively — both matched
their original script's p-values to ~1e-11.

Same story for `rm_anova`'s `lmm_interaction` engine (23 vs. 36):
```json
{ "test": "rm_anova", "id": "rm_anova_colleague_cesd", "measures": ["hippo", "hippo_thick"],
  "engine": "lmm_interaction", "factors": ["hemisphere", "time_f", "group3"],
  "continuous_moderator": "ads_score" }
```

And for `age_interaction` (new module, replacing 19/20 for volume and
34/35 for thickness):
```json
{ "test": "age_interaction", "id": "age_interaction_pairwise", "measures": ["hippo", "hippo_thick"],
  "mode": "pairwise", "factorial_axes": ["social", "duration"] },
{ "test": "age_interaction", "id": "age_interaction_threeway", "measures": ["hippo", "hippo_thick"],
  "mode": "threeway_intervention", "moderators": "all" }
```
`mode: "pairwise"` tests age x each factorial axis (2-way terms only, no
3-way — real output from `age_interaction_pairwise_hippo_thick`):
```
roi,n_obs,n_subjects,age_x_social_p,age_x_duration_p,...
CA1,170,85,0.894,0.119,...
```
`mode: "threeway_intervention"` tests the full intervention x age x
third-variable interaction, once per third variable (sex, plus every
`design.moderators.baseline_cols` column):
```
third_var,roi,n_obs,interaction_p,interaction_p_fdr,significant
sex,CA1,254,0.215,0.430,FALSE
```

**The general pattern**: if you're adding a new *measure* that's
conceptually a volume or a thickness (or genuinely something else — define
a custom profile), you very likely don't need a new test module at all.
Add the measure, list it under whichever `analyses[].measures` you want to
test, done. A new module is only needed for a genuinely new *statistical
design* (a new formula shape), not a new brain region or a new unit of
measurement.

---

## 6. Troubleshooting

| Symptom | Likely cause | Where to look |
|---|---|---|
| `INVALID CONFIG: ... required property` | Missing required schema field | The JSON pointer path in the error message |
| `INVALID CONFIG: /design/contrasts: group level(s) [...] appear in both...` | A group value listed in both `intervention` and `control` | `design.contrasts` |
| `no rows matched roi_set` | ROI names in the config don't match the tidy file's `region` column | Check `tidy/<measure>.tsv`'s `region` column directly, or set `roi_set: "all"` temporarily to see what's actually there |
| `tidy file not found for measure '...'` | Extract stage hasn't run yet, or wrote to a different `output.root` | Run `--stage extract` first; check `output.root` matches between runs |
| Extraction returns 0 rows / stats files look like `/annex/objects/...` | DataLad content not fetched — a pointer file, not real data | Run `--stage ingest` (not `--dry-run`) before extract |
| `unknown measure '...'` in an `analyses[]` entry | Typo, or the measure wasn't defined in `measures[]` | `analyses[].measures` must be a subset of `measures[].name` |
| `design.contrasts.intervention/control not set...required for this test` | Running `lmm`/`change_ancova`/`moderator`/`loso`/`multivariate` on a config with no `contrasts` block | These tests need an intervention/control contrast; a pure observational/cross-sectional study can't use them as-is |
| `SKIP: test '...' has no module yet` | You listed a `test` value in `analyses[]` that matches the schema but isn't implemented (e.g. `pointwise`, `hierarchical`) | See the "not yet ported" list in `FLEX_PIPELINE.md`; use the original `scripts/analysis/NN_*.R` script directly for it |

---

## 7. Where to go next

- **Concepts** (profiles, why each test module exists, parity methodology): [`FLEX_PIPELINE.md`](FLEX_PIPELINE.md)
- **Field-by-field schema reference**: [`configs/flex/README.md`](../configs/flex/README.md)
- **What each test actually computes**: the header comment at the top of each `scripts/flex/R/test_*.R` file — they carry the same rationale as the original scripts they replace
- **A minimal starting point**: `configs/flex/study.TEMPLATE.json`
