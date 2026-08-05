# flex study configs

One JSON file per study, validated against `study.schema.json`. See
`../../docs/FLEX_PIPELINE.md` for the concepts (profiles, test modules,
parity verification) and a porting walkthrough. This file is the field
reference.

- `study.schema.json` — the JSON Schema every config is validated against
  (`python scripts/flex/config.py <path>`, or `run_study.py --validate-only`).
- `study.TEMPLATE.json` — minimal starting point for a new study.
- `study.pk01.json` — fully worked example (the current study).

## Field reference

### `study`
| Field | Required | Description |
|---|---|---|
| `id` | yes | Short identifier, used in default output paths |
| `label` | no | Human-readable description |

### `datasets.<name>`
| Field | Required | Description |
|---|---|---|
| `path` | yes | Directory containing `sub-*_ses-*.long.sub-*` timepoint dirs |
| `datalad` | no | `true` if this is a DataLad sibling to `datalad get` from (triggers the ingest stage) |
| `dir_pattern` | no | Override regex (named groups `base`, `ses`) for the timepoint directory naming |
| `subject_id_map` | no | Path to a short_id/long_id/match_type TSV (see `scripts/build_hipsta_subject_id_map.py`) if this dataset's own directory names don't match `design.participants`' subject IDs |

### `sessions`
| Field | Required | Description |
|---|---|---|
| `include` | no | Restrict to these session values; omit to include every session found |
| `baseline` | no | Session used as the baseline/eTIV-reference timepoint; defaults to the earliest session present |

### `measures[]`
| Field | Required | Description |
|---|---|---|
| `name` | yes | Identifier used in `analyses[].measures` and the output `tidy/<name>.tsv` filename |
| `dataset` | yes | Key into `datasets` |
| `source` | yes | One of `aseg`, `aparc`, `aparc.a2009s`, `hippo-amygdala`, `thalamus`, `brainstem`, `hipsta` |
| `metric` | no | Column to extract for `aseg`/`aparc`/`aparc.a2009s` (default: `Volume_mm3` / `ThickAvg`) |
| `profile` | yes | `volume`, `thickness`, `volume_midline`, or a key in `profiles` |
| `region_filter` | no | `{prefix}` and/or `{include}` to split a multi-region source (e.g. `hippo-amygdala` into separate `hippo`/`amygdala` measures) |
| `roi_set` | no | List of region names, or `"all"` (default) |
| `roi_column` | no | Column identifying the ROI in the tidy file (default: `region`) |
| `region_groups` | no | Path to a `region\tcluster` TSV for `test_multivariate.R`'s per-cluster mode (28's design) |
| `pointwise` | no | For `source: "hipsta"`: also emit a grid-point-level companion tidy file |

### `profiles.<name>` (optional, extends the three built-ins)
`{transform: "log"|"identity", etiv_covariate: true|false, hemisphere: "pooled"|"none"}`

### `design`
| Field | Required | Description |
|---|---|---|
| `participants` | yes | Path to the participants TSV |
| `roles.subject` | yes | Column holding the subject ID |
| `roles.group` | yes | Column holding the group assignment |
| `roles.age`, `roles.sex` | no | Column names, if used as covariates |
| `contrasts.intervention`, `contrasts.control` | no | Group values pooled into each side of the intervention contrast — required by any test that derives one (`lmm`, `change_ancova`, `multivariate`, `loso`, `moderator`); omit entirely for non-intervention designs |
| `factorial.<axis>.<level>` | no | Group values assigned to each level of a named axis (crossed by `test_factorial.R`, or used as `test_rm_anova.R`'s `between`) |
| `covariates[]` | no | `{col, type: "numeric"|"factor", transform: "none"|"z"|"rank_z", as}` — applied uniformly across every test |
| `moderators.file` | no | Subject-level moderators TSV |
| `moderators.baseline_cols` | no | Columns in it to test as candidate moderators |

### `analyses[]`
| Field | Required | Description |
|---|---|---|
| `test` | yes | One of `lmm`, `change_ancova`, `multivariate`, `moderator`, `factorial`, `loso`, `rm_anova` (others are schema-valid but not yet implemented — see docs/FLEX_PIPELINE.md) |
| `measures` | yes | Measure names to run this test against |
| `id` | no | Output-folder discriminator (`<id>_<measure>`) when multiple analyses share one test type against one measure (e.g. two differently-configured `rm_anova` runs) |
| *(test-specific)* | no | e.g. `moderators` (moderator), `factorial_axes` (factorial), `engine`/`between`/`interaction_factors`/`factors`/`continuous_moderator` (rm_anova), `keep_hemisphere_separate` (multivariate) — see each module's header comment in `scripts/flex/R/test_*.R` |

### `output`
| Field | Required | Description |
|---|---|---|
| `root` | yes | Output directory root (`tidy/`, `results/`, `logs/` go under it) |
| `alpha` | no | Significance threshold (default 0.05) |
| `fdr_within` | no | `"measure"` (default) or `"analysis"` |
