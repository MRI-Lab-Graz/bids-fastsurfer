# Configuration Files Guide

This directory contains configuration files for the FastSurfer longitudinal analysis pipeline.

## Quick Start

### For Single-Model Analysis (run_fslmer_univariate.sh)

**Minimal example:**
```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ factor(tp) * group",
  "outdir": "results/hippo_analysis"
}
```

**Full template with all options:**
See `fslmer_univariate.TEMPLATE.json`

### For Meta-Pipeline (run_meta_pipeline.sh)

**Minimal example:**
```json
{
  "subjects_dir": "/path/to/derivatives/fastsurfer",
  "participants": "configs/participants.tsv",
  "out_root": "results/my_study",
  "models": ["group5_factor_time"],
  "link_long": true
}
```

**Full template with all options:**
See `meta_pipeline.TEMPLATE.json`

## Available Templates

| File | Purpose | When to Use |
|------|---------|------------|
| `fslmer_univariate.TEMPLATE.json` | Complete reference with ALL options | Learning all available parameters |
| `fslmer_univariate.minimal.json` | Minimal working example | Quick start for simple analyses |
| `fslmer_univariate.example.json` | Basic example | Standard single-ROI analysis |
| `meta_pipeline.TEMPLATE.json` | Complete meta-pipeline reference | Learning full pipeline workflow |
| `meta_pipeline.minimal.json` | Minimal meta-pipeline example | Quick start for full pipeline |
| `meta_pipeline.example.json` | Standard pipeline example | Typical intervention study |
| `meta_pipeline.inline.pk01.json` | Inline model definitions | Defining models directly in config |

## Existing Study Configurations

| File | Description |
|------|-------------|
| `fslmer_group5_factor_time.json` | Factorial time × group design |
| `fslmer_group5_decomposed.json` | Decomposed intervention effects |
| `fslmer_group_beh.json` | Behavioral intervention model |
| `fslmer_group5.json` | Five-group comparison |
| `fslmer_group_time.json` | Time × group interaction |
| `fslmer_test_limbic.json` | Limbic region test case |
| `fslmer_test_single_roi.json` | Single ROI test case |

## Configuration Formats

### Format 1: Single Model Config (fslmer_univariate)

Used with `run_fslmer_univariate.sh`:

```json
{
  "qdec": "path/to/qdec.table.dat",
  "aseg": "path/to/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ factor(tp) * group + age + sex",
  "random_effects": "RIRS",
  "add_baseline": true,
  "outdir": "results/output"
}
```

### Format 2: Meta-Pipeline with Model Keys

Used with `run_meta_pipeline.sh`:

```json
{
  "subjects_dir": "/path/to/fastsurfer",
  "participants": "configs/participants.tsv",
  "out_root": "results/study",
  "models": ["group5_factor_time", "group5_decomposed"],
  "effects": {
    "group5_factor_time": ["tp2", "tp3", "sex", "age"]
  }
}
```

### Format 3: Inline Model Definitions

Define models directly in meta config:

```json
{
  "subjects_dir": "/path/to/fastsurfer",
  "participants": "configs/participants.tsv",
  "out_root": "results/study",
  "models": [
    {
      "id": "simple_model",
      "formula": "~ factor(tp)",
      "zcols": "RI"
    },
    {
      "id": "full_model",
      "formula": "~ factor(tp) * group + age",
      "zcols": "RIRS"
    }
  ]
}
```

## Key Parameters Explained

### ROI Selection (choose ONE)

- **`roi`**: Single region name
  - Example: `"Left.Hippocampus"`
  - Use for focused analysis

- **`region_pattern`**: Regex pattern
  - Example: `"Hippocampus|Amygdala"` (limbic regions)
  - Example: `"^Left\\."` (all left hemisphere)
  - Use for multiple related regions

- **`all_regions`**: Boolean
  - Set to `true` for whole-brain analysis
  - Excludes summary measures by default

### Formula Syntax

The `formula` parameter uses R formula syntax:

```r
~ tp                           # Simple time effect (continuous)
~ factor(tp)                   # Time as categorical factor
~ tp * group                   # Time × group interaction
~ factor(tp) * group_5 + age   # Factorial with covariates
~ tp + I(tp^2)                 # Quadratic time trend
~ s(tp)                        # GAM smooth (requires engine: "gam")
```

**Important:**
- Use `factor(tp)` to treat time as categorical (recommended for discrete timepoints)
- Use `tp` alone for continuous time trends
- Interactions: `*` expands to main effects + interaction (a*b = a + b + a:b)

### Random Effects

- **`RI`**: Random intercept only
  - Each subject has different baseline
  - Assumes parallel trajectories
  
- **`RS`**: Random slope only
  - Each subject has different rate of change
  - Assumes same baseline
  
- **`RIRS`**: Random intercept and slope (recommended)
  - Each subject has different baseline AND rate
  - Most flexible, best for longitudinal data
  
- **`"1,2"`**: Legacy index format
  - Deprecated, use named format above

### Effect Names

For summarization, you can use human-friendly names that are auto-converted:

| Friendly Name | Converts To | Meaning |
|--------------|-------------|---------|
| `tp2` | `factor(tp)2` | Time point 2 effect |
| `tp3` | `factor(tp)3` | Time point 3 effect |
| `sex` | `sexM` | Male sex effect |
| `age` | `age` | Age covariate |
| `tp2:smallgroup_2w` | `factor(tp)2:group_5smallgroup_2w` | Interaction |

## Common Use Cases

### 1. Simple Time Effect

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ factor(tp)",
  "random_effects": "RIRS",
  "outdir": "results/time_effect"
}
```

### 2. Group × Time Interaction

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ factor(tp) * group",
  "random_effects": "RIRS",
  "outdir": "results/group_by_time"
}
```

### 3. Full Model with Covariates

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ factor(tp) * group_5 + age + sex + EstimatedTotalIntraCranialVol",
  "random_effects": "RIRS",
  "add_baseline": true,
  "outdir": "results/full_model"
}
```

### 4. Multiple Regions (Limbic System)

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "region_pattern": "Hippocampus|Amygdala",
  "formula": "~ factor(tp) * group",
  "random_effects": "RIRS",
  "outdir": "results/limbic_analysis"
}
```

### 5. Whole-Brain Analysis

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "all_regions": true,
  "formula": "~ factor(tp)",
  "random_effects": "RI",
  "outdir": "results/whole_brain"
}
```

## Testing Your Configuration

### 1. Validate JSON Syntax

```bash
# Check if JSON is valid
python -m json.tool configs/your_config.json
```

### 2. Dry Run

```bash
# See what would be executed without running
bash scripts/run_meta_pipeline.sh --config configs/your_config.json --dry-run
```

### 3. Print Available Columns

```bash
# See what data columns are available
Rscript scripts/fslmer_univariate.R \
  --qdec results/qdec.table.dat \
  --aseg results/aseg.long.table \
  --print-cols
```

## Tips and Best Practices

### 1. Start Small
- Test with a single ROI before running whole-brain
- Use minimal configs first, add complexity incrementally

### 2. Use Absolute Paths
- Recommended for `subjects_dir`, `bids_root`
- Relative paths work from project root

### 3. Comment Your Configs
- JSON doesn't support comments natively
- Use `"_comment"` keys (these are ignored by the scripts)

### 4. Version Control
- Keep configs in git
- Name descriptively: `intervention_study_v2.json`

### 5. Baseline Correction
- Use `add_baseline: true` for individual difference correction
- Add `baseline_value` to formula when enabled

### 6. Model Comparison
- Run multiple models with different formulas
- Compare AIC/BIC values in results

## Troubleshooting

### Error: "Config file not found"
**Solution:** Check path is correct, use absolute or project-relative path

### Error: "ROI not found in aseg table"
**Solution:** 
- Check exact column name with `--print-cols`
- ROI names are case-sensitive
- Use dots not underscores: `Left.Hippocampus` not `Left_Hippocampus`

### Error: "Formula parse error"
**Solution:**
- Ensure proper R syntax
- Use `factor(tp)` for categorical time
- Check variable names match QDEC columns

### Error: "Random effects specification invalid"
**Solution:**
- Use: `"RI"`, `"RS"`, or `"RIRS"`
- Or legacy format: `"1,2"` or `[1,2]`

### Effect not found in coefficients
**Solution:**
- Check time levels: you might have tp2, tp3 but not tp4
- Verify factor() wrapping in formula
- Use human-friendly names that auto-convert

## Additional Resources

- **Full documentation:** See `docs/common_functions_reference.md`
- **Refactoring overview:** See `REFACTORING_SUMMARY.md`
- **Script help:** Run with `--help` flag
  ```bash
  bash scripts/run_fslmer_univariate.sh --help
  bash scripts/run_meta_pipeline.sh --help
  ```
- **R script details:**
  ```bash
  Rscript scripts/fslmer_univariate.R --help
  ```

## Examples Directory

See `docs/fslmer_univariate_examples.md` for detailed analysis examples and interpretation guides.

## Questions?

For configuration questions:
1. Check the TEMPLATE files for parameter descriptions
2. Review existing configs for working examples
3. Run with `--dry-run` to preview execution
4. Use `--print-cols` to explore available data
