# Configuration Quick Reference Card

## Essential Parameters

### run_fslmer_univariate.sh

```json
{
  "qdec": "path/to/qdec.table.dat",     // REQUIRED: Subject metadata
  "aseg": "path/to/aseg.long.table",    // REQUIRED: Volume data
  
  "roi": "Left.Hippocampus",            // Single region OR
  "region_pattern": "Hippocampus|Amygdala",  // Multiple regions OR
  "all_regions": true,                  // All regions
  
  "formula": "~ factor(tp) * group + age + sex",  // R formula
  "random_effects": "RIRS",             // RI, RS, or RIRS
  "add_baseline": true,                 // Baseline correction
  
  "outdir": "results/output",           // Output directory
  "alpha": 0.05                         // Significance level
}
```

### run_meta_pipeline.sh

```json
{
  "subjects_dir": "/path/to/fastsurfer",     // REQUIRED
  "participants": "configs/participants.tsv", // REQUIRED
  "out_root": "results/study",               // Output root
  
  "models": ["group5_factor_time"],     // Model keys/paths
  "link_long": true,                    // Create .long symlinks
  
  "effects": {                          // Effects to summarize
    "group5_factor_time": ["tp2", "tp3", "sex"]
  }
}
```

## Formula Cheat Sheet

| Pattern | Meaning | Example |
|---------|---------|---------|
| `~ x` | Main effect | `~ tp` |
| `~ x + y` | Additive effects | `~ tp + group` |
| `~ x * y` | Interaction + main | `~ tp * group` |
| `~ x:y` | Interaction only | `~ tp:group` |
| `~ factor(x)` | Categorical | `~ factor(tp)` |
| `~ I(x^2)` | Polynomial | `~ tp + I(tp^2)` |
| `~ s(x)` | Smooth (GAM) | `~ s(tp)` |

## Random Effects

| Code | Meaning | Use When |
|------|---------|----------|
| `"RI"` | Random intercept | Different baselines, same slope |
| `"RS"` | Random slope | Same baseline, different slopes |
| `"RIRS"` | Both | Different baselines AND slopes ⭐ |

⭐ **Recommended for longitudinal data**

## ROI Selection Patterns

```json
// Single ROI
"roi": "Left.Hippocampus"

// Limbic regions
"region_pattern": "Hippocampus|Amygdala"

// All ventricles
"region_pattern": "Ventricle"

// Left hemisphere only
"region_pattern": "^Left\\."

// Everything except summary measures
"all_regions": true

// Everything including summaries
"all_regions": true,
"include_summary": true
```

## Effect Name Shortcuts

| Shortcut | Expands To | Notes |
|----------|------------|-------|
| `tp2` | `factor(tp)2` | Time point 2 |
| `tp3` | `factor(tp)3` | Time point 3 |
| `sex` | `sexM` | Male effect |
| `age` | `age` | Age covariate |
| `tp2:smallgroup_2w` | `factor(tp)2:group_5smallgroup_2w` | Interaction |

## Common Patterns

### Time Only
```json
{
  "formula": "~ factor(tp)",
  "random_effects": "RI"
}
```

### Time × Group
```json
{
  "formula": "~ factor(tp) * group",
  "random_effects": "RIRS"
}
```

### Full Model
```json
{
  "formula": "~ factor(tp) * group_5 + age + sex + EstimatedTotalIntraCranialVol",
  "random_effects": "RIRS",
  "add_baseline": true
}
```

### GAM Smooth
```json
{
  "engine": "gam",
  "formula": "~ s(tp) + group",
  "gam_k": 4
}
```

## File Locations

```
configs/
├── fslmer_univariate.TEMPLATE.json    # Complete reference
├── fslmer_univariate.minimal.json     # Quick start
├── meta_pipeline.TEMPLATE.json        # Pipeline reference
└── meta_pipeline.minimal.json         # Pipeline quick start
```

## Command Line

```bash
# Single model
bash scripts/run_fslmer_univariate.sh --config configs/my_config.json

# Meta pipeline
bash scripts/run_meta_pipeline.sh --config configs/meta_config.json

# Dry run (preview only)
bash scripts/run_meta_pipeline.sh --config config.json --dry-run

# Override settings
bash scripts/run_meta_pipeline.sh \
  --config config.json \
  --out-root results/new_run \
  --models group5_factor_time
```

## Quick Checks

```bash
# Validate JSON
python -m json.tool config.json

# See available columns
Rscript scripts/fslmer_univariate.R \
  --qdec qdec.dat --aseg aseg.table --print-cols

# Help
bash scripts/run_fslmer_univariate.sh --help
bash scripts/run_meta_pipeline.sh --help
```

## Pro Tips

✅ **Do:**
- Use `factor(tp)` for discrete timepoints
- Use `RIRS` random effects for longitudinal data
- Add baseline correction for individual differences
- Test with single ROI before whole-brain
- Use dry-run to preview execution

❌ **Don't:**
- Mix ROI selection methods (choose one)
- Forget factor() for categorical time
- Use underscores in ROI names (use dots)
- Skip validation of JSON syntax

## Need More?

- Full templates: `configs/*.TEMPLATE.json`
- Complete guide: `configs/README.md`
- Examples: `docs/fslmer_univariate_examples.md`
- Functions: `docs/common_functions_reference.md`
