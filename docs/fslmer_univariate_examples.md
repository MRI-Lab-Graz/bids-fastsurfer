# fslmer_univariate.R - Usage Examples and Documentation

## Overview

The `fslmer_univariate.R` script performs flexible univariate linear mixed-effects (LME) analysis for subcortical/cortical volume data from FreeSurfer/FastSurfer longitudinal processing. It supports single ROI or multi-region analysis with various statistical models.

## Quick Start

### 1. Basic Single ROI Analysis

```bash
# Analyze Left Hippocampus volume changes over time
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --roi Left.Hippocampus \
  --formula '~ tp' \
  --engine glm \
  --outdir results/hippo_basic
```

### 2. Group × Time Interaction Analysis

```bash
# Test for different trajectories between groups
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --roi Left.Hippocampus \
  --formula '~ tp*group_5' \
  --engine glm \
  --outdir results/hippo_group_interaction
```

### 3. Multi-Region Analysis with Pattern Matching

```bash
# Analyze all hippocampus and amygdala regions
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --region-pattern 'Hippocampus|Amygdala' \
  --formula '~ tp*group_5' \
  --engine glm \
  --outdir results/limbic_analysis
```

### 4. Comprehensive Analysis of All Brain Regions

```bash
# Analyze all subcortical regions (excludes summary measures)
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --all-regions \
  --formula '~ tp' \
  --engine glm \
  --outdir results/whole_brain_longitudinal
```

### 5. Using Configuration Files

```bash
# Use pre-configured analysis
bash scripts/run_fslmer_univariate.sh \
  --config configs/fslmer_test_single_roi.json \
  --engine glm
```

## Configuration File Examples

### Single ROI Configuration
File: `configs/fslmer_test_single_roi.json`
```json
{
  "qdec": "results/running_intervention_inline/qdec.table.dat",
  "aseg": "results/running_intervention_inline/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ tp*group_5",
  "random-effects": "RIRS",
  "outdir": "test_fslmer_single_roi",
  "time_col": "tp",
  "id_col": "Measure.volume",
  "save_merged": true
}
```

### Multi-Region Configuration  
File: `configs/fslmer_test_limbic.json`
```json
{
  "qdec": "results/running_intervention_inline/qdec.table.dat",
  "aseg": "results/running_intervention_inline/aseg.long.table",
  "region_pattern": "Hippocampus|Amygdala",
  "formula": "~ tp*group_5",
  "random-effects": "RIRS",
  "outdir": "test_fslmer_limbic",
  "time_col": "tp",
  "id_col": "Measure.volume",
  "save_merged": true
}
```

## Advanced Examples

### 1. Controlling for Covariates

```bash
# Include age, sex, and brain volume as covariates
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --roi Left.Hippocampus \
  --formula '~ tp*group_5 + age + sex + EstimatedTotalIntraCranialVol' \
  --engine glm \
  --outdir results/hippo_with_covariates
```

### 2. Using Baseline Values as Covariates

```bash
# Add baseline volume as covariate
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --roi Left.Hippocampus \
  --add-baseline \
  --formula '~ baseline_value + tp*group_5' \
  --engine glm \
  --outdir results/hippo_baseline_adjusted
```

### 3. GAM Analysis with Smooth Time Effects

```bash
# Use GAM for non-linear time effects
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --roi Left.Hippocampus \
  --engine gam \
  --formula '~ s(tp) + group_5' \
  --outdir results/hippo_gam_smooth
```

### 4. Quadratic Time Effects

```bash
# Model quadratic change over time
bash scripts/run_fslmer_univariate.sh \
  --qdec results/prep_long/qdec.table.dat \
  --aseg results/prep_long/aseg.long.table \
  --roi Left.Hippocampus \
  --formula '~ tp + I(tp^2) + group_5' \
  --engine glm \
  --outdir results/hippo_quadratic
```

## Understanding Output Files

### Single ROI Analysis Output
- `fit.rds`: R object containing the fitted model
- `used_config.json`: Configuration used for the analysis

### Multi-Region Analysis Output
- `lme_coefficients.csv`: Summary table of all coefficients across ROIs
- `individual_rois/`: Directory with individual model files per ROI
- `used_config.json`: Configuration used for the analysis

## Common Use Cases

### 1. Intervention Study Analysis
```bash
# Test intervention effect over time
--formula '~ tp*interv_group + baseline_covariates'
```

### 2. Disease Progression Study
```bash
# Model disease progression with age effects
--formula '~ tp*diagnosis + age + sex + education'
```

### 3. Development Study
```bash
# Non-linear developmental trajectories
--engine gam --formula '~ s(age) + sex'
```

## Data Requirements

### qdec.table.dat Format
- Tab-separated values (TSV)
- Required columns: `fsid`, `fsid_base`, `tp` (time point)
- Additional covariates as needed (e.g., `group_5`, `age`, `sex`)

### aseg.long.table Format  
- Whitespace-separated values
- Required: `Measure.volume` column with FreeSurfer IDs
- Format: `subject_session.long.subject_base`
- ROI columns with volume measurements

## Troubleshooting

### Common Issues and Solutions

1. **Missing fslmer package**: Use `--engine glm` or `--engine gam`
2. **No ROIs found**: Check ROI names with `--print-cols`
3. **Formula errors**: Test formula syntax in R: `as.formula("~ tp*group")`
4. **File not found**: Verify paths and run `prep_long.py` first

### Getting Help

```bash
# Show comprehensive help
bash scripts/run_fslmer_univariate.sh --help

# Show detailed R argument documentation  
Rscript scripts/fslmer_univariate.R --help

# Preview data columns
Rscript scripts/fslmer_univariate.R --print-cols \
  --qdec your_qdec.dat --aseg your_aseg.table
```

## Environment Setup

The script automatically uses the project's micromamba environment if available, otherwise falls back to system R. To set up the recommended environment:

```bash
# Install required packages and environment
bash scripts/install.sh
```