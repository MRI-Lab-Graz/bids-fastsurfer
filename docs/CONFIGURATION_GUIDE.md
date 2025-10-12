# Configuration Files Summary

This document provides an overview of all configuration templates and documentation created for the FastSurfer longitudinal analysis pipeline.

## 📁 Files Created

### Configuration Templates

| File | Lines | Purpose | Use When |
|------|-------|---------|----------|
| **fslmer_univariate.TEMPLATE.json** | Comprehensive | All available parameters with detailed comments | Learning all options, creating custom configs |
| **fslmer_univariate.minimal.json** | Minimal | Quick start template | Simple single-ROI analysis |
| **meta_pipeline.TEMPLATE.json** | Comprehensive | All meta-pipeline parameters | Full pipeline with multiple models |
| **meta_pipeline.minimal.json** | Minimal | Quick start meta-pipeline | Simple end-to-end workflow |

### Documentation

| File | Purpose | Audience |
|------|---------|----------|
| **configs/README.md** | Complete configuration guide | All users |
| **docs/CONFIG_QUICK_REFERENCE.md** | Quick reference card | Quick lookups |

## 🎯 Quick Selection Guide

### "I want to analyze a single brain region"
→ Use **fslmer_univariate.minimal.json**

```bash
cp configs/fslmer_univariate.minimal.json configs/my_analysis.json
# Edit my_analysis.json with your paths and settings
bash scripts/run_fslmer_univariate.sh --config configs/my_analysis.json
```

### "I want to run a full end-to-end pipeline"
→ Use **meta_pipeline.minimal.json**

```bash
cp configs/meta_pipeline.minimal.json configs/my_study.json
# Edit my_study.json with your paths
bash scripts/run_meta_pipeline.sh --config configs/my_study.json
```

### "I want to see all available options"
→ Check **TEMPLATE.json** files

```bash
# For single model analysis
less configs/fslmer_univariate.TEMPLATE.json

# For meta-pipeline
less configs/meta_pipeline.TEMPLATE.json
```

### "I need quick syntax help"
→ Use **CONFIG_QUICK_REFERENCE.md**

```bash
less docs/CONFIG_QUICK_REFERENCE.md
```

## 📊 Configuration Comparison

### fslmer_univariate configs

| Feature | TEMPLATE | minimal | example |
|---------|----------|---------|---------|
| All parameters | ✅ | ❌ | ❌ |
| Detailed comments | ✅ | ❌ | ⚠️ |
| Ready to use | ⚠️ | ✅ | ✅ |
| Learning tool | ✅ | ❌ | ⚠️ |
| Production use | ❌ | ✅ | ✅ |

### meta_pipeline configs

| Feature | TEMPLATE | minimal | example | inline |
|---------|----------|---------|---------|--------|
| All parameters | ✅ | ❌ | ⚠️ | ⚠️ |
| Detailed comments | ✅ | ❌ | ⚠️ | ⚠️ |
| Model keys | ✅ | ✅ | ✅ | ❌ |
| Inline models | ✅ | ❌ | ❌ | ✅ |
| Effect lists | ✅ | ❌ | ✅ | ✅ |
| Ready to use | ⚠️ | ✅ | ✅ | ✅ |

## 🔑 Key Parameters by Category

### Input Files (Required)
```json
{
  "qdec": "path/to/qdec.table.dat",
  "aseg": "path/to/aseg.long.table"
}
```

### ROI Selection (Choose ONE)
```json
{
  "roi": "Left.Hippocampus",              // Single region
  "region_pattern": "Hippocampus|Amygdala", // Multiple
  "all_regions": true                      // All regions
}
```

### Statistical Model
```json
{
  "formula": "~ factor(tp) * group + age",
  "random_effects": "RIRS",
  "add_baseline": true
}
```

### Output
```json
{
  "outdir": "results/my_analysis",
  "save_merged": false
}
```

### Pipeline-Specific (meta_pipeline)
```json
{
  "subjects_dir": "/path/to/fastsurfer",
  "participants": "configs/participants.tsv",
  "out_root": "results/study",
  "models": ["group5_factor_time"],
  "link_long": true
}
```

## 📖 Documentation Structure

```
fast-surfer-workshop/
├── configs/
│   ├── README.md                        # Complete guide
│   ├── fslmer_univariate.TEMPLATE.json  # All options reference
│   ├── fslmer_univariate.minimal.json   # Quick start
│   ├── fslmer_univariate.example.json   # Working example
│   ├── meta_pipeline.TEMPLATE.json      # Pipeline reference
│   ├── meta_pipeline.minimal.json       # Pipeline quick start
│   └── meta_pipeline.example.json       # Pipeline example
│
├── docs/
│   ├── CONFIG_QUICK_REFERENCE.md        # Quick reference card
│   ├── common_functions_reference.md    # Function library docs
│   └── fslmer_univariate_examples.md    # Analysis examples
│
└── REFACTORING_SUMMARY.md               # Code refactoring overview
```

## 🚀 Common Workflows

### 1. Hippocampus Time Effect

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "roi": "Left.Hippocampus",
  "formula": "~ factor(tp)",
  "random_effects": "RIRS",
  "outdir": "results/hippo_time"
}
```

### 2. Limbic Group × Time

```json
{
  "qdec": "results/qdec.table.dat",
  "aseg": "results/aseg.long.table",
  "region_pattern": "Hippocampus|Amygdala",
  "formula": "~ factor(tp) * group",
  "random_effects": "RIRS",
  "outdir": "results/limbic_interaction"
}
```

### 3. Full Pipeline Multiple Models

```json
{
  "subjects_dir": "/data/fastsurfer",
  "participants": "configs/participants.tsv",
  "out_root": "results/intervention",
  "models": ["group5_factor_time", "group5_decomposed"],
  "link_long": true,
  "qc": true,
  "effects": {
    "group5_factor_time": ["tp2", "tp3", "sex", "age"]
  }
}
```

## 🎓 Learning Path

### Beginner
1. Read `configs/README.md` - Overview and basics
2. Use `fslmer_univariate.minimal.json` - First analysis
3. Check `CONFIG_QUICK_REFERENCE.md` - Quick syntax help

### Intermediate
1. Explore `fslmer_univariate.TEMPLATE.json` - All options
2. Try `meta_pipeline.minimal.json` - Full pipeline
3. Read `docs/fslmer_univariate_examples.md` - Detailed examples

### Advanced
1. Study `meta_pipeline.TEMPLATE.json` - Advanced features
2. Review `common_functions_reference.md` - Extend scripts
3. Check `REFACTORING_SUMMARY.md` - Architecture understanding

## 🔧 Validation and Testing

### 1. JSON Syntax Check
```bash
python -m json.tool configs/your_config.json
```

### 2. Dry Run
```bash
bash scripts/run_meta_pipeline.sh --config config.json --dry-run
```

### 3. Print Available Data
```bash
Rscript scripts/fslmer_univariate.R \
  --qdec data.dat --aseg table.txt --print-cols
```

### 4. Help Text
```bash
bash scripts/run_fslmer_univariate.sh --help
bash scripts/run_meta_pipeline.sh --help
Rscript scripts/fslmer_univariate.R --help
```

## 💡 Pro Tips

### Configuration
- ✅ Start with minimal templates
- ✅ Add complexity incrementally
- ✅ Test with single ROI before whole-brain
- ✅ Use version control for configs
- ✅ Document your choices in comments

### Formulas
- ✅ Use `factor(tp)` for discrete timepoints
- ✅ Use `RIRS` random effects for longitudinal data
- ✅ Include relevant covariates (age, sex, eTIV)
- ✅ Add baseline correction when appropriate

### Workflow
- ✅ Validate JSON before running
- ✅ Use dry-run to preview
- ✅ Check outputs incrementally
- ✅ Keep configs with results for reproducibility

## 🆘 Getting Help

### By Topic

| Question | See |
|----------|-----|
| "What parameters are available?" | TEMPLATE.json files |
| "How do I write formulas?" | CONFIG_QUICK_REFERENCE.md |
| "What's the full workflow?" | configs/README.md |
| "How do I customize?" | TEMPLATE.json + README.md |
| "Common errors?" | configs/README.md (Troubleshooting) |

### By File Type

| File Type | Purpose | Documentation |
|-----------|---------|---------------|
| `.TEMPLATE.json` | Reference | Inline comments |
| `.minimal.json` | Quick start | See README.md |
| `.example.json` | Working example | See README.md |

## 📝 File Purposes Summary

### Templates (Learn from these)
- **Purpose**: Show all available options with explanations
- **Use**: Reference when creating custom configs
- **Edit**: No, copy and modify

### Minimal (Start with these)
- **Purpose**: Quick start with essential parameters only
- **Use**: First analyses, simple workflows
- **Edit**: Yes, customize for your data

### Examples (Copy these)
- **Purpose**: Real working configurations
- **Use**: Production analyses, as-is or customized
- **Edit**: Yes, adapt to your needs

## 🎯 Next Steps

After setting up your configuration:

1. **Validate syntax**
   ```bash
   python -m json.tool your_config.json
   ```

2. **Test with dry-run**
   ```bash
   bash scripts/run_meta_pipeline.sh --config your_config.json --dry-run
   ```

3. **Run analysis**
   ```bash
   bash scripts/run_meta_pipeline.sh --config your_config.json
   ```

4. **Check outputs**
   ```bash
   ls -lh results/your_output_dir/
   ```

## 📚 Related Documentation

- `REFACTORING_SUMMARY.md` - Code architecture
- `docs/common_functions_reference.md` - Function library
- `docs/fslmer_univariate_examples.md` - Analysis examples
- Script help: `--help` flags on all scripts

---

**Questions?** Check the configs/README.md or run scripts with `--help`
