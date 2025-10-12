# Configuration Files - Complete Index

## 📋 All Configuration Files

### ✨ New Template Files (Just Created)

| File | Type | Lines | Purpose |
|------|------|-------|---------|
| `fslmer_univariate.TEMPLATE.json` | Template | Full | Complete reference with ALL available parameters |
| `fslmer_univariate.minimal.json` | Starter | 14 | Quick start - essential parameters only |
| `meta_pipeline.TEMPLATE.json` | Template | Full | Complete meta-pipeline reference |
| `meta_pipeline.minimal.json` | Starter | 9 | Quick start - minimal pipeline config |

### 📚 New Documentation Files

| File | Purpose |
|------|---------|
| `configs/README.md` | Complete configuration guide with examples |
| `docs/CONFIG_QUICK_REFERENCE.md` | Quick reference card for syntax |
| `docs/CONFIGURATION_GUIDE.md` | Overview and learning path |

### 📁 Existing Configuration Files

| File | Model Type | Description |
|------|------------|-------------|
| `fslmer_univariate.example.json` | Single model | Basic univariate example |
| `fslmer_group5_factor_time.json` | Factorial | Time × 5-group factorial design |
| `fslmer_group5_decomposed.json` | Decomposed | Decomposed intervention effects |
| `fslmer_group5.json` | Group comparison | Five-group comparison |
| `fslmer_group_time.json` | Interaction | Time × group interaction |
| `fslmer_group_beh.json` | Behavioral | Behavioral intervention model |
| `fslmer_test_limbic.json` | Test case | Limbic region test |
| `fslmer_test_single_roi.json` | Test case | Single ROI test |
| `meta_pipeline.example.json` | Pipeline | Standard meta-pipeline |
| `meta_pipeline.inline.pk01.json` | Pipeline | Inline model definitions |
| `meta_pipeline.pk01.json` | Pipeline | PK01 study configuration |
| `analyse_qdec.example.json` | QDEC | QDEC analysis example |

## 🎯 Quick Navigation

### "I'm new to this"
1. Start → `configs/README.md`
2. Then → `fslmer_univariate.minimal.json`
3. Reference → `docs/CONFIG_QUICK_REFERENCE.md`

### "I need all options"
1. Single model → `fslmer_univariate.TEMPLATE.json`
2. Pipeline → `meta_pipeline.TEMPLATE.json`
3. Functions → `docs/common_functions_reference.md`

### "I want examples"
1. Simple → `fslmer_univariate.minimal.json`
2. Standard → `fslmer_univariate.example.json`
3. Complex → `fslmer_group5_factor_time.json`
4. Pipeline → `meta_pipeline.example.json`

### "I need help"
1. FAQ → `configs/README.md` (bottom)
2. Quick ref → `docs/CONFIG_QUICK_REFERENCE.md`
3. Full guide → `docs/CONFIGURATION_GUIDE.md`

## 📖 Documentation Hierarchy

```
1. Quick Start
   └─ fslmer_univariate.minimal.json
   └─ meta_pipeline.minimal.json

2. Learning
   └─ configs/README.md (start here)
   └─ docs/CONFIG_QUICK_REFERENCE.md
   └─ docs/CONFIGURATION_GUIDE.md

3. Reference
   └─ fslmer_univariate.TEMPLATE.json
   └─ meta_pipeline.TEMPLATE.json
   └─ docs/common_functions_reference.md

4. Examples
   └─ fslmer_univariate.example.json
   └─ fslmer_group5_*.json
   └─ meta_pipeline.example.json

5. Architecture
   └─ REFACTORING_SUMMARY.md
   └─ scripts/common_functions.sh
```

## 🔍 Find By Task

### Task: Single ROI Analysis
**Config:** `fslmer_univariate.minimal.json`
**Docs:** `configs/README.md` → "Common Use Cases" → #1

### Task: Multiple Regions
**Config:** `fslmer_univariate.TEMPLATE.json` (see region_pattern)
**Example:** `fslmer_test_limbic.json`
**Docs:** `configs/README.md` → "ROI Selection Patterns"

### Task: Whole Brain
**Config:** `fslmer_univariate.TEMPLATE.json` (see all_regions)
**Docs:** `configs/README.md` → "Common Use Cases" → #5

### Task: Full Pipeline
**Config:** `meta_pipeline.minimal.json`
**Example:** `meta_pipeline.example.json`
**Docs:** `configs/README.md` → "Format 2"

### Task: Inline Models
**Example:** `meta_pipeline.inline.pk01.json`
**Docs:** `configs/README.md` → "Format 3"

### Task: Custom Formulas
**Reference:** `docs/CONFIG_QUICK_REFERENCE.md` → "Formula Cheat Sheet"
**Docs:** `configs/README.md` → "Formula Syntax"

### Task: Effect Summarization
**Template:** `meta_pipeline.TEMPLATE.json` (see effects section)
**Example:** `meta_pipeline.example.json`
**Docs:** `configs/README.md` → "Effect Names"

## 📊 Feature Matrix

| Feature | minimal | example | TEMPLATE | Notes |
|---------|---------|---------|----------|-------|
| All parameters | ❌ | ⚠️ | ✅ | TEMPLATE shows everything |
| Comments/docs | ⚠️ | ⚠️ | ✅ | TEMPLATE heavily documented |
| Ready to use | ✅ | ✅ | ❌ | Edit paths in minimal/example |
| ROI selection | ✅ | ✅ | ✅ | All support single/multi/all |
| Random effects | ✅ | ✅ | ✅ | New: named format (RIRS) |
| Baseline correction | ✅ | ❌ | ✅ | add_baseline parameter |
| Advanced options | ❌ | ❌ | ✅ | GAM, contrasts, etc. |
| Effect lists | ❌ | ❌ | ✅ | For meta-pipeline |

## 🎓 By Experience Level

### Beginner
**Start here:**
- `configs/README.md` - Read first
- `fslmer_univariate.minimal.json` - First config
- `docs/CONFIG_QUICK_REFERENCE.md` - Quick help

**Commands:**
```bash
cp configs/fslmer_univariate.minimal.json configs/my_first.json
# Edit my_first.json
bash scripts/run_fslmer_univariate.sh --config configs/my_first.json
```

### Intermediate
**Explore:**
- `fslmer_univariate.TEMPLATE.json` - All options
- `meta_pipeline.minimal.json` - Pipeline start
- `fslmer_group5_factor_time.json` - Complex example

**Try:**
```bash
# Multi-region analysis
cp configs/fslmer_test_limbic.json configs/my_limbic.json

# Full pipeline
cp configs/meta_pipeline.minimal.json configs/my_pipeline.json
```

### Advanced
**Reference:**
- `meta_pipeline.TEMPLATE.json` - All pipeline options
- `docs/common_functions_reference.md` - Extend scripts
- `REFACTORING_SUMMARY.md` - Architecture

**Customize:**
```bash
# Inline model definitions
cp configs/meta_pipeline.inline.pk01.json configs/custom_models.json

# Create new scripts using common_functions.sh
```

## 🔗 Cross-References

### configs/README.md references:
- `fslmer_univariate.TEMPLATE.json` - Parameter details
- `meta_pipeline.TEMPLATE.json` - Pipeline options
- `docs/CONFIG_QUICK_REFERENCE.md` - Quick syntax
- `docs/fslmer_univariate_examples.md` - Analysis examples

### docs/CONFIG_QUICK_REFERENCE.md references:
- `configs/*.TEMPLATE.json` - Complete references
- `configs/README.md` - Detailed guide
- `docs/common_functions_reference.md` - Functions

### docs/CONFIGURATION_GUIDE.md references:
- `configs/README.md` - Main guide
- `*.TEMPLATE.json` - All parameters
- `REFACTORING_SUMMARY.md` - Architecture

## 🎯 Common Paths

### Path 1: Quick Single Analysis
```
1. Read: docs/CONFIG_QUICK_REFERENCE.md
2. Copy: fslmer_univariate.minimal.json
3. Edit: Your paths
4. Run: bash scripts/run_fslmer_univariate.sh --config your.json
```

### Path 2: Learn All Options
```
1. Read: configs/README.md
2. Study: fslmer_univariate.TEMPLATE.json
3. Reference: docs/CONFIG_QUICK_REFERENCE.md
4. Practice: Modify example configs
```

### Path 3: Full Pipeline
```
1. Read: configs/README.md (Meta-Pipeline section)
2. Copy: meta_pipeline.minimal.json
3. Reference: meta_pipeline.TEMPLATE.json
4. Run: bash scripts/run_meta_pipeline.sh --config your.json
```

### Path 4: Advanced Customization
```
1. Study: All TEMPLATE.json files
2. Read: docs/common_functions_reference.md
3. Check: REFACTORING_SUMMARY.md
4. Extend: Create custom scripts
```

## 📦 File Sizes & Lines

```
Configs:
- fslmer_univariate.TEMPLATE.json:   ~150 lines (heavily commented)
- fslmer_univariate.minimal.json:     14 lines (essential only)
- meta_pipeline.TEMPLATE.json:       ~250 lines (comprehensive)
- meta_pipeline.minimal.json:          9 lines (bare minimum)

Documentation:
- configs/README.md:                 ~600 lines (complete guide)
- docs/CONFIG_QUICK_REFERENCE.md:    ~200 lines (reference card)
- docs/CONFIGURATION_GUIDE.md:       ~350 lines (overview)
- docs/common_functions_reference.md: ~400 lines (functions)
```

## ✅ Validation Status

All new config files validated:
- ✅ `fslmer_univariate.minimal.json` - Valid JSON
- ✅ `meta_pipeline.minimal.json` - Valid JSON
- ✅ All TEMPLATE files - Valid JSON with comments as keys
- ✅ Syntax checked with `python -m json.tool`

## 🚀 Getting Started Checklist

- [ ] Read `configs/README.md` introduction
- [ ] Review `docs/CONFIG_QUICK_REFERENCE.md`
- [ ] Copy appropriate minimal.json file
- [ ] Edit paths in your config
- [ ] Validate JSON syntax
- [ ] Run with --dry-run first
- [ ] Execute actual analysis
- [ ] Check outputs

## 📞 Need Help?

1. **Syntax questions:** `docs/CONFIG_QUICK_REFERENCE.md`
2. **Parameter details:** `*.TEMPLATE.json` files
3. **Usage examples:** `configs/README.md`
4. **Functions:** `docs/common_functions_reference.md`
5. **Architecture:** `REFACTORING_SUMMARY.md`
6. **Script help:** `--help` flags

---

**Last Updated:** October 12, 2025
**Templates Created:** 4 configs + 3 documentation files
**Status:** All validated ✅
