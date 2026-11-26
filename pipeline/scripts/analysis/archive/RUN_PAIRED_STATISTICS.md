# Paired Statistical Testing - Usage Guide

## Overview

The `07_paired_statistical_testing.py` script performs proper within-subject paired statistical testing for drug vs placebo comparisons. It accounts for the crossover design where each subject receives both conditions.

## What This Script Does

### Proper Paired Analysis
1. **Computes individual-subject metrics**
   - Participation coefficients for each subject in each condition
   - Gateway/broadcaster classifications per subject
   - Workspace identification per subject

2. **Performs paired statistical tests**
   - Paired t-tests for participation coefficient changes
   - McNemar's test for role changes (gateway/broadcaster)
   - Multiple comparison correction (FDR or Bonferroni)

3. **Generates comprehensive outputs**
   - ROI-level statistics with p-values, effect sizes, CIs
   - Network-level summaries
   - Significance masks
   - Markdown reports

### Key Difference from Previous Analysis

**Previous (06b_compare_conditions_gateways.py)**:
- Averages all subjects → single group brain
- Compares group averages (unpaired logic)
- Descriptive only, no significance testing

**New (07_paired_statistical_testing.py)**:
- Computes metrics for each subject individually
- Takes paired differences within subjects
- Tests if mean difference ≠ 0 (proper paired logic)
- Provides p-values, confidence intervals, effect sizes

---

## Usage

### Basic Usage

```bash
# LSD vs Placebo
python3 scripts/analysis/07_paired_statistical_testing.py \
    --dataset london_lsd \
    --drug_task LSD

# MDMA vs Placebo
python3 scripts/analysis/07_paired_statistical_testing.py \
    --dataset basel_mmm \
    --drug_task MDMA
```

### With Options

```bash
# Use Bonferroni correction instead of FDR
python3 scripts/analysis/07_paired_statistical_testing.py \
    --dataset basel_mmm \
    --drug_task MPH \
    --correction bonferroni

# Change workspace threshold
python3 scripts/analysis/07_paired_statistical_testing.py \
    --dataset london_lsd \
    --drug_task LSD \
    --threshold_percentile 95
```

### All Options

| Option | Description | Default |
|--------|-------------|---------|
| `--dataset` | Dataset name (required) | - |
| `--drug_task` | Drug task name (required) | - |
| `--threshold_percentile` | Workspace threshold | 90 |
| `--correction` | Multiple comparison method | fdr_bh |
| `--project_dir` | Project directory | Auto-detect |

**Correction methods**:
- `fdr_bh`: False Discovery Rate (Benjamini-Hochberg) - **recommended**
- `bonferroni`: Bonferroni correction (more conservative)
- `none`: No correction (not recommended)

---

## Run All Drugs

### Batch Script for All Basel MMM Drugs

```bash
#!/bin/bash
# Run paired statistics for all Basel MMM drugs

cd /home/s3648540/data_pi-michielvanelkm/Parsa/GWT

for drug in MDMA MPH moda; do
    echo "======================================"
    echo "Running paired statistics: $drug vs placebo"
    echo "======================================"

    python3 scripts/analysis/07_paired_statistical_testing.py \
        --dataset basel_mmm \
        --drug_task $drug

    echo ""
done

echo "All Basel MMM analyses complete!"
```

### Run LSD

```bash
python3 scripts/analysis/07_paired_statistical_testing.py \
    --dataset london_lsd \
    --drug_task LSD
```

---

## Output Files

Results are saved to: `results/group_analysis/{dataset}/paired_stats_{drug}_vs_{control}/`

### Files Generated

| File | Description |
|------|-------------|
| `*_synergy_pc_stats.csv` | ROI-level statistics for synergy PC changes |
| `*_redundancy_pc_stats.csv` | ROI-level statistics for redundancy PC changes |
| `*_synergy_pc_network_summary.csv` | Network-level summary for synergy PC |
| `*_redundancy_pc_network_summary.csv` | Network-level summary for redundancy PC |
| `*_metadata.yaml` | Analysis metadata and summary counts |
| `*_results.md` | Markdown report with full analysis log |

### ROI-Level Statistics Columns

**In `*_synergy_pc_stats.csv` and `*_redundancy_pc_stats.csv`**:

| Column | Description |
|--------|-------------|
| `ROI` | ROI ID (1-454) |
| `Network` | Network ID (0-7) |
| `Network_Name` | Network name (Visual, DMN, etc.) |
| `mean_diff` | Mean change: drug - control |
| `std_diff` | Standard deviation of differences |
| `se` | Standard error |
| `t_stat` | Paired t-statistic |
| `p_value` | Uncorrected p-value |
| `p_corrected` | FDR/Bonferroni corrected p-value |
| `significant` | True if p_corrected < 0.05 |
| `ci_lower` | 95% CI lower bound |
| `ci_upper` | 95% CI upper bound |
| `cohen_d` | Effect size (Cohen's d for paired data) |

### Network Summary Columns

| Column | Description |
|--------|-------------|
| `Network` | Network name |
| `Network_ID` | Network ID |
| `N_ROIs` | Number of ROIs in network |
| `Mean_Change` | Average change across network ROIs |
| `Median_Change` | Median change |
| `N_Significant` | Number of ROIs with significant changes |
| `Pct_Significant` | Percentage of ROIs significant |
| `Mean_Cohen_d` | Average effect size |
| `Max_Change` | Largest positive change |
| `Min_Change` | Largest negative change |
| `Min_p_corrected` | Smallest corrected p-value |

---

## Interpreting Results

### Participation Coefficient Changes

**Example**: Visual network under LSD
```
mean_diff = 0.036
p_corrected = 0.001
cohen_d = 1.2
ci_lower = 0.020, ci_upper = 0.052
```

**Interpretation**:
> "Visual network redundancy participation coefficient significantly increased
> under LSD (mean Δ = 0.036, 95% CI [0.020, 0.052], t(18) = 4.5, p_FDR = 0.001,
> Cohen's d = 1.2), indicating a large and reliable enhancement of visual
> broadcasting capacity."

### Effect Size Guidelines (Cohen's d)

| |d| | Interpretation |
|------|----------------|
| < 0.2 | Negligible |
| 0.2 - 0.5 | Small |
| 0.5 - 0.8 | Medium |
| > 0.8 | Large |

### Statistical Significance

- **p_corrected < 0.05**: Significant (survives multiple comparison correction)
- **p_corrected < 0.01**: Highly significant
- **p_corrected < 0.001**: Very highly significant

---

## Example Workflow

### Complete Analysis Pipeline

```bash
#!/bin/bash
# Complete gateway/broadcaster analysis with statistics

DATASET="basel_mmm"
DRUG="MDMA"

# 1. Run descriptive comparison (if not already done)
python3 scripts/analysis/06b_compare_conditions_gateways.py \
    --dataset $DATASET \
    --drug_task $DRUG

# 2. Run paired statistical testing
python3 scripts/analysis/07_paired_statistical_testing.py \
    --dataset $DATASET \
    --drug_task $DRUG

# 3. Check results
echo "Descriptive results:"
ls results/group_analysis/$DATASET/comparison_${DRUG}_vs_placebo/

echo ""
echo "Statistical results:"
ls results/group_analysis/$DATASET/paired_stats_${DRUG}_vs_placebo/

# 4. View network summaries
echo ""
echo "Network-level statistical summary:"
cat results/group_analysis/$DATASET/paired_stats_${DRUG}_vs_placebo/${DATASET}_${DRUG}_vs_placebo_redundancy_pc_network_summary.csv | column -t -s,
```

---

## Comparison: Descriptive vs Statistical

### What You Get from Both Scripts

| Metric | 06b (Descriptive) | 07 (Statistical) |
|--------|-------------------|------------------|
| Network counts | ✓ Yes | ✗ No |
| Role changes | ✓ Yes | ✗ No |
| PC changes | ✓ Mean only | ✓ Mean + CI + p-value |
| Significance testing | ✗ No | ✓ Yes |
| Effect sizes | ✗ No | ✓ Yes (Cohen's d) |
| Confidence intervals | ✗ No | ✓ Yes |
| Individual variability | ✗ No | ✓ Implicit in stats |

**Recommendation**: Run **both** scripts
- Use 06b for network count changes and role classification patterns
- Use 07 for statistical inference and confidence intervals

---

## Troubleshooting

### Error: "Subject mismatch between conditions"

**Cause**: Control and drug subjects don't match perfectly

**Solution**: Check your data with:
```bash
python3 scripts/utils/analysis_loader.py --dataset YOUR_DATASET
```

### Error: "Failed to load data"

**Cause**: PhiID matrices not computed for one or both conditions

**Solution**:
1. Check which matrices exist:
   ```bash
   ls results/phiid_matrices/ | grep YOUR_DRUG
   ls results/phiid_matrices/ | grep placebo
   ```
2. Run PhiID computation if needed

### Warning: "No significant changes found"

**Cause**: Either real (no effect) or underpowered (small n)

**Solution**:
- Check descriptive results (06b) - are changes present?
- If yes but not significant → may need larger sample or less conservative correction
- Try `--correction none` to see uncorrected p-values (for exploration only)

---

## Technical Details

### Paired T-Test Formula

For each ROI:
```
differences_i = drug_PC_i - control_PC_i  (for each subject i)
mean_diff = mean(differences)
SE = std(differences) / sqrt(n)
t = mean_diff / SE
df = n - 1
```

### Cohen's d for Paired Data

```
d = mean_diff / std_diff
```

Where `std_diff` is the standard deviation of the paired differences (not pooled SD).

### False Discovery Rate (FDR)

Controls the expected proportion of false positives among rejected hypotheses.

**Benefits**:
- More powerful than Bonferroni
- Appropriate for exploratory analysis
- Standard in neuroimaging

**Formula** (Benjamini-Hochberg):
1. Order p-values: p(1) ≤ p(2) ≤ ... ≤ p(m)
2. Find largest i where p(i) ≤ (i/m) * α
3. Reject all H(1) through H(i)

---

## Citation

When using this analysis, cite:

**Methods**:
> "We performed within-subject paired t-tests on individual-level participation
> coefficients, controlling the false discovery rate at α = 0.05 using the
> Benjamini-Hochberg procedure (Benjamini & Hochberg, 1995). Effect sizes were
> computed as Cohen's d for paired samples."

**References**:
- Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289-300.
- Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). Erlbaum.

---

## Questions?

For issues or questions:
1. Check the markdown report: `*_results.md`
2. Verify data pairing in metadata: `*_metadata.yaml`
3. Compare with descriptive results from 06b script

**Expected runtime**: 2-5 minutes per drug comparison
