#!/usr/bin/env python3
"""
Correlate Per-Subject Role Switches with Subjective Effects
============================================================

This script correlates per-subject gateway/broadcaster role switch metrics
with subjective effects measured by the 11ASC questionnaire.

Requires:
- Run 07_paired_statistical_testing.py first (generates per-subject metrics)
- Subjective data file in data/{dataset}/{dataset}_subjective.csv

Metrics correlated:
1. Gateway-to-broadcaster switches (count)
2. Broadcaster-to-gateway switches (count)
3. Net gateway change (count)
4. Mean PC gradient change in workspace (continuous)
5. Total role switches (count)

Usage:
    python 08_correlate_role_switches_subjective.py --dataset london_lsd --drug LSD
    python 08_correlate_role_switches_subjective.py --dataset berlin_ketamine --drug ketamine

Author: GWT Analysis Pipeline
Date: January 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
import argparse
import yaml
import warnings
warnings.filterwarnings('ignore')

# Constants
N_ROIS = 454

# Base directory (project root)
BASE_DIR = Path(__file__).resolve().parent.parent.parent


# 11ASC dimension names
ASC_DIMENSIONS = [
    'Unity', 'Spiritual', 'Bliss', 'Insight', 'Disembodiment',
    'Control & Cognition', 'Anxiety', 'Complex Imagery',
    'Elementary Imagery', 'Audio/Visual', 'Meaning Perception'
]

ASC_SHORT = [
    'Unity', 'Spirit', 'Bliss', 'Insight', 'Disemb',
    'Cont&Cog', 'Anxiety', 'Complex', 'Element',
    'Aud/Vis', 'Meaning'
]

# Metrics to analyze (from role switch CSV)
ROLE_SWITCH_METRICS = [
    'gateway_to_broadcaster',
    'broadcaster_to_gateway',
    'total_role_switches',
    'net_gateway_change',
    'net_broadcaster_change',
    'mean_synergy_pc_change_workspace',
    'mean_redundancy_pc_change_workspace',
    'pc_gradient_change_workspace'
]

METRIC_LABELS = {
    'gateway_to_broadcaster': 'G->B Switches',
    'broadcaster_to_gateway': 'B->G Switches',
    'total_role_switches': 'Total Switches',
    'net_gateway_change': 'Net Gateway Change',
    'net_broadcaster_change': 'Net Broadcaster Change',
    'mean_synergy_pc_change_workspace': 'Mean Syn PC Change',
    'mean_redundancy_pc_change_workspace': 'Mean Red PC Change',
    'pc_gradient_change_workspace': 'PC Gradient Change'
}


def get_dataset_config(dataset_id: str) -> dict:
    """Load dataset configuration from YAML file"""
    config_path = BASE_DIR / 'configs' / 'datasets' / f'{dataset_id}.yaml'
    if not config_path.exists():
        raise FileNotFoundError(f"Dataset config not found: {config_path}")

    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def get_control_task(config: dict) -> str:
    """Get the control/placebo task name from config"""
    for task_name, task_info in config.get('tasks', {}).items():
        if task_info.get('is_control', False):
            return task_name
    raise ValueError("No control task found in dataset config")


def load_role_switch_metrics(stats_dir: Path, prefix: str) -> pd.DataFrame:
    """
    Load per-subject role switch metrics from script 07 output.

    Parameters
    ----------
    stats_dir : Path
        Directory containing paired stats output
    prefix : str
        File prefix (e.g., 'london_lsd_LSD_vs_PCB')

    Returns
    -------
    metrics_df : pd.DataFrame
        Per-subject role switch metrics
    """
    metrics_file = stats_dir / f'{prefix}_role_switch_metrics.csv'

    if not metrics_file.exists():
        raise FileNotFoundError(
            f"Role switch metrics not found: {metrics_file}\n"
            f"Please run script 07 with --save_per_subject first:\n"
            f"  python scripts/analysis/07_paired_statistical_testing.py --dataset <dataset> --drug_task <drug>"
        )

    return pd.read_csv(metrics_file)


def load_subjective_data(data_dir: Path, dataset_id: str) -> pd.DataFrame:
    """Load 11ASC subjective effects data"""
    # Try both spellings
    subjective_file = data_dir / f'{dataset_id}_subjective.csv'
    if not subjective_file.exists():
        subjective_file = data_dir / f'{dataset_id}_subejctive.csv'

    if not subjective_file.exists():
        raise FileNotFoundError(
            f"Subjective data not found. Tried:\n"
            f"  - {data_dir / f'{dataset_id}_subjective.csv'}\n"
            f"  - {data_dir / f'{dataset_id}_subejctive.csv'}"
        )

    return pd.read_csv(subjective_file)


def match_subjects(role_switch_df: pd.DataFrame, subjective_df: pd.DataFrame) -> tuple:
    """
    Match subjects between role switch metrics and subjective data.

    Returns
    -------
    matched_metrics : pd.DataFrame
        Role switch metrics for matched subjects
    matched_subjective : pd.DataFrame
        Subjective data for matched subjects (same order)
    subject_ids : list
        List of matched subject IDs
    """
    # Get subject IDs from role switch data
    role_switch_subjects = set(role_switch_df['subject'].astype(int))

    # Get subject IDs from subjective data
    if 'Subject' in subjective_df.columns:
        subjective_subjects = set(subjective_df['Subject'].astype(int))
    else:
        raise ValueError("Subjective data must have 'Subject' column")

    # Find common subjects
    common_subjects = sorted(list(role_switch_subjects & subjective_subjects))

    if len(common_subjects) == 0:
        raise ValueError(
            f"No common subjects found!\n"
            f"Role switch subjects: {sorted(role_switch_subjects)}\n"
            f"Subjective subjects: {sorted(subjective_subjects)}"
        )

    # Filter and reorder both dataframes
    matched_metrics = role_switch_df[role_switch_df['subject'].astype(int).isin(common_subjects)].copy()
    matched_metrics = matched_metrics.sort_values('subject').reset_index(drop=True)

    subjective_indexed = subjective_df.set_index('Subject')
    matched_subjective = subjective_indexed.loc[common_subjects].reset_index()

    return matched_metrics, matched_subjective, common_subjects


def calculate_correlations(metrics_df: pd.DataFrame, subjective_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate correlations between role switch metrics and subjective effects.

    Returns
    -------
    correlations_df : pd.DataFrame
        Correlation results with columns:
        Metric, ASC_Dimension, ASC_Short, Pearson_r, Pearson_p, Spearman_rho, Spearman_p, N
    """
    results = []

    # Get subjective columns (excluding Subject)
    subjective_cols = [c for c in subjective_df.columns if c != 'Subject']

    for metric in ROLE_SWITCH_METRICS:
        if metric not in metrics_df.columns:
            continue

        metric_values = metrics_df[metric].values

        for i, asc_col in enumerate(subjective_cols):
            asc_values = subjective_df[asc_col].values

            # Remove NaN values
            valid_mask = ~np.isnan(metric_values) & ~np.isnan(asc_values)
            n_valid = np.sum(valid_mask)

            if n_valid < 3:
                results.append({
                    'Metric': metric,
                    'Metric_Label': METRIC_LABELS.get(metric, metric),
                    'ASC_Dimension': ASC_DIMENSIONS[i] if i < len(ASC_DIMENSIONS) else asc_col,
                    'ASC_Short': ASC_SHORT[i] if i < len(ASC_SHORT) else asc_col[:8],
                    'Pearson_r': np.nan,
                    'Pearson_p': np.nan,
                    'Spearman_rho': np.nan,
                    'Spearman_p': np.nan,
                    'N': n_valid
                })
                continue

            # Calculate correlations
            r, p_pearson = pearsonr(metric_values[valid_mask], asc_values[valid_mask])
            rho, p_spearman = spearmanr(metric_values[valid_mask], asc_values[valid_mask])

            results.append({
                'Metric': metric,
                'Metric_Label': METRIC_LABELS.get(metric, metric),
                'ASC_Dimension': ASC_DIMENSIONS[i] if i < len(ASC_DIMENSIONS) else asc_col,
                'ASC_Short': ASC_SHORT[i] if i < len(ASC_SHORT) else asc_col[:8],
                'Pearson_r': r,
                'Pearson_p': p_pearson,
                'Spearman_rho': rho,
                'Spearman_p': p_spearman,
                'N': n_valid
            })

    df = pd.DataFrame(results)

    # Add significance column
    df['Significance'] = df['Pearson_p'].apply(
        lambda p: '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    )

    return df


def create_correlation_heatmap(correlations_df: pd.DataFrame, output_path: Path,
                               dataset_id: str, drug_task: str, control_task: str, n_subjects: int):
    """Create a heatmap showing all correlations"""

    # Pivot for heatmap
    pivot_r = correlations_df.pivot(index='ASC_Short', columns='Metric_Label', values='Pearson_r')
    pivot_p = correlations_df.pivot(index='ASC_Short', columns='Metric_Label', values='Pearson_p')

    # Reorder columns
    col_order = [METRIC_LABELS[m] for m in ROLE_SWITCH_METRICS if METRIC_LABELS[m] in pivot_r.columns]
    pivot_r = pivot_r[col_order]
    pivot_p = pivot_p[col_order]

    # Reorder rows
    row_order = ASC_SHORT[:len(pivot_r)]
    pivot_r = pivot_r.reindex(row_order)
    pivot_p = pivot_p.reindex(row_order)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10))

    # Heatmap of correlation coefficients
    sns.heatmap(pivot_r, annot=True, fmt='.2f', cmap='RdBu_r',
                vmin=-0.8, vmax=0.8, center=0,
                ax=ax1, cbar_kws={'label': 'Pearson r'})
    ax1.set_title('Correlation Coefficients (Pearson r)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Role Switch Metric', fontsize=11)
    ax1.set_ylabel('11ASC Dimension', fontsize=11)
    plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')

    # Add significance stars
    for i, row in enumerate(row_order):
        for j, col in enumerate(col_order):
            p = pivot_p.loc[row, col]
            r = pivot_r.loc[row, col]
            if pd.notna(p):
                if p < 0.001:
                    sig_text = '***'
                elif p < 0.01:
                    sig_text = '**'
                elif p < 0.05:
                    sig_text = '*'
                else:
                    sig_text = ''

                if sig_text:
                    ax1.text(j + 0.5, i + 0.7, sig_text, ha='center', va='center',
                            color='white' if abs(r) > 0.4 else 'black',
                            fontsize=9, fontweight='bold')

    # Heatmap of p-values
    sns.heatmap(pivot_p, annot=True, fmt='.3f', cmap='YlOrRd_r',
                vmin=0, vmax=0.1,
                ax=ax2, cbar_kws={'label': 'p-value'})
    ax2.set_title('Statistical Significance (p-values)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Role Switch Metric', fontsize=11)
    ax2.set_ylabel('11ASC Dimension', fontsize=11)
    plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')

    # Title
    dataset_label = dataset_id.replace('_', ' ').title()
    plt.suptitle(f'{dataset_label} (n={n_subjects}): Role Switch Metrics ({drug_task.upper()}) vs 11ASC',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    print(f"\nSaved heatmap to: {output_path}")
    plt.close()


def create_scatter_plots(metrics_df: pd.DataFrame, subjective_df: pd.DataFrame,
                         correlations_df: pd.DataFrame, output_path: Path,
                         dataset_id: str, drug_task: str, control_task: str, n_subjects: int):
    """Create scatter plots for significant correlations"""

    # Find significant correlations
    sig_corrs = correlations_df[correlations_df['Pearson_p'] < 0.05].copy()
    sig_corrs = sig_corrs.sort_values('Pearson_r', key=lambda x: abs(x), ascending=False)

    if len(sig_corrs) == 0:
        print("No significant correlations found for scatter plots")
        return

    # Create figure
    n_plots = min(9, len(sig_corrs))
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_plots > 1 else [axes]

    subjective_cols = [c for c in subjective_df.columns if c != 'Subject']

    for idx, (_, row) in enumerate(sig_corrs.head(n_plots).iterrows()):
        ax = axes[idx]

        metric = row['Metric']
        asc_short = row['ASC_Short']

        # Find ASC column
        asc_idx = ASC_SHORT.index(asc_short) if asc_short in ASC_SHORT else 0
        asc_col = subjective_cols[asc_idx] if asc_idx < len(subjective_cols) else subjective_cols[0]

        x = metrics_df[metric].values
        y = subjective_df[asc_col].values

        # Remove NaN
        valid_mask = ~np.isnan(x) & ~np.isnan(y)
        x_valid = x[valid_mask]
        y_valid = y[valid_mask]

        # Scatter plot
        ax.scatter(x_valid, y_valid, alpha=0.6, s=50, color='steelblue',
                  edgecolor='black', linewidth=0.5)

        # Regression line
        z = np.polyfit(x_valid, y_valid, 1)
        p = np.poly1d(z)
        x_line = np.linspace(x_valid.min(), x_valid.max(), 100)
        ax.plot(x_line, p(x_line), "r-", alpha=0.8, linewidth=2)

        # Labels
        ax.set_title(f"{row['Metric_Label']} x {asc_short}\n" +
                    f"r = {row['Pearson_r']:.2f}, p = {row['Pearson_p']:.3f}",
                    fontsize=10, fontweight='bold')
        ax.set_xlabel(row['Metric_Label'], fontsize=9)
        ax.set_ylabel(f'{asc_short} Score', fontsize=9)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')

    dataset_label = dataset_id.replace('_', ' ').title()
    plt.suptitle(f'{dataset_label} (n={n_subjects}): Significant Correlations - Role Switches vs 11ASC',
                fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    print(f"Saved scatter plots to: {output_path}")
    plt.close()


def create_barplot_summary(correlations_df: pd.DataFrame, output_path: Path,
                           dataset_id: str, drug_task: str, control_task: str, n_subjects: int):
    """Create bar plots showing correlation strengths for each metric"""

    # Get unique metrics
    metrics = correlations_df['Metric'].unique()
    n_metrics = len(metrics)

    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    for idx, metric in enumerate(metrics):
        if idx >= len(axes):
            break

        ax = axes[idx]
        metric_data = correlations_df[correlations_df['Metric'] == metric]

        r_values = metric_data['Pearson_r'].values
        p_values = metric_data['Pearson_p'].values
        labels = metric_data['ASC_Short'].values

        # Bar colors based on significance
        colors = []
        for p in p_values:
            if p < 0.001:
                colors.append('darkred')
            elif p < 0.01:
                colors.append('red')
            elif p < 0.05:
                colors.append('orange')
            else:
                colors.append('lightgray')

        bars = ax.bar(range(len(labels)), r_values, color=colors, edgecolor='black', linewidth=1)

        ax.set_title(METRIC_LABELS.get(metric, metric), fontsize=11, fontweight='bold')
        ax.set_xlabel('11ASC Dimension', fontsize=9)
        ax.set_ylabel('Pearson r', fontsize=9)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax.set_ylim([-0.8, 0.8])
        ax.grid(True, alpha=0.3, axis='y')

    # Hide unused subplots
    for idx in range(n_metrics, len(axes)):
        axes[idx].axis('off')

    dataset_label = dataset_id.replace('_', ' ').title()
    plt.suptitle(f'{dataset_label} (n={n_subjects}): Role Switch Metrics ({drug_task.upper()}) vs 11ASC by Metric',
                fontsize=14, fontweight='bold')

    # Add legend
    fig.text(0.5, 0.02, 'Colors: darkred=p<0.001  red=p<0.01  orange=p<0.05  gray=ns',
            ha='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray'))

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    print(f"Saved bar plot summary to: {output_path}")
    plt.close()


def create_markdown_report(correlations_df: pd.DataFrame, metrics_df: pd.DataFrame,
                           n_subjects: int, output_path: Path,
                           dataset_id: str, drug_task: str, control_task: str):
    """Create a comprehensive markdown report"""
    from datetime import datetime

    md = []

    # Header
    md.append("# Role Switch Metrics and Subjective Effects Correlation Analysis")
    md.append(f"## {dataset_id.replace('_', ' ').title()} Dataset - 11ASC Questionnaire")
    md.append(f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md.append(f"**Dataset:** {dataset_id}")
    md.append(f"**Comparison:** {drug_task} vs {control_task}")
    md.append(f"**Valid Subjects:** {n_subjects}")
    md.append("\n---\n")

    # Executive Summary
    md.append("## Executive Summary")
    md.append("\nThis analysis examines the relationship between per-subject role switch metrics ")
    md.append("(gateway <-> broadcaster transitions) and subjective effects measured by the ")
    md.append("11-dimensional Altered States of Consciousness (11ASC) questionnaire.\n")

    # Role switch summary
    md.append("## Role Switch Metrics Summary")
    md.append("\n| Metric | Mean | Std | Min | Max |")
    md.append("|--------|------|-----|-----|-----|")

    for metric in ROLE_SWITCH_METRICS:
        if metric in metrics_df.columns:
            values = metrics_df[metric]
            md.append(f"| {METRIC_LABELS.get(metric, metric)} | {values.mean():.2f} | {values.std():.2f} | {values.min():.2f} | {values.max():.2f} |")

    # Significant correlations
    md.append("\n## Significant Correlations (p < 0.05)")

    sig_corrs = correlations_df[correlations_df['Pearson_p'] < 0.05].copy()
    sig_corrs = sig_corrs.sort_values('Pearson_p')

    if len(sig_corrs) > 0:
        md.append(f"\n**Found {len(sig_corrs)} significant correlations:**\n")
        md.append("| Metric | ASC Dimension | Pearson r | p-value | Sig |")
        md.append("|--------|---------------|-----------|---------|-----|")

        for _, row in sig_corrs.iterrows():
            md.append(f"| {row['Metric_Label']} | {row['ASC_Short']} | {row['Pearson_r']:.3f} | {row['Pearson_p']:.4f} | {row['Significance']} |")
    else:
        md.append("\nNo significant correlations found at p < 0.05")

    # Full correlation table
    md.append("\n## Full Correlation Table")
    md.append("\n| Metric | ASC | r | p | rho | p_spear | Sig |")
    md.append("|--------|-----|---|---|-----|---------|-----|")

    for _, row in correlations_df.iterrows():
        md.append(f"| {row['Metric_Label'][:15]} | {row['ASC_Short']} | {row['Pearson_r']:.2f} | {row['Pearson_p']:.3f} | {row['Spearman_rho']:.2f} | {row['Spearman_p']:.3f} | {row['Significance']} |")

    # Interpretation
    md.append("\n## Interpretation")
    md.append("\n### Key Questions Addressed:")
    md.append("1. Do subjects with more gateway->broadcaster switches report different subjective effects?")
    md.append("2. Is the PC gradient change (synergy vs redundancy shift) related to altered states?")
    md.append("3. Do subjects with net increases in gateways show different consciousness patterns?")

    # Methods
    md.append("\n## Methods")
    md.append("\n### Metrics Analyzed:")
    for metric in ROLE_SWITCH_METRICS:
        md.append(f"- **{METRIC_LABELS.get(metric, metric)}**: {metric}")

    md.append("\n### Statistical Analysis:")
    md.append("- Pearson and Spearman correlations computed")
    md.append("- Multiple comparison corrections not applied (exploratory analysis)")

    # Save
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(md))

    print(f"\nSaved markdown report to: {output_path}")


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Correlate role switch metrics with subjective effects (11ASC)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python 08_correlate_role_switches_subjective.py --dataset london_lsd --drug LSD
    python 08_correlate_role_switches_subjective.py --dataset berlin_ketamine --drug ketamine
        """
    )
    parser.add_argument('--dataset', type=str, required=True,
                       help='Dataset ID (e.g., london_lsd, berlin_ketamine)')
    parser.add_argument('--drug', type=str, required=True,
                       help='Drug condition name (e.g., LSD, ketamine)')
    parser.add_argument('--control', type=str, default=None,
                       help='Control condition (default: auto-detect from config)')

    return parser.parse_args()


def main():
    args = parse_args()
    dataset_id = args.dataset
    drug_task = args.drug

    print("=" * 80)
    print("ROLE SWITCH METRICS AND SUBJECTIVE EFFECTS CORRELATION ANALYSIS")
    print(f"Dataset: {dataset_id}")
    print(f"Drug: {drug_task}")
    print("=" * 80)

    # Load config
    config = get_dataset_config(dataset_id)
    control_task = args.control if args.control else get_control_task(config)

    print(f"\nControl task: {control_task}")

    # Setup paths
    stats_dir = BASE_DIR / 'results' / 'group_analysis' / dataset_id / f'paired_stats_{drug_task}_vs_{control_task}'
    data_dir = BASE_DIR.parent / 'data' / dataset_id
    output_dir = BASE_DIR / 'results' / 'figures' / dataset_id / 'role_switch_subjective_correlations'
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = f'{dataset_id}_{drug_task}_vs_{control_task}'

    print(f"\nLoading role switch metrics from: {stats_dir}")

    # Load role switch metrics
    metrics_df = load_role_switch_metrics(stats_dir, prefix)
    print(f"  Loaded metrics for {len(metrics_df)} subjects")

    # Load subjective data
    print(f"\nLoading subjective data from: {data_dir}")
    subjective_df = load_subjective_data(data_dir, dataset_id)
    print(f"  Loaded data for {len(subjective_df)} subjects")

    # Match subjects
    print("\nMatching subjects...")
    matched_metrics, matched_subjective, subject_ids = match_subjects(metrics_df, subjective_df)
    n_subjects = len(subject_ids)
    print(f"  Matched {n_subjects} subjects: {subject_ids}")

    # Calculate correlations
    print("\nCalculating correlations...")
    correlations_df = calculate_correlations(matched_metrics, matched_subjective)

    # File prefix for output
    file_prefix = f'{dataset_id}_{drug_task}_n{n_subjects}'

    # Save correlation table
    table_path = output_dir / f'{file_prefix}_role_switch_subjective_table.csv'
    correlations_df.to_csv(table_path, index=False)
    print(f"Saved correlation table to: {table_path}")

    # Create visualizations
    print("\nCreating visualizations...")

    # Heatmap
    heatmap_path = output_dir / f'{file_prefix}_role_switch_subjective_heatmap.png'
    create_correlation_heatmap(correlations_df, heatmap_path, dataset_id, drug_task, control_task, n_subjects)

    # Bar plots
    barplot_path = output_dir / f'{file_prefix}_role_switch_subjective_barplot.png'
    create_barplot_summary(correlations_df, barplot_path, dataset_id, drug_task, control_task, n_subjects)

    # Scatter plots
    scatter_path = output_dir / f'{file_prefix}_role_switch_subjective_scatters.png'
    create_scatter_plots(matched_metrics, matched_subjective, correlations_df, scatter_path,
                        dataset_id, drug_task, control_task, n_subjects)

    # Markdown report
    report_path = output_dir / f'{file_prefix}_role_switch_subjective_report.md'
    create_markdown_report(correlations_df, matched_metrics, n_subjects, report_path,
                          dataset_id, drug_task, control_task)

    # Print summary
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)

    sig_count = len(correlations_df[correlations_df['Pearson_p'] < 0.05])
    print(f"\nSignificant correlations (p < 0.05): {sig_count}")

    if sig_count > 0:
        print("\nTop significant correlations:")
        sig_df = correlations_df[correlations_df['Pearson_p'] < 0.05].sort_values('Pearson_p')
        for _, row in sig_df.head(5).iterrows():
            print(f"  {row['Metric_Label']} x {row['ASC_Short']}: r={row['Pearson_r']:.3f}, p={row['Pearson_p']:.4f}")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
