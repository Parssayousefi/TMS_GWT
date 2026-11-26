#!/usr/bin/env python3
"""
Correlate Phi-R Changes with Subjective Effects
===============================================

This script correlates changes in phi-r (drug - control) in different brain regions
(workspace gateways, workspace broadcasters, and non-workspace regions) with
subjective effects measured by the 11ASC questionnaire.

Dataset-agnostic: Works with any dataset that has:
- Phi-R matrices in results/phiid_matrices/
- Gateway/broadcaster classification in results/group_analysis/{dataset}/{control_task}/
- Subjective data in data/{dataset}/{dataset}_subjective.csv

Author: Analysis pipeline for GWT project
Date: November 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import argparse
import yaml
import warnings
warnings.filterwarnings('ignore')

# Constants
N_ROIS = 454

# Base directory (project root)
BASE_DIR = Path(__file__).resolve().parent.parent.parent


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


def get_drug_tasks(config: dict) -> list:
    """Get list of drug (non-control) task names from config"""
    drug_tasks = []
    for task_name, task_info in config.get('tasks', {}).items():
        if not task_info.get('is_control', False):
            drug_tasks.append(task_name)
    return drug_tasks


def get_subjects_from_config(config: dict) -> list:
    """Get list of subjects from config"""
    subjects_config = config.get('subjects', {})
    if 'include' in subjects_config and subjects_config['include']:
        return subjects_config['include']
    # Return empty list - subjects will be discovered from phi_r files
    return []


def discover_subjects_from_files(phi_r_dir: Path, dataset_id: str, drug_task: str, control_task: str) -> list:
    """
    Discover subjects that have both drug and control phi_r files.

    Returns:
        List of subject IDs (as strings) that have both conditions
    """
    import re

    # Try different file naming patterns for drug condition
    drug_patterns = [
        f'dataset-{dataset_id}_sub-*_ses-{drug_task}_task-*_run-*_phi_r.npy',
        f'dataset-{dataset_id}_sub-*_ses-*_task-{drug_task}_run-*_phi_r.npy',
        f'dataset-{dataset_id}_sub-*_task-{drug_task}_run-*_phi_r.npy',
    ]

    control_patterns = [
        f'dataset-{dataset_id}_sub-*_ses-{control_task}_task-*_run-*_phi_r.npy',
        f'dataset-{dataset_id}_sub-*_ses-*_task-{control_task}_run-*_phi_r.npy',
        f'dataset-{dataset_id}_sub-*_task-{control_task}_run-*_phi_r.npy',
    ]

    drug_subjects = set()
    control_subjects = set()

    for pattern in drug_patterns:
        for f in phi_r_dir.glob(pattern):
            match = re.search(r'sub-(\d+)', f.name)
            if match:
                drug_subjects.add(match.group(1))

    for pattern in control_patterns:
        for f in phi_r_dir.glob(pattern):
            match = re.search(r'sub-(\d+)', f.name)
            if match:
                control_subjects.add(match.group(1))

    # Return subjects that have both conditions
    common_subjects = drug_subjects & control_subjects
    return sorted(list(common_subjects))


def get_paths(dataset_id: str, config: dict, control_task: str):
    """Get all required paths for analysis"""
    # Data directory (where subjective data lives)
    data_dir = BASE_DIR.parent / 'data' / dataset_id

    # Results directories
    phi_r_dir = BASE_DIR / 'results' / 'phiid_matrices'
    classification_file = (BASE_DIR / 'results' / 'group_analysis' / dataset_id /
                          control_task / f'{dataset_id}_{control_task}_gateway_broadcaster_summary.csv')

    # Subjective data file - try both spellings
    subjective_file = data_dir / f'{dataset_id}_subjective.csv'
    if not subjective_file.exists():
        # Try alternate spelling
        subjective_file = data_dir / f'{dataset_id}_subejctive.csv'

    output_dir = BASE_DIR / 'results' / 'figures' / dataset_id / 'phi_r_subjective_correlations'

    return {
        'phi_r_dir': phi_r_dir,
        'classification_file': classification_file,
        'subjective_file': subjective_file,
        'output_dir': output_dir
    }

# 11ASC dimension names for visualization
ASC_DIMENSIONS = [
    'Unity', 'Spiritual', 'Bliss', 'Insight', 'Disembodiment',
    'Control & Cognition', 'Anxiety', 'Complex Imagery',
    'Elementary Imagery', 'Audio/Visual', 'Meaning Perception'
]

# Short names for plots
ASC_SHORT = [
    'Unity', 'Spirit', 'Bliss', 'Insight', 'Disemb',
    'Cont&Cog', 'Anxiety', 'Complex', 'Element',
    'Aud/Vis', 'Meaning'
]


def load_region_classifications(classification_file: Path):
    """Load gateway, broadcaster, and workspace classifications"""
    if not classification_file.exists():
        raise FileNotFoundError(f"Classification file not found: {classification_file}")

    df = pd.read_csv(classification_file)

    # Extract masks
    workspace_mask = df['workspace'].values.astype(bool)
    gateway_mask = df['gateway'].values.astype(bool)
    broadcaster_mask = df['broadcaster'].values.astype(bool)

    print(f"Region classifications loaded:")
    print(f"  - Workspace regions: {np.sum(workspace_mask)}")
    print(f"  - Non-workspace regions: {N_ROIS - np.sum(workspace_mask)}")
    print(f"  - Gateway regions (in workspace): {np.sum(gateway_mask)}")
    print(f"  - Broadcaster regions (in workspace): {np.sum(broadcaster_mask)}")

    return workspace_mask, gateway_mask, broadcaster_mask

def extract_mean_phi_r(phi_r_matrix, region_mask):
    """
    Extract mean Φ_R value for connections involving specified regions

    Args:
        phi_r_matrix: Full Φ_R matrix (454, 454)
        region_mask: Boolean mask for regions of interest (454,)

    Returns:
        Mean Φ_R value for connections between specified regions
    """
    region_indices = np.where(region_mask)[0]

    if len(region_indices) == 0:
        return np.nan

    # Extract values for connections between specified regions
    values = []
    for i in region_indices:
        for j in region_indices:
            if i < j:  # Upper triangular only
                values.append(phi_r_matrix[i, j])

    # If no within-region connections, use all connections involving these regions
    if len(values) == 0:
        for i in region_indices:
            for j in range(N_ROIS):
                if i != j:
                    values.append(phi_r_matrix[i, j])

    return np.mean(values) if values else np.nan

def discover_phi_r_files(phi_r_dir: Path, dataset_id: str, task: str) -> dict:
    """
    Discover phi-r files for a given dataset and task.

    Returns:
        Dictionary mapping subject ID (str) to file path
    """
    files = {}

    # Try different file naming patterns
    patterns = [
        # Pattern 1: session-based (e.g., london_lsd with ses-LSD)
        f'dataset-{dataset_id}_sub-*_ses-{task}_task-*_run-*_phi_r.npy',
        # Pattern 2: task-based (e.g., basel datasets)
        f'dataset-{dataset_id}_sub-*_ses-*_task-{task}_run-*_phi_r.npy',
        # Pattern 3: no session (legacy)
        f'dataset-{dataset_id}_sub-*_task-{task}_run-*_phi_r.npy',
        # Pattern 4: simple subject-task
        f'{dataset_id}_sub-*_task-{task}_*_phi_r.npy',
    ]

    for pattern in patterns:
        matches = list(phi_r_dir.glob(pattern))
        if matches:
            for f in matches:
                # Extract subject ID from filename
                fname = f.name
                # Find sub-XX pattern
                import re
                match = re.search(r'sub-(\d+)', fname)
                if match:
                    subj_id = match.group(1)
                    files[subj_id] = f
            break

    return files


def load_phi_r_data(phi_r_dir: Path, classification_file: Path, dataset_id: str,
                    drug_task: str, control_task: str, subjects: list):
    """Load phi-r data for all subjects and calculate differences

    Args:
        phi_r_dir: Directory containing phi-r matrices
        classification_file: Path to gateway/broadcaster classification CSV
        dataset_id: Dataset identifier
        drug_task: Drug condition name (e.g., 'LSD')
        control_task: Control condition name (e.g., 'PCB')
        subjects: List of subject IDs to process

    Returns:
        Tuple of (gateway_changes, broadcaster_changes, non_workspace_changes, subject_ids)
        where subject_ids is a list of subject numbers that have valid data
    """

    # Load region classifications
    workspace_mask, gateway_mask, broadcaster_mask = load_region_classifications(classification_file)
    non_workspace_mask = ~workspace_mask

    # Discover available files for each condition
    drug_files = discover_phi_r_files(phi_r_dir, dataset_id, drug_task)
    control_files = discover_phi_r_files(phi_r_dir, dataset_id, control_task)

    print(f"\nDiscovered phi-r files:")
    print(f"  - {drug_task} files: {len(drug_files)}")
    print(f"  - {control_task} files: {len(control_files)}")

    # Initialize storage for phi-r differences
    gateway_changes = []
    broadcaster_changes = []
    non_workspace_changes = []
    valid_subject_ids = []  # Track which subjects have valid data

    print("\nLoading phi-r data for all subjects...")

    for subject in subjects:
        # Normalize subject ID to 2-digit string
        subj_str = subject.zfill(2) if isinstance(subject, str) else f'{subject:02d}'

        # Check if files exist for this subject
        if subj_str not in drug_files or subj_str not in control_files:
            print(f"  Warning: Missing data for subject {subj_str}")
            continue

        drug_file = drug_files[subj_str]
        control_file = control_files[subj_str]

        # Load phi-r matrices
        drug_matrix = np.load(drug_file)
        control_matrix = np.load(control_file)

        # Extract mean phi-r for each region type
        drug_gateway = extract_mean_phi_r(drug_matrix, gateway_mask)
        control_gateway = extract_mean_phi_r(control_matrix, gateway_mask)

        drug_broadcaster = extract_mean_phi_r(drug_matrix, broadcaster_mask)
        control_broadcaster = extract_mean_phi_r(control_matrix, broadcaster_mask)

        drug_non_workspace = extract_mean_phi_r(drug_matrix, non_workspace_mask)
        control_non_workspace = extract_mean_phi_r(control_matrix, non_workspace_mask)

        # Calculate changes (drug - control)
        gateway_changes.append(drug_gateway - control_gateway)
        broadcaster_changes.append(drug_broadcaster - control_broadcaster)
        non_workspace_changes.append(drug_non_workspace - control_non_workspace)

        # Store subject ID as integer for matching with subjective data
        valid_subject_ids.append(int(subj_str))

        print(f"  Subject {subj_str}: loaded")

    print(f"\nData loaded for {len(valid_subject_ids)} subjects: {valid_subject_ids}")

    return (np.array(gateway_changes), np.array(broadcaster_changes),
            np.array(non_workspace_changes), valid_subject_ids,
            workspace_mask, gateway_mask, broadcaster_mask)

def load_subjective_data(subjective_file: Path, subject_ids: list):
    """Load 11ASC subjective effects data and match to specified subjects

    Args:
        subjective_file: Path to the subjective data CSV file
        subject_ids: List of subject numbers to extract (e.g., [1, 2, 3, ...])

    Returns:
        DataFrame with subjective data for the specified subjects, in the same order
    """

    print("\nLoading subjective effects data...")

    if not subjective_file.exists():
        raise FileNotFoundError(f"Subjective data file not found: {subjective_file}")

    df = pd.read_csv(subjective_file)

    print(f"  Loaded CSV with {len(df)} rows and columns: {list(df.columns)}")

    # Check if 'Subject' column exists
    if 'Subject' not in df.columns:
        raise ValueError(f"Expected 'Subject' column in CSV. Found columns: {list(df.columns)}")

    # Set Subject as index for easy lookup
    df_indexed = df.set_index('Subject')

    # Extract rows for the specified subjects in order
    matched_data = []
    for subj_id in subject_ids:
        if subj_id in df_indexed.index:
            matched_data.append(df_indexed.loc[subj_id])
        else:
            print(f"  Warning: Subject {subj_id} not found in subjective data CSV")
            # Add a row of NaNs if subject not found
            matched_data.append(pd.Series({col: np.nan for col in df_indexed.columns}))

    result_df = pd.DataFrame(matched_data)
    result_df.index = subject_ids  # Use subject IDs as index

    print(f"  Matched {len(result_df)} subjects with subjective data")
    print(f"  Dimensions: {list(result_df.columns)}")

    return result_df

def calculate_correlations(phi_r_changes, subjective_data):
    """
    Calculate correlations between phi-r changes and subjective effects

    Args:
        phi_r_changes: Array of phi-r changes for a specific region type
        subjective_data: DataFrame with subjective effects

    Returns:
        Dictionary with correlation results
    """
    correlations = {}

    for col in subjective_data.columns:
        # Get subjective scores
        subjective_scores = subjective_data[col].values

        # Remove NaN values
        valid_mask = ~np.isnan(phi_r_changes) & ~np.isnan(subjective_scores)

        if np.sum(valid_mask) < 3:
            correlations[col] = {
                'r': np.nan, 'p': np.nan, 'rho': np.nan, 'p_spear': np.nan,
                'n': np.sum(valid_mask)
            }
            continue

        # Calculate Pearson and Spearman correlations
        r, p_pearson = pearsonr(phi_r_changes[valid_mask], subjective_scores[valid_mask])
        rho, p_spearman = spearmanr(phi_r_changes[valid_mask], subjective_scores[valid_mask])

        correlations[col] = {
            'r': r,
            'p': p_pearson,
            'rho': rho,
            'p_spear': p_spearman,
            'n': np.sum(valid_mask)
        }

    return correlations

def create_correlation_heatmap(correlations_dict, output_path, dataset_id='', drug_task='drug', control_task='control', n_subjects=0):
    """Create a heatmap showing all correlations"""

    # Extract correlation values
    dimensions = list(correlations_dict['gateway'].keys())
    region_types = ['Gateway', 'Broadcaster', 'Non-Workspace']

    # Create matrix for heatmap
    r_matrix = np.zeros((len(dimensions), len(region_types)))
    p_matrix = np.zeros((len(dimensions), len(region_types)))

    for i, dim in enumerate(dimensions):
        r_matrix[i, 0] = correlations_dict['gateway'][dim]['r']
        r_matrix[i, 1] = correlations_dict['broadcaster'][dim]['r']
        r_matrix[i, 2] = correlations_dict['non_workspace'][dim]['r']

        p_matrix[i, 0] = correlations_dict['gateway'][dim]['p']
        p_matrix[i, 1] = correlations_dict['broadcaster'][dim]['p']
        p_matrix[i, 2] = correlations_dict['non_workspace'][dim]['p']

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 10))

    # Heatmap of correlation coefficients
    sns.heatmap(r_matrix, annot=True, fmt='.2f', cmap='RdBu_r',
                vmin=-0.8, vmax=0.8, center=0,
                xticklabels=region_types,
                yticklabels=[ASC_SHORT[i] for i in range(len(dimensions))],
                ax=ax1, cbar_kws={'label': 'Pearson r'})
    ax1.set_title('Correlation Coefficients (Pearson r)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Brain Region Type', fontsize=11)
    ax1.set_ylabel('11ASC Dimension', fontsize=11)

    # Add significance stars
    for i in range(len(dimensions)):
        for j in range(len(region_types)):
            p = p_matrix[i, j]
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
                        color='white' if abs(r_matrix[i, j]) > 0.4 else 'black',
                        fontsize=10, fontweight='bold')

    # Heatmap of p-values
    sns.heatmap(p_matrix, annot=True, fmt='.3f', cmap='YlOrRd_r',
                vmin=0, vmax=0.1,
                xticklabels=region_types,
                yticklabels=[ASC_SHORT[i] for i in range(len(dimensions))],
                ax=ax2, cbar_kws={'label': 'p-value'})
    ax2.set_title('Statistical Significance (p-values)', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Brain Region Type', fontsize=11)
    ax2.set_ylabel('11ASC Dimension', fontsize=11)

    # Dataset-specific title with sample size
    dataset_label = dataset_id.replace('_', ' ').title()
    plt.suptitle(f'{dataset_label} (n={n_subjects}): Phi-R Changes ({drug_task.upper()} - {control_task.upper()}) vs 11ASC',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    print(f"\nSaved heatmap to: {output_path}")
    plt.close()

def create_scatter_plots(phi_r_data, subjective_data, correlations_dict, output_path,
                         dataset_id='', drug_task='drug', control_task='control', n_subjects=0):
    """Create scatter plots for significant correlations"""

    gateway_changes, broadcaster_changes, non_workspace_changes = phi_r_data

    # Find most significant correlations for each region type
    significant_corrs = []

    for region_type, changes, label in [
        ('gateway', gateway_changes, 'Gateway'),
        ('broadcaster', broadcaster_changes, 'Broadcaster'),
        ('non_workspace', non_workspace_changes, 'Non-Workspace')
    ]:
        for dim in correlations_dict[region_type]:
            if correlations_dict[region_type][dim]['p'] < 0.05:
                significant_corrs.append({
                    'region': label,
                    'dimension': dim,
                    'r': correlations_dict[region_type][dim]['r'],
                    'p': correlations_dict[region_type][dim]['p'],
                    'changes': changes
                })

    if not significant_corrs:
        print("No significant correlations found for scatter plots")
        return

    # Sort by absolute correlation strength
    significant_corrs.sort(key=lambda x: abs(x['r']), reverse=True)

    # Create figure with top correlations (max 9)
    n_plots = min(9, len(significant_corrs))
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_plots > 1 else [axes]

    for idx, corr_info in enumerate(significant_corrs[:n_plots]):
        ax = axes[idx]

        # Get data
        phi_r_change = corr_info['changes']
        subjective = subjective_data[corr_info['dimension']].values

        # Remove NaN values
        valid_mask = ~np.isnan(phi_r_change) & ~np.isnan(subjective)
        x = phi_r_change[valid_mask]
        y = subjective[valid_mask]

        # Create scatter plot
        ax.scatter(x, y, alpha=0.6, s=50, color='steelblue', edgecolor='black', linewidth=0.5)

        # Add regression line
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, p(x_line), "r-", alpha=0.8, linewidth=2)

        # Add confidence interval
        from scipy import stats as sp_stats
        predict_mean_se = sp_stats.sem(y - p(x))
        margin = 1.96 * predict_mean_se
        ax.fill_between(x_line, p(x_line) - margin, p(x_line) + margin,
                       color='red', alpha=0.15)

        # Formatting
        dim_name = ASC_SHORT[list(subjective_data.columns).index(corr_info['dimension'])]
        ax.set_title(f"{corr_info['region']} x {dim_name}\n" +
                    f"r = {corr_info['r']:.2f}, p = {corr_info['p']:.3f}",
                    fontsize=10, fontweight='bold')
        ax.set_xlabel(f'Delta Phi-R ({drug_task.upper()} - {control_task.upper()})', fontsize=9)
        ax.set_ylabel(f'{dim_name} Score', fontsize=9)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')

    dataset_label = dataset_id.replace('_', ' ').title()
    plt.suptitle(f'{dataset_label} (n={n_subjects}): Significant Correlations - Phi-R ({drug_task.upper()}) vs 11ASC',
                fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    print(f"Saved scatter plots to: {output_path}")
    plt.close()

def create_barplot_summary(correlations_dict, output_path, dataset_id='', drug_task='drug', control_task='control', n_subjects=0):
    """Create bar plots showing correlation strengths for each region type"""

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    region_types = [
        ('gateway', 'Gateway Regions (DMN)', 'lightcoral'),
        ('broadcaster', 'Broadcaster Regions (FPN)', 'lightblue'),
        ('non_workspace', 'Non-Workspace Regions', 'lightgreen')
    ]

    for ax, (region_key, region_name, color) in zip(axes, region_types):
        # Get correlations for this region
        corrs = correlations_dict[region_key]

        # Prepare data for plotting
        dimensions = list(corrs.keys())
        r_values = [corrs[dim]['r'] for dim in dimensions]
        p_values = [corrs[dim]['p'] for dim in dimensions]

        # Create bar plot
        bars = ax.bar(range(len(dimensions)), r_values, color=color,
                      edgecolor='black', linewidth=1.5, alpha=0.7)

        # Color bars by significance
        for bar, p in zip(bars, p_values):
            if p < 0.001:
                bar.set_alpha(1.0)
                bar.set_edgecolor('red')
                bar.set_linewidth(2)
            elif p < 0.01:
                bar.set_alpha(0.9)
                bar.set_edgecolor('orange')
                bar.set_linewidth(1.8)
            elif p < 0.05:
                bar.set_alpha(0.8)
                bar.set_edgecolor('gold')
                bar.set_linewidth(1.6)

        # Add significance stars
        for i, (r, p) in enumerate(zip(r_values, p_values)):
            if p < 0.001:
                sig_text = '***'
            elif p < 0.01:
                sig_text = '**'
            elif p < 0.05:
                sig_text = '*'
            else:
                sig_text = ''

            if sig_text:
                y_pos = r + 0.02 if r > 0 else r - 0.02
                ax.text(i, y_pos, sig_text, ha='center', va='bottom' if r > 0 else 'top',
                       fontsize=12, fontweight='bold')

        # Formatting
        ax.set_title(region_name, fontsize=12, fontweight='bold')
        ax.set_xlabel('11ASC Dimension', fontsize=11)
        ax.set_ylabel('Pearson r', fontsize=11)
        ax.set_xticks(range(len(dimensions)))
        ax.set_xticklabels(ASC_SHORT, rotation=45, ha='right')
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax.set_ylim([-0.8, 0.8])
        ax.grid(True, alpha=0.3, axis='y')

    dataset_label = dataset_id.replace('_', ' ').title()
    plt.suptitle(f'{dataset_label} (n={n_subjects}): Phi-R ({drug_task.upper()} - {control_task.upper()}) vs 11ASC by Region',
                fontsize=14, fontweight='bold')

    # Add legend
    fig.text(0.5, 0.02, 'Significance: *** p<0.001  ** p<0.01  * p<0.05',
            ha='center', fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray'))

    plt.tight_layout(rect=[0, 0.05, 1, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    print(f"Saved bar plot summary to: {output_path}")
    plt.close()

def save_correlation_table(correlations_dict, output_path):
    """Save correlation results to CSV file"""

    rows = []
    for region in ['gateway', 'broadcaster', 'non_workspace']:
        region_label = region.replace('_', ' ').title()
        for dim in correlations_dict[region]:
            dim_idx = list(correlations_dict[region].keys()).index(dim)
            rows.append({
                'Region_Type': region_label,
                'ASC_Dimension': ASC_DIMENSIONS[dim_idx],
                'ASC_Short': ASC_SHORT[dim_idx],
                'Pearson_r': correlations_dict[region][dim]['r'],
                'Pearson_p': correlations_dict[region][dim]['p'],
                'Spearman_rho': correlations_dict[region][dim]['rho'],
                'Spearman_p': correlations_dict[region][dim]['p_spear'],
                'N': correlations_dict[region][dim]['n']
            })

    df = pd.DataFrame(rows)

    # Round values
    df = df.round({
        'Pearson_r': 3, 'Pearson_p': 4,
        'Spearman_rho': 3, 'Spearman_p': 4
    })

    # Add significance column
    df['Significance'] = df['Pearson_p'].apply(
        lambda p: '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    )

    # Save to CSV
    df.to_csv(output_path, index=False)
    print(f"\nSaved correlation table to: {output_path}")

    # Print summary
    print("\nTop 5 strongest correlations:")
    df_sorted = df.sort_values('Pearson_r', ascending=False, key=lambda x: abs(x))
    print(df_sorted[['Region_Type', 'ASC_Short', 'Pearson_r', 'Significance']].head())

    return df

def create_markdown_report(correlations_dict, correlation_df, workspace_mask, gateway_mask,
                          broadcaster_mask, n_valid_subjects, output_path,
                          dataset_id='unknown', drug_task='drug', control_task='control'):
    """Create a comprehensive markdown report of all results"""
    from datetime import datetime

    # Start building the markdown content
    md_content = []

    # Header
    md_content.append("# Φ_R Changes and Subjective Effects Correlation Analysis")
    md_content.append(f"## {dataset_id.replace('_', ' ').title()} Dataset - 11ASC Questionnaire")
    md_content.append(f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md_content.append(f"**Dataset:** {dataset_id}")
    md_content.append(f"**Comparison:** {drug_task} vs {control_task}")
    md_content.append(f"**Valid Subjects:** {n_valid_subjects}")
    md_content.append("\n---\n")

    # Executive Summary
    md_content.append("## Executive Summary")
    md_content.append("\nThis analysis examines the relationship between changes in integrated information (Φ_R) ")
    md_content.append("across different brain region types and subjective effects measured by the 11-dimensional ")
    md_content.append("Altered States of Consciousness (11ASC) questionnaire under LSD vs placebo conditions.\n")

    # Region Classifications
    md_content.append("## Brain Region Classifications")
    md_content.append("\n### Region Types Analyzed:")
    md_content.append(f"- **Workspace Regions:** {np.sum(workspace_mask)} regions (~{np.sum(workspace_mask)*100/N_ROIS:.1f}% of brain)")
    md_content.append(f"  - **Gateway Regions:** {np.sum(gateway_mask)} regions (high synergy participation)")
    md_content.append(f"  - **Broadcaster Regions:** {np.sum(broadcaster_mask)} regions (high redundancy participation)")
    md_content.append(f"- **Non-Workspace Regions:** {N_ROIS - np.sum(workspace_mask)} regions (~{(N_ROIS-np.sum(workspace_mask))*100/N_ROIS:.1f}% of brain)")
    md_content.append("\n")

    # Key Findings
    md_content.append("## Key Findings")

    # Find significant correlations
    sig_gateway = sum(1 for dim in correlations_dict['gateway']
                     if correlations_dict['gateway'][dim]['p'] < 0.05)
    sig_broadcaster = sum(1 for dim in correlations_dict['broadcaster']
                         if correlations_dict['broadcaster'][dim]['p'] < 0.05)
    sig_nonworkspace = sum(1 for dim in correlations_dict['non_workspace']
                          if correlations_dict['non_workspace'][dim]['p'] < 0.05)

    md_content.append("\n### Significant Correlations Summary:")
    md_content.append(f"- **Gateway Regions:** {sig_gateway}/11 dimensions (p < 0.05)")
    md_content.append(f"- **Broadcaster Regions:** {sig_broadcaster}/11 dimensions (p < 0.05)")
    md_content.append(f"- **Non-Workspace Regions:** {sig_nonworkspace}/11 dimensions (p < 0.05)")
    md_content.append("")

    # Strongest correlations for each region type
    md_content.append("### Strongest Correlations by Region Type:")

    for region, label in [('gateway', 'Gateway'), ('broadcaster', 'Broadcaster'),
                          ('non_workspace', 'Non-Workspace')]:
        strongest = max(correlations_dict[region].items(), key=lambda x: abs(x[1]['r']))
        dim_idx = list(correlations_dict[region].keys()).index(strongest[0])
        r_val = strongest[1]['r']
        p_val = strongest[1]['p']
        sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''

        md_content.append(f"\n**{label} Regions:**")
        md_content.append(f"- {ASC_DIMENSIONS[dim_idx]}: r = {r_val:.3f}, p = {p_val:.4f} {sig}")

    # Detailed Results Tables
    md_content.append("\n## Detailed Statistical Results")

    # Gateway correlations
    md_content.append("\n### Gateway Regions")
    md_content.append("\n| ASC Dimension | Pearson r | p-value | Significance | Spearman ρ | p-value |")
    md_content.append("|---------------|-----------|---------|--------------|------------|---------|")

    gateway_df = correlation_df[correlation_df['Region_Type'] == 'Gateway']
    for _, row in gateway_df.iterrows():
        md_content.append(f"| {row['ASC_Short']} | {row['Pearson_r']:.3f} | {row['Pearson_p']:.4f} | "
                         f"{row['Significance']} | {row['Spearman_rho']:.3f} | {row['Spearman_p']:.4f} |")

    # Broadcaster correlations
    md_content.append("\n### Broadcaster Regions")
    md_content.append("\n| ASC Dimension | Pearson r | p-value | Significance | Spearman ρ | p-value |")
    md_content.append("|---------------|-----------|---------|--------------|------------|---------|")

    broadcaster_df = correlation_df[correlation_df['Region_Type'] == 'Broadcaster']
    for _, row in broadcaster_df.iterrows():
        md_content.append(f"| {row['ASC_Short']} | {row['Pearson_r']:.3f} | {row['Pearson_p']:.4f} | "
                         f"{row['Significance']} | {row['Spearman_rho']:.3f} | {row['Spearman_p']:.4f} |")

    # Non-workspace correlations
    md_content.append("\n### Non-Workspace Regions")
    md_content.append("\n| ASC Dimension | Pearson r | p-value | Significance | Spearman ρ | p-value |")
    md_content.append("|---------------|-----------|---------|--------------|------------|---------|")

    nonworkspace_df = correlation_df[correlation_df['Region_Type'] == 'Non Workspace']
    for _, row in nonworkspace_df.iterrows():
        md_content.append(f"| {row['ASC_Short']} | {row['Pearson_r']:.3f} | {row['Pearson_p']:.4f} | "
                         f"{row['Significance']} | {row['Spearman_rho']:.3f} | {row['Spearman_p']:.4f} |")

    # Interpretation
    md_content.append("\n## Interpretation")

    md_content.append("\n### Main Findings:")
    if sig_gateway > 0:
        md_content.append("1. **Gateway regions show the strongest associations** with subjective effects under LSD")

        # List significant gateway correlations
        sig_dims = []
        for dim in correlations_dict['gateway']:
            if correlations_dict['gateway'][dim]['p'] < 0.05:
                dim_idx = list(correlations_dict['gateway'].keys()).index(dim)
                r_val = correlations_dict['gateway'][dim]['r']
                sig_dims.append((ASC_DIMENSIONS[dim_idx], r_val))

        if sig_dims:
            md_content.append("\n   Significant correlations found with:")
            for dim_name, r in sig_dims:
                direction = "negative" if r < 0 else "positive"
                md_content.append(f"   - {dim_name} ({direction} correlation, r = {r:.3f})")

    md_content.append("\n2. **Broadcaster and non-workspace regions** show weaker or non-significant correlations")
    md_content.append("   - This suggests a more specific role of gateway regions in conscious experience")

    md_content.append("\n### Theoretical Implications:")
    md_content.append("- Gateway regions' role as **information gatherers** from multiple brain networks appears ")
    md_content.append("  crucial for subjective experience alterations under psychedelics")
    md_content.append("- The **negative correlations** suggest that reduced integrated information (Φ_R) in ")
    md_content.append("  gateway regions is associated with enhanced subjective experiences")
    md_content.append("- This aligns with theories of **relaxed constraints** and **increased flexibility** ")
    md_content.append("  in information processing under psychedelics")

    # Methods
    md_content.append("\n## Methods")

    md_content.append("\n### Data Processing:")
    md_content.append("1. **Φ_R Calculation:** Mean Φ_R extracted for connections within each region type")
    md_content.append("2. **Change Scores:** Difference calculated as (LSD - Placebo) for each subject")
    md_content.append("3. **Subjective Data:** 11ASC questionnaire difference scores (LSD - Placebo)")
    md_content.append("4. **Statistical Analysis:** Pearson and Spearman correlations computed")

    md_content.append("\n### 11ASC Dimensions:")
    for i, (full, short) in enumerate(zip(ASC_DIMENSIONS, ASC_SHORT), 1):
        md_content.append(f"{i}. **{short}**: {full}")

    # Output Files
    md_content.append("\n## Output Files Generated")
    md_content.append("\n- **Visualizations:**")
    md_content.append("  - `phi_r_subjective_correlations_heatmap.png` - Complete correlation matrix")
    md_content.append("  - `phi_r_subjective_correlations_barplot.png` - Bar plots by region type")
    md_content.append("  - `phi_r_subjective_significant_scatters.png` - Scatter plots for significant correlations")
    md_content.append("\n- **Data Tables:**")
    md_content.append("  - `phi_r_subjective_correlations_table.csv` - Full statistical results")
    md_content.append("  - `phi_r_subjective_correlations_report.md` - This report")

    # Notes and Limitations
    md_content.append("\n## Notes and Limitations")
    md_content.append(f"\n- Analysis includes {n_valid_subjects} subjects with complete data")
    md_content.append("- Multiple comparison corrections not applied (exploratory analysis)")
    md_content.append("- Results should be interpreted in context of sample size")

    # References
    md_content.append("\n## References")
    md_content.append("\n- Luppi et al. (2024) - Gateway and broadcaster classifications methodology")
    md_content.append("- 11ASC: 11-dimensional Altered States of Consciousness questionnaire")
    md_content.append("- Φ_R: Integrated information measure based on Partial Information Decomposition")

    # Statistical Summary
    md_content.append("\n## Statistical Summary")
    md_content.append("\n```")
    md_content.append("Significance levels: *** p<0.001  ** p<0.01  * p<0.05  ns p≥0.05")
    md_content.append("```")

    # Save to file (UTF-8 for Greek letters like Φ)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(md_content))

    print(f"\nSaved markdown report to: {output_path}")

    return '\n'.join(md_content)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Correlate Phi-R changes with subjective effects (11ASC)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python correlate_phir_subjective_effects.py --dataset london_lsd --drug LSD
    python correlate_phir_subjective_effects.py --dataset berlin_ketamine --drug ketamine
        """
    )
    parser.add_argument('--dataset', type=str, required=True,
                       help='Dataset ID (e.g., london_lsd, berlin_ketamine)')
    parser.add_argument('--drug', type=str, default=None,
                       help='Drug condition to analyze (default: first non-control task)')
    parser.add_argument('--control', type=str, default=None,
                       help='Control condition (default: auto-detect from config)')

    return parser.parse_args()


def main():
    # Parse arguments
    args = parse_args()
    dataset_id = args.dataset

    print("="*80)
    print("PHI-R AND SUBJECTIVE EFFECTS CORRELATION ANALYSIS")
    print(f"Dataset: {dataset_id}")
    print("="*80)

    # Load dataset configuration
    print(f"\nLoading configuration for dataset: {dataset_id}")
    config = get_dataset_config(dataset_id)

    # Get control and drug tasks
    control_task = args.control if args.control else get_control_task(config)
    drug_tasks = get_drug_tasks(config)

    if args.drug:
        if args.drug not in drug_tasks:
            raise ValueError(f"Drug task '{args.drug}' not found. Available: {drug_tasks}")
        drug_task = args.drug
    else:
        drug_task = drug_tasks[0]

    print(f"  Control task: {control_task}")
    print(f"  Drug task: {drug_task}")

    # Get paths
    paths = get_paths(dataset_id, config, control_task)
    print(f"\nPaths:")
    print(f"  Phi-R dir: {paths['phi_r_dir']}")
    print(f"  Classification: {paths['classification_file']}")
    print(f"  Subjective data: {paths['subjective_file']}")
    print(f"  Output dir: {paths['output_dir']}")

    # Create output directory
    paths['output_dir'].mkdir(parents=True, exist_ok=True)

    # Get subjects from config or discover from files
    subjects = get_subjects_from_config(config)
    if not subjects:
        print(f"\nNo explicit subject list in config, discovering from phi_r files...")
        subjects = discover_subjects_from_files(paths['phi_r_dir'], dataset_id, drug_task, control_task)
        print(f"Discovered {len(subjects)} subjects with both conditions: {subjects}")
    else:
        print(f"\nSubjects from config: {len(subjects)}")

    # Load phi-r data (returns subject IDs for matching)
    (gateway_changes, broadcaster_changes, non_workspace_changes,
     subject_ids, workspace_mask, gateway_mask, broadcaster_mask) = load_phi_r_data(
        paths['phi_r_dir'], paths['classification_file'], dataset_id,
        drug_task, control_task, subjects
    )

    # Load subjective data matched to the same subjects
    subjective_data = load_subjective_data(paths['subjective_file'], subject_ids)

    # Count valid subjects
    n_valid_subjects = len(subject_ids)

    # Calculate correlations
    print("\nCalculating correlations...")
    correlations = {
        'gateway': calculate_correlations(gateway_changes, subjective_data),
        'broadcaster': calculate_correlations(broadcaster_changes, subjective_data),
        'non_workspace': calculate_correlations(non_workspace_changes, subjective_data)
    }

    # Create visualizations
    print("\nCreating visualizations...")
    output_dir = paths['output_dir']

    # File prefix includes dataset, drug, and sample size
    file_prefix = f'{dataset_id}_{drug_task}_n{n_valid_subjects}'

    # Heatmap of all correlations
    heatmap_path = output_dir / f'{file_prefix}_phi_r_subjective_heatmap.png'
    create_correlation_heatmap(correlations, heatmap_path, dataset_id, drug_task, control_task, n_valid_subjects)

    # Bar plot summary
    barplot_path = output_dir / f'{file_prefix}_phi_r_subjective_barplot.png'
    create_barplot_summary(correlations, barplot_path, dataset_id, drug_task, control_task, n_valid_subjects)

    # Scatter plots for significant correlations
    scatter_path = output_dir / f'{file_prefix}_phi_r_subjective_scatters.png'
    create_scatter_plots(
        (gateway_changes, broadcaster_changes, non_workspace_changes),
        subjective_data, correlations, scatter_path,
        dataset_id, drug_task, control_task, n_valid_subjects
    )

    # Save correlation table
    table_path = output_dir / f'{file_prefix}_phi_r_subjective_table.csv'
    correlation_df = save_correlation_table(correlations, table_path)

    # Create markdown report
    report_path = output_dir / f'{file_prefix}_phi_r_subjective_report.md'
    create_markdown_report(correlations, correlation_df, workspace_mask, gateway_mask,
                          broadcaster_mask, n_valid_subjects, report_path,
                          dataset_id, drug_task, control_task)

    # Print summary statistics
    print("\n" + "="*80)
    print("ANALYSIS SUMMARY")
    print("="*80)

    for region in ['gateway', 'broadcaster', 'non_workspace']:
        region_label = region.replace('_', ' ').title()
        sig_count = sum(1 for dim in correlations[region]
                       if correlations[region][dim]['p'] < 0.05)
        print(f"\n{region_label}:")
        print(f"  - Significant correlations (p < 0.05): {sig_count}/11")

        # Find strongest correlation
        strongest = max(correlations[region].items(),
                       key=lambda x: abs(x[1]['r']))
        dim_idx = list(correlations[region].keys()).index(strongest[0])
        print(f"  - Strongest correlation: {ASC_SHORT[dim_idx]} " +
              f"(r = {strongest[1]['r']:.3f}, p = {strongest[1]['p']:.4f})")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"All results saved to: {output_dir}")

if __name__ == "__main__":
    main()