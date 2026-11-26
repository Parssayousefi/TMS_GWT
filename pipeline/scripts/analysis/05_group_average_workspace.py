#!/usr/bin/env python3
"""
Group-Level Analysis: Identify Synergistic Workspace
=====================================================

This script replicates the analysis from Luppi et al. (2024) to:
1. Compute group-average synergy and redundancy matrices
2. Identify the synergistic workspace based on node strength
3. Save results for subsequent gateway/broadcaster analysis
4. Optionally compare workspaces across multiple datasets

Supports:
- Single dataset:   --dataset london_lsd
- Multiple datasets: --dataset london_lsd berlin_ketamine --compare

Based on Methods section: "Brain parcellation" and "Gradient of redundancy-to-synergy"
"""

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
import yaml
import sys
import argparse
from scipy import stats as scipy_stats

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from analysis_loader import (
    load_baseline_group_data,
    load_group_all_measures,
    get_control_task_name
)
from dataset_handler import load_dataset_handler

# Analysis parameters
N_ROIS = 454  # Schaefer 400 + Tian 54


# =============================================================================
# Core computation functions (Luppi et al. 2024 methodology)
# =============================================================================

def compute_node_strength(matrix):
    """
    Compute node strength (sum of connections for each node).

    Following Luppi et al.: "we calculated the nodal strength of each brain
    region as the sum of all its interactions"

    NaN-aware: PhiID pairs involving zero-variance ROIs are stored as NaN
    (not 0) so they are excluded via nansum rather than biasing rankings
    toward the low end. Rows that are entirely NaN (ROI with no valid edges
    in any subject) yield NaN strength so the ROI drops out of rankings
    downstream rather than getting a spurious rank of 1.
    """
    all_nan = np.all(np.isnan(matrix), axis=1)
    strength = np.nansum(matrix, axis=1)
    strength[all_nan] = np.nan
    return strength


def identify_workspace(synergy_strength, redundancy_strength):
    """
    Identify synergistic workspace regions.

    Step 1 of Luppi et al. (2024)'s two-step operationalization:
      Step 1 (here): rank the 454 ROIs by **node strength** on synergy vs
        redundancy — workspace = ROIs where synergy rank > redundancy rank
        (pure sign-of-rank-difference test, no percentile threshold).
      Step 2 (in `06_identify_gateways_broadcasters.py`): within the workspace
        only, rank ROIs by **participation coefficient** using the 8-module
        partition (7 Yeo + 1 subcortical) and split into gateways (positive
        sign) vs broadcasters (negative sign).

    Luppi et al. 2023: "Regions with positive gradient values (meaning they
    rank higher in synergy than redundancy relative to other brain regions)
    are identified as belonging to the synergistic global workspace"

    Note: The gradient is based on RELATIVE rankings, not absolute values.
    A positive gradient does not necessarily mean synergy > redundancy in
    absolute terms; it means the region's balance favours synergy relative
    to the rest of the brain.

    Parameters
    ----------
    synergy_strength : array (454,)
    redundancy_strength : array (454,)

    Returns
    -------
    workspace_mask : array (454,) of bool
    gradient : array (454,) — synergy rank minus redundancy rank
    synergy_ranks : array (454,)
    redundancy_ranks : array (454,)
    """
    # Rank regions by strength (higher strength = higher rank).
    # NaN-safe: rankdata with nan_policy='omit' leaves NaN strengths as NaN rank,
    # ranking only the finite ROIs. Argsort would otherwise put NaN at the end
    # and give zero-variance ROIs a spuriously high rank.
    synergy_ranks = scipy_stats.rankdata(synergy_strength, nan_policy='omit')
    redundancy_ranks = scipy_stats.rankdata(redundancy_strength, nan_policy='omit')

    # Gradient: synergy rank - redundancy rank (NaN - NaN = NaN for excluded ROIs)
    gradient = synergy_ranks - redundancy_ranks

    # Workspace = regions where synergy rank > redundancy rank (gradient > 0).
    # NaN comparisons are False, so NaN ROIs are correctly excluded from workspace.
    workspace_mask = gradient > 0

    return workspace_mask, gradient, synergy_ranks, redundancy_ranks


# =============================================================================
# Single-dataset workspace computation
# =============================================================================

def compute_workspace_for_dataset(dataset_name, project_dir=None, task=None, synergy_suffix='synergy'):
    """
    Compute the synergistic workspace for one dataset's specified condition.

    Parameters
    ----------
    dataset_name : str
        Dataset name
    project_dir : Path, optional
        Project directory
    task : str, optional
        Task/condition to analyze. If None, uses control task.
    synergy_suffix : str
        File suffix for synergy matrices (default: 'synergy').
        Use 'synergy_total' for total source synergy (Eq. 5).

    Returns a dict with all results needed for saving and comparison.
    """
    handler = load_dataset_handler(dataset_name, project_dir=project_dir)

    # Use provided task or default to control task
    if task is None:
        task = get_control_task_name(dataset_name)

    # Load group data for specified task
    results_dir = handler.project_dir / 'results' / 'phiid_matrices'
    all_matrices, loaded_subjects = load_group_all_measures(
        dataset_name,
        task,
        results_dir=results_dir,
        synergy_suffix=synergy_suffix
    )

    # Group averaging: arithmetic mean across subjects (Luppi 2024).
    # PhiID atoms are MI-based information-theoretic quantities — Fisher z / atanh
    # are inappropriate (they apply to correlations, not MI). Do NOT z-transform.
    # NaN-aware: zero-variance ROI pairs are stored as NaN per-subject; nanmean
    # averages only subjects with valid values at each pair.
    group_synergy = np.nanmean(all_matrices['synergy'], axis=0)
    group_redundancy = np.nanmean(all_matrices['redundancy'], axis=0)
    group_phi = np.nanmean(all_matrices['phi'], axis=0)
    group_phi_r = np.nanmean(all_matrices['phi_r'], axis=0)

    # Node strengths
    synergy_strength = compute_node_strength(group_synergy)
    redundancy_strength = compute_node_strength(group_redundancy)

    # Workspace identification
    workspace_mask, gradient, synergy_ranks, redundancy_ranks = identify_workspace(
        synergy_strength, redundancy_strength
    )

    return {
        'dataset': dataset_name,
        'control_task': task,  # Note: kept key name for compatibility, but can be any task
        'handler': handler,
        'subjects': loaded_subjects,
        'all_matrices': all_matrices,
        'group_synergy': group_synergy,
        'group_redundancy': group_redundancy,
        'group_phi': group_phi,
        'group_phi_r': group_phi_r,
        'synergy_strength': synergy_strength,
        'redundancy_strength': redundancy_strength,
        'synergy_ranks': synergy_ranks,
        'redundancy_ranks': redundancy_ranks,
        'workspace_mask': workspace_mask,
        'gradient': gradient,
    }


def save_single_dataset_results(result, base_results_dir, synergy_label=''):
    """
    Save workspace results for a single dataset.
    Output directory: base_results_dir/{dataset}/{control_task}[_total_synergy]/

    Parameters
    ----------
    synergy_label : str
        Appended to directory and file prefix (e.g. '_total_synergy').
    """
    dataset = result['dataset']
    control_task = result['control_task']
    dir_name = f"{control_task}{synergy_label}"
    analysis_dir = base_results_dir / dataset / dir_name
    analysis_dir.mkdir(exist_ok=True, parents=True)

    prefix = f"{dataset}_{control_task}{synergy_label}"

    # Group-average matrices
    np.save(analysis_dir / f'{prefix}_group_avg_synergy.npy', result['group_synergy'])
    np.save(analysis_dir / f'{prefix}_group_avg_redundancy.npy', result['group_redundancy'])
    np.save(analysis_dir / f'{prefix}_group_avg_phi.npy', result['group_phi'])
    np.save(analysis_dir / f'{prefix}_group_avg_phi_r.npy', result['group_phi_r'])

    # Node strengths and ranks
    np.save(analysis_dir / f'{prefix}_synergy_node_strength.npy', result['synergy_strength'])
    np.save(analysis_dir / f'{prefix}_redundancy_node_strength.npy', result['redundancy_strength'])
    np.save(analysis_dir / f'{prefix}_synergy_ranks.npy', result['synergy_ranks'])
    np.save(analysis_dir / f'{prefix}_redundancy_ranks.npy', result['redundancy_ranks'])

    # Workspace identification
    np.save(analysis_dir / f'{prefix}_workspace_mask.npy', result['workspace_mask'])
    np.save(analysis_dir / f'{prefix}_synergy_redundancy_gradient.npy', result['gradient'])

    # Metadata
    metadata = {
        'dataset': dataset,
        'task': control_task,
        'n_subjects': len(result['subjects']),
        'subjects': result['subjects'],
        'n_rois': N_ROIS,
        'n_workspace_regions': int(result['workspace_mask'].sum()),
        'workspace_method': 'gradient-based (synergy rank > redundancy rank)',
        'reference': 'Luppi et al. (2023, 2024)',
        'analysis_date': pd.Timestamp.now().isoformat(),
        'output_directory': str(analysis_dir),
        'file_prefix': prefix,
    }
    with open(analysis_dir / f'{prefix}_workspace_metadata.yaml', 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)

    # Summary CSV - ROI uses 1-454 to match atlas convention (see CLAUDE.md)
    summary_df = pd.DataFrame({
        'ROI': np.arange(1, N_ROIS + 1),
        'synergy_strength': result['synergy_strength'],
        'redundancy_strength': result['redundancy_strength'],
        'synergy_rank': result['synergy_ranks'],
        'redundancy_rank': result['redundancy_ranks'],
        'syn_red_gradient': result['gradient'],
        'workspace_member': result['workspace_mask'],
    })
    summary_df.to_csv(analysis_dir / f'{prefix}_workspace_summary.csv', index=False)

    return analysis_dir


def print_single_dataset_report(result):
    """Print a console report for one dataset."""
    dataset = result['dataset']
    control_task = result['control_task']
    n_subj = len(result['subjects'])
    n_ws = int(result['workspace_mask'].sum())

    lines = []
    lines.append("=" * 70)
    lines.append(f"WORKSPACE: {dataset} (control={control_task}, N={n_subj})")
    lines.append("=" * 70)

    lines.append(f"\nGroup-average synergy:")
    lines.append(f"  Mean: {np.nanmean(result['group_synergy']):.6f}, "
                 f"Range: [{np.nanmin(result['group_synergy']):.6f}, {np.nanmax(result['group_synergy']):.6f}]")
    lines.append(f"Group-average redundancy:")
    lines.append(f"  Mean: {np.nanmean(result['group_redundancy']):.6f}, "
                 f"Range: [{np.nanmin(result['group_redundancy']):.6f}, {np.nanmax(result['group_redundancy']):.6f}]")

    # Node-strength rankings: NaN (zero-variance ROIs) ignored via nan_policy='omit'
    syn_order = np.argsort(result['synergy_strength'])  # NaN sorts to end
    red_order = np.argsort(result['redundancy_strength'])
    # Drop trailing NaN entries so "top 5" reflects real ROIs
    syn_finite = syn_order[np.isfinite(result['synergy_strength'][syn_order])]
    red_finite = red_order[np.isfinite(result['redundancy_strength'][red_order])]
    lines.append(f"\nSynergy node strength — top 5: {syn_finite[-5:][::-1]}")
    lines.append(f"Redundancy node strength — top 5: {red_finite[-5:][::-1]}")

    lines.append(f"\nWorkspace: {n_ws} regions ({100*n_ws/N_ROIS:.1f}% of {N_ROIS})")
    lines.append(f"Gradient range: [{np.nanmin(result['gradient'])}, {np.nanmax(result['gradient'])}]")
    grad_order = np.argsort(result['gradient'])
    grad_finite = grad_order[np.isfinite(result['gradient'][grad_order])]
    lines.append(f"Most synergistic (top 5): {grad_finite[-5:][::-1]}")
    lines.append(f"Most redundant (bottom 5): {grad_finite[:5]}")

    report = '\n'.join(lines)
    print(report)
    return report


# =============================================================================
# Cross-dataset comparison
# =============================================================================

def dice_coefficient(mask_a, mask_b):
    """
    Dice similarity coefficient: 2*|A & B| / (|A| + |B|).
    Ranges from 0 (no overlap) to 1 (identical).
    """
    intersection = np.sum(mask_a & mask_b)
    return 2.0 * intersection / (np.sum(mask_a) + np.sum(mask_b))


def _compute_workspace_from_matrices(synergy_matrices, redundancy_matrices):
    """
    Given subject-level matrices, compute group workspace.
    Used by the permutation test.

    Parameters
    ----------
    synergy_matrices : array (n_subjects, 454, 454)
    redundancy_matrices : array (n_subjects, 454, 454)

    Returns
    -------
    workspace_mask : array (454,) of bool
    gradient : array (454,)
    """
    group_syn = np.nanmean(synergy_matrices, axis=0)
    group_red = np.nanmean(redundancy_matrices, axis=0)
    syn_str = compute_node_strength(group_syn)
    red_str = compute_node_strength(group_red)
    mask, grad, _, _ = identify_workspace(syn_str, red_str)
    return mask, grad


def permutation_test_dice(result_a, result_b, n_permutations=10000, seed=42):
    """
    Permutation test for workspace overlap between two datasets.

    Null hypothesis: the two datasets' workspaces are no more similar than
    expected if subjects were randomly assigned to two groups of the same
    sizes.

    Procedure:
    1. Pool all control subjects from both datasets.
    2. For each permutation, randomly assign pooled subjects into two groups
       matching the original dataset sizes.
    3. Recompute workspace for each shuffled group.
    4. Compute Dice between the two shuffled workspaces.
    5. p-value = proportion of null Dice values >= observed Dice.
    """
    rng = np.random.default_rng(seed)

    # Observed Dice
    observed_dice = dice_coefficient(result_a['workspace_mask'], result_b['workspace_mask'])

    # Pool subject-level matrices
    syn_a = result_a['all_matrices']['synergy']
    syn_b = result_b['all_matrices']['synergy']
    red_a = result_a['all_matrices']['redundancy']
    red_b = result_b['all_matrices']['redundancy']

    pooled_syn = np.concatenate([syn_a, syn_b], axis=0)
    pooled_red = np.concatenate([red_a, red_b], axis=0)

    n_a = syn_a.shape[0]
    n_total = pooled_syn.shape[0]

    null_dice = np.zeros(n_permutations)
    for i in range(n_permutations):
        perm = rng.permutation(n_total)
        idx_a = perm[:n_a]
        idx_b = perm[n_a:]

        mask_a, _ = _compute_workspace_from_matrices(pooled_syn[idx_a], pooled_red[idx_a])
        mask_b, _ = _compute_workspace_from_matrices(pooled_syn[idx_b], pooled_red[idx_b])

        null_dice[i] = dice_coefficient(mask_a, mask_b)

    p_value = np.mean(null_dice >= observed_dice)

    return observed_dice, p_value, null_dice


def compare_workspaces(results, base_results_dir, n_permutations=10000):
    """
    Compare workspace masks and gradients across multiple datasets.

    Computes for every pair of datasets:
    - Dice coefficient of workspace masks
    - Spearman correlation of gradient vectors
    - Permutation test of Dice (tests whether overlap exceeds chance)

    Also computes the conjunction workspace (regions in workspace for ALL datasets).

    Saves results to: base_results_dir/baseline_workspace_comparison/
    """
    comparison_dir = base_results_dir / 'baseline_workspace_comparison'
    comparison_dir.mkdir(exist_ok=True, parents=True)

    dataset_names = [r['dataset'] for r in results]
    n_datasets = len(results)

    print("\n" + "=" * 70)
    print("CROSS-DATASET WORKSPACE COMPARISON")
    print("=" * 70)
    print(f"Datasets: {', '.join(dataset_names)}")
    print(f"Permutations: {n_permutations}")

    # --- Pairwise comparisons ---
    pairwise_stats = []
    for (i, j) in combinations(range(n_datasets), 2):
        ra, rb = results[i], results[j]
        name_a, name_b = ra['dataset'], rb['dataset']

        print(f"\n--- {name_a} vs {name_b} ---")

        # Dice
        observed_dice, perm_p, null_dice = permutation_test_dice(
            ra, rb, n_permutations=n_permutations
        )
        print(f"  Dice coefficient: {observed_dice:.4f}")
        print(f"  Permutation p-value: {perm_p:.4f} ({n_permutations} permutations)")

        # Spearman correlation of gradients
        spearman_rho, spearman_p = scipy_stats.spearmanr(ra['gradient'], rb['gradient'])
        print(f"  Spearman rho (gradient): {spearman_rho:.4f} (p={spearman_p:.2e})")

        # Overlap counts
        overlap = np.sum(ra['workspace_mask'] & rb['workspace_mask'])
        union = np.sum(ra['workspace_mask'] | rb['workspace_mask'])
        only_a = np.sum(ra['workspace_mask'] & ~rb['workspace_mask'])
        only_b = np.sum(~ra['workspace_mask'] & rb['workspace_mask'])
        print(f"  Overlap: {overlap} ROIs, Union: {union} ROIs")
        print(f"  Only in {name_a}: {only_a}, Only in {name_b}: {only_b}")

        pairwise_stats.append({
            'dataset_a': name_a,
            'dataset_b': name_b,
            'n_subjects_a': len(ra['subjects']),
            'n_subjects_b': len(rb['subjects']),
            'dice_coefficient': float(observed_dice),
            'dice_permutation_p': float(perm_p),
            'n_permutations': n_permutations,
            'spearman_rho': float(spearman_rho),
            'spearman_p': float(spearman_p),
            'n_overlap': int(overlap),
            'n_union': int(union),
            'n_only_a': int(only_a),
            'n_only_b': int(only_b),
            'workspace_size_a': int(ra['workspace_mask'].sum()),
            'workspace_size_b': int(rb['workspace_mask'].sum()),
        })

        # Save null distribution for this pair
        np.save(
            comparison_dir / f'null_dice_{name_a}_vs_{name_b}.npy',
            null_dice
        )

    # --- Conjunction workspace (intersection of all) ---
    all_masks = np.stack([r['workspace_mask'] for r in results], axis=0)
    conjunction_mask = np.all(all_masks, axis=0)
    n_conjunction = int(conjunction_mask.sum())

    print(f"\n--- Conjunction workspace (all {n_datasets} datasets) ---")
    print(f"  Regions in workspace across ALL datasets: {n_conjunction} ({100*n_conjunction/N_ROIS:.1f}%)")

    # --- Pooled workspace (combined subject-level data) ---
    print(f"\n--- Pooled workspace (all subjects combined) ---")

    # Pool all subject-level matrices
    all_synergy = [r['all_matrices']['synergy'] for r in results]
    all_redundancy = [r['all_matrices']['redundancy'] for r in results]
    pooled_synergy_matrices = np.concatenate(all_synergy, axis=0)
    pooled_redundancy_matrices = np.concatenate(all_redundancy, axis=0)

    n_pooled_subjects = pooled_synergy_matrices.shape[0]
    print(f"  Total pooled subjects: {n_pooled_subjects}")

    # Compute pooled group averages (NaN-aware: zero-variance pairs excluded)
    pooled_group_synergy = np.nanmean(pooled_synergy_matrices, axis=0)
    pooled_group_redundancy = np.nanmean(pooled_redundancy_matrices, axis=0)

    # Compute pooled node strengths
    pooled_synergy_strength = compute_node_strength(pooled_group_synergy)
    pooled_redundancy_strength = compute_node_strength(pooled_group_redundancy)

    # Compute pooled workspace
    pooled_workspace_mask, pooled_gradient, pooled_syn_ranks, pooled_red_ranks = identify_workspace(
        pooled_synergy_strength, pooled_redundancy_strength
    )

    n_pooled_workspace = int(pooled_workspace_mask.sum())
    print(f"  Pooled workspace: {n_pooled_workspace} regions ({100*n_pooled_workspace/N_ROIS:.1f}%)")
    print(f"  Pooled gradient range: [{np.nanmin(pooled_gradient)}, {np.nanmax(pooled_gradient)}]")

    # Save pooled results (for use by script 06 with --pooled flag)
    np.save(comparison_dir / 'pooled_group_avg_synergy.npy', pooled_group_synergy)
    np.save(comparison_dir / 'pooled_group_avg_redundancy.npy', pooled_group_redundancy)
    np.save(comparison_dir / 'pooled_workspace_mask.npy', pooled_workspace_mask)
    np.save(comparison_dir / 'pooled_gradient.npy', pooled_gradient)
    np.save(comparison_dir / 'pooled_synergy_strength.npy', pooled_synergy_strength)
    np.save(comparison_dir / 'pooled_redundancy_strength.npy', pooled_redundancy_strength)
    np.save(comparison_dir / 'pooled_synergy_ranks.npy', pooled_syn_ranks)
    np.save(comparison_dir / 'pooled_redundancy_ranks.npy', pooled_red_ranks)
    print(f"  [OK] Saved pooled group-average matrices and workspace")

    # --- Per-ROI summary CSV ---
    # ROI column uses 1-454 to match atlas convention (see CLAUDE.md: ROI Indexing Convention)
    summary_data = {'ROI': np.arange(1, N_ROIS + 1)}
    for r in results:
        ds = r['dataset']
        summary_data[f'{ds}_gradient'] = r['gradient']
        summary_data[f'{ds}_workspace'] = r['workspace_mask']
        summary_data[f'{ds}_synergy_strength'] = r['synergy_strength']
        summary_data[f'{ds}_redundancy_strength'] = r['redundancy_strength']
    summary_data['conjunction_workspace'] = conjunction_mask
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(comparison_dir / 'workspace_comparison_summary.csv', index=False)

    # --- Save stats YAML ---
    stats_output = {
        'datasets': dataset_names,
        'n_rois': N_ROIS,
        'pairwise_comparisons': pairwise_stats,
        'conjunction_workspace_size': n_conjunction,
        'conjunction_workspace_pct': round(100 * n_conjunction / N_ROIS, 1),
        'pooled_n_subjects': n_pooled_subjects,
        'pooled_workspace_size': n_pooled_workspace,
        'pooled_workspace_pct': round(100 * n_pooled_workspace / N_ROIS, 1),
        'analysis_date': pd.Timestamp.now().isoformat(),
        'reference': 'Luppi et al. (2023, 2024)',
    }
    with open(comparison_dir / 'workspace_comparison_stats.yaml', 'w') as f:
        yaml.dump(stats_output, f, default_flow_style=False, sort_keys=False)

    # --- Save conjunction mask ---
    np.save(comparison_dir / 'conjunction_workspace_mask.npy', conjunction_mask)

    # --- Save metadata ---
    metadata = {
        'datasets_compared': dataset_names,
        'n_datasets': n_datasets,
        'per_dataset_info': [
            {
                'dataset': r['dataset'],
                'control_task': r['control_task'],
                'n_subjects': len(r['subjects']),
                'workspace_size': int(r['workspace_mask'].sum()),
            }
            for r in results
        ],
        'pooled_info': {
            'n_total_subjects': n_pooled_subjects,
            'workspace_size': n_pooled_workspace,
            'workspace_pct': round(100 * n_pooled_workspace / N_ROIS, 1),
        },
        'n_permutations': n_permutations,
        'analysis_date': pd.Timestamp.now().isoformat(),
        'output_directory': str(comparison_dir),
        'pooled_files': [
            'pooled_group_avg_synergy.npy',
            'pooled_group_avg_redundancy.npy',
            'pooled_workspace_mask.npy',
            'pooled_gradient.npy',
            'pooled_synergy_strength.npy',
            'pooled_redundancy_strength.npy',
        ],
    }
    with open(comparison_dir / 'workspace_comparison_metadata.yaml', 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False, sort_keys=False)

    # --- Markdown report ---
    md_lines = [
        "# Baseline Workspace Comparison",
        "",
        f"**Datasets**: {', '.join(dataset_names)}",
        f"**Date**: {pd.Timestamp.now().isoformat()}",
        f"**Method**: Luppi et al. (2023, 2024) — gradient-based workspace identification",
        "",
        "## Per-Dataset Summary",
        "",
        "| Dataset | Control task | N subjects | Workspace size |",
        "|---------|-------------|------------|----------------|",
    ]
    for r in results:
        n_ws = int(r['workspace_mask'].sum())
        md_lines.append(
            f"| {r['dataset']} | {r['control_task']} | {len(r['subjects'])} | "
            f"{n_ws} ({100*n_ws/N_ROIS:.1f}%) |"
        )

    md_lines.extend([
        "",
        "## Pairwise Comparisons",
        "",
        "| Pair | Dice | Dice perm. p | Spearman rho | Spearman p | Overlap | Union |",
        "|------|------|-------------|-------------|-----------|---------|-------|",
    ])
    for ps in pairwise_stats:
        md_lines.append(
            f"| {ps['dataset_a']} vs {ps['dataset_b']} | "
            f"{ps['dice_coefficient']:.4f} | {ps['dice_permutation_p']:.4f} | "
            f"{ps['spearman_rho']:.4f} | {ps['spearman_p']:.2e} | "
            f"{ps['n_overlap']} | {ps['n_union']} |"
        )

    md_lines.extend([
        "",
        "## Conjunction Workspace",
        "",
        f"Regions in workspace across **all** {n_datasets} datasets: "
        f"**{n_conjunction}** ({100*n_conjunction/N_ROIS:.1f}% of {N_ROIS} ROIs)",
        "",
        "## Pooled Workspace",
        "",
        f"Workspace computed from **pooled** subject data (N={n_pooled_subjects} total subjects):",
        f"**{n_pooled_workspace}** regions ({100*n_pooled_workspace/N_ROIS:.1f}% of {N_ROIS} ROIs)",
        "",
        "Pooled files saved for use with `06_identify_gateways_broadcasters.py --pooled`:",
        "- `pooled_group_avg_synergy.npy`",
        "- `pooled_group_avg_redundancy.npy`",
        "- `pooled_workspace_mask.npy`",
        "- `pooled_gradient.npy`",
        "",
        "## Interpretation",
        "",
        "- **Dice coefficient**: Spatial overlap of binary workspace masks (0=none, 1=identical)",
        "- **Permutation test**: Tests whether observed Dice exceeds chance (shuffling subjects across datasets)",
        "- **Spearman rho**: Correlation of the full 454-element synergy-redundancy gradient vectors",
        "- **Conjunction**: ROIs consistently identified as workspace members across all datasets",
        "- **Pooled**: Workspace computed by combining all subjects, giving more statistical power",
    ])

    with open(comparison_dir / 'workspace_comparison_results.md', 'w') as f:
        f.write('\n'.join(md_lines))

    print(f"\nComparison results saved to: {comparison_dir}")
    return comparison_dir


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Group-level workspace identification (Luppi et al. 2024)'
    )
    parser.add_argument(
        '--dataset', required=True, nargs='+',
        help='One or more dataset names (e.g., london_lsd berlin_ketamine)'
    )
    parser.add_argument(
        '--compare', action='store_true',
        help='Run cross-dataset workspace comparison (requires 2+ datasets)'
    )
    parser.add_argument(
        '--n_permutations', type=int, default=10000,
        help='Number of permutations for Dice significance test (default: 10000)'
    )
    parser.add_argument(
        '--project_dir', help='Project directory (optional)'
    )
    parser.add_argument(
        '--task', help='Task/condition to analyze (default: control task). E.g., LSD, ketamine, MDMA'
    )
    parser.add_argument(
        '--synergy_type', choices=['persistent', 'total'], default='persistent',
        help='Which synergy measure to use (default: persistent = sts atom). '
             '"total" loads *_synergy_total.npy files (Eq. 5, source synergy).'
    )

    args = parser.parse_args()

    datasets = args.dataset

    # Map synergy_type flag to file suffix and output label
    if args.synergy_type == 'total':
        synergy_suffix = 'synergy_total'
        synergy_label = '_total_synergy'
        print(f"NOTE: Using TOTAL source synergy (Eq. 5) files (*_synergy_total.npy)")
    else:
        synergy_suffix = 'synergy'
        synergy_label = ''

    # Resolve project dir from first dataset handler
    first_handler = load_dataset_handler(datasets[0], project_dir=args.project_dir)
    base_results_dir = first_handler.project_dir / 'results' / 'group_analysis'

    # --- Step 1: Compute workspace for each dataset ---
    all_results = []
    for ds_name in datasets:
        task_str = f" (task={args.task})" if args.task else ""
        print(f"\nComputing workspace for: {ds_name}{task_str}")
        try:
            result = compute_workspace_for_dataset(ds_name, project_dir=args.project_dir, task=args.task, synergy_suffix=synergy_suffix)
        except Exception as e:
            print(f"ERROR processing {ds_name}: {e}")
            print(f"  Check that PhiID matrices exist for {ds_name}")
            continue

        # Print report
        report = print_single_dataset_report(result)

        # Save per-dataset results
        out_dir = save_single_dataset_results(result, base_results_dir, synergy_label=synergy_label)
        print(f"\nSaved to: {out_dir}")

        # Save markdown report
        prefix = f"{ds_name}_{result['control_task']}{synergy_label}"
        with open(out_dir / f'{prefix}_workspace_results.md', 'w') as f:
            f.write(report)

        all_results.append(result)

    if len(all_results) == 0:
        print("\nNo datasets processed successfully. Exiting.")
        return

    # --- Step 2: Cross-dataset comparison ---
    if args.compare:
        if len(all_results) < 2:
            print("\nWARNING: --compare requires at least 2 datasets. Skipping comparison.")
        else:
            compare_workspaces(
                all_results, base_results_dir,
                n_permutations=args.n_permutations
            )

    # --- Summary ---
    print("\n" + "=" * 70)
    print("ALL DONE")
    print("=" * 70)
    for r in all_results:
        n_ws = int(r['workspace_mask'].sum())
        print(f"  {r['dataset']}: {n_ws} workspace regions "
              f"({100*n_ws/N_ROIS:.1f}%), N={len(r['subjects'])} subjects")
    if args.compare and len(all_results) >= 2:
        print(f"\n  Comparison saved to: {base_results_dir / 'baseline_workspace_comparison'}")

    print("\nNext steps:")
    for r in all_results:
        print(f"  python scripts/analysis/06_identify_gateways_broadcasters.py --dataset {r['dataset']}")


if __name__ == '__main__':
    main()
