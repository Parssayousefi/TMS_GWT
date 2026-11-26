#!/usr/bin/env python3
"""
Network-Based Statistics (NBS) for Drug vs Control Comparison
==============================================================

This script replicates the NBS analysis from Luppi et al. (2024):
1. Compare Φ_R between drug and control conditions
2. Use permutation testing to control family-wise error
3. Identify connected components of significant changes

**UPDATED**: Now dataset-aware - works with any drug/control comparison

Based on Methods: "Network-based statistic" and "Testing for shared effects"

Key specifications from paper:
- F-threshold: 9 (two-sided, alpha=0.05)
- Permutations: 10,000 (we'll use 5,000 for computational feasibility)
- Within-subject comparison (paired design)

IMPORTANT NOTE (2025-10-13):
- Fixed critical bug: NBS now correctly uses EDGE COUNT (extent statistic)
  instead of NODE COUNT for component size calculation
- This matches Zalesky et al. (2010) original NBS implementation
- Previous results with critical value ~314 were based on node counts (incorrect)
- After this fix, expect critical value ~50-200 edges (correct)
"""

import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import sys
import argparse
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from dataset_handler import load_dataset_handler
from analysis_loader import load_drug_comparison_data, get_control_task_name

# Parameters
N_ROIS = 454
DEFAULT_N_PERMUTATIONS = 5000  # Paper uses 10,000
DEFAULT_F_THRESHOLD = 9  # From paper

def load_subject_phi_r(dataset_handler, subject, task, run, session=None):
    """Load Φ_R matrix for a single subject and task"""
    results_dir = dataset_handler.project_dir / 'results' / 'phiid_matrices'

    # Build filename
    filename_parts = [f'sub-{subject}']
    if session:
        filename_parts.append(f'ses-{session}')
    filename_parts.extend([f'task-{task}', f'run-{run}', 'phi_r.npy'])
    filename = '_'.join(filename_parts)

    phi_r_file = results_dir / filename

    if not phi_r_file.exists():
        return None

    return np.load(phi_r_file)

def compute_paired_f_statistic(drug_matrices, control_matrices):
    """
    Compute F-statistic for paired (within-subject) comparison

    Parameters
    ----------
    drug_matrices : array (n_subjects, n_rois, n_rois)
    control_matrices : array (n_subjects, n_rois, n_rois)

    Returns
    -------
    f_stats : array (n_rois, n_rois)
        F-statistic for each edge
    """
    n_subjects = drug_matrices.shape[0]

    # Paired differences: drug - control
    diff = drug_matrices - control_matrices  # (n_subjects, n_rois, n_rois)

    # Mean and variance of differences
    mean_diff = np.mean(diff, axis=0)
    var_diff = np.var(diff, axis=0, ddof=1)

    # F-statistic (from paired t-test)
    # F = t^2, where t = mean_diff / (std_diff / sqrt(n))
    # F = (mean_diff^2 * n) / var_diff
    with np.errstate(divide='ignore', invalid='ignore'):
        f_stats = (mean_diff ** 2 * n_subjects) / var_diff
        f_stats[np.isnan(f_stats)] = 0
        f_stats[np.isinf(f_stats)] = 0

    return f_stats

def threshold_and_find_components(f_matrix, threshold, return_size=True, statistic='extent'):
    """
    Threshold F-statistic matrix and find connected components

    Parameters
    ----------
    f_matrix : array (n_rois, n_rois)
        F-statistic matrix
    threshold : float
        Threshold for significance
    return_size : bool
        If True, return size of largest component; else return all component info
    statistic : str
        'extent' = count edges (default, standard NBS)
        'intensity' = sum of F-statistics (alternative NBS statistic)

    Returns
    -------
    max_component_size : int or float (if return_size=True)
        Size of largest connected component (in edges for 'extent', sum-F for 'intensity')
    OR
    n_components, labels, sizes : tuple (if return_size=False)
    """
    # Threshold matrix (keeping edges above threshold)
    thresholded = (np.abs(f_matrix) > threshold).astype(int)

    # Make symmetric (shouldn't be necessary, but just in case)
    thresholded = np.maximum(thresholded, thresholded.T)

    # Find connected components
    graph = csr_matrix(thresholded)
    n_components, labels = connected_components(
        csgraph=graph, directed=False, return_labels=True
    )

    if n_components == 0:
        if return_size:
            return 0
        else:
            return 0, labels, []

    # CRITICAL FIX: Calculate component statistics based on EDGES, not nodes
    # This is the correct implementation of NBS (Zalesky et al. 2010)
    component_stats = []
    for comp_idx in range(n_components):
        comp_mask = (labels == comp_idx)
        nodes_in_comp = np.where(comp_mask)[0]

        # Count edges or sum F-statistics in this component
        n_edges = 0
        sum_f = 0.0
        for i in nodes_in_comp:
            for j in nodes_in_comp:
                if i < j and np.abs(f_matrix[i, j]) > threshold:
                    n_edges += 1
                    sum_f += np.abs(f_matrix[i, j])

        if statistic == 'extent':
            component_stats.append(n_edges)  # Count of edges (standard NBS)
        elif statistic in ['intensity', 'mass']:
            component_stats.append(sum_f)  # Sum of F-statistics (alternative NBS)
        else:
            raise ValueError(f"Unknown statistic: {statistic}")

    if return_size:
        return max(component_stats) if component_stats else 0
    else:
        return n_components, labels, component_stats

def permutation_test_nbs(drug_matrices, control_matrices, f_threshold, n_permutations):
    """
    NBS permutation test for paired comparison

    Parameters
    ----------
    drug_matrices : array (n_subjects, n_rois, n_rois)
    control_matrices : array (n_subjects, n_rois, n_rois)
    f_threshold : float
        Threshold for forming connected components
    n_permutations : int
        Number of permutations

    Returns
    -------
    null_distribution : array (n_permutations,)
        Null distribution of maximum component sizes
    """
    n_subjects = drug_matrices.shape[0]

    print(f"Running {n_permutations} permutations...")
    null_distribution = np.zeros(n_permutations)

    for perm_idx in range(n_permutations):
        # Print progress every 100 permutations
        if perm_idx % 100 == 0:
            print(f"  Permutation {perm_idx}/{n_permutations} ({100*perm_idx/n_permutations:.1f}%)")
        # Random sign flips for paired design
        # For each subject, randomly swap drug and control
        sign_flips = np.random.choice([-1, 1], size=n_subjects)

        # Create permuted data
        permuted_drug = []
        permuted_control = []
        for subj_idx in range(n_subjects):
            if sign_flips[subj_idx] == 1:
                # Keep original
                permuted_drug.append(drug_matrices[subj_idx])
                permuted_control.append(control_matrices[subj_idx])
            else:
                # Swap
                permuted_drug.append(control_matrices[subj_idx])
                permuted_control.append(drug_matrices[subj_idx])

        permuted_drug = np.array(permuted_drug)
        permuted_control = np.array(permuted_control)

        # Compute F-statistics for permuted data
        perm_f_stats = compute_paired_f_statistic(permuted_drug, permuted_control)

        # Find largest component
        max_comp_size = threshold_and_find_components(
            perm_f_stats, f_threshold, return_size=True
        )

        null_distribution[perm_idx] = max_comp_size

    return null_distribution

def identify_significant_components(
    observed_f_stats,
    f_threshold,
    null_distribution,
    alpha=0.05,
    print_func=print
):
    """
    Identify significant connected components

    Parameters
    ----------
    observed_f_stats : array (n_rois, n_rois)
    f_threshold : float
    null_distribution : array (n_permutations,)
    alpha : float
        Significance level
    print_func : callable
        Function to use for printing (default: print)

    Returns
    -------
    sig_components : dict
        Information about significant components
    """
    # Find components in observed data
    n_components, labels, component_sizes = threshold_and_find_components(
        observed_f_stats, f_threshold, return_size=False
    )

    # Critical value from null distribution
    critical_value = np.percentile(null_distribution, (1 - alpha) * 100)

    print_func(f"\nNull distribution:")
    print_func(f"  Mean: {null_distribution.mean():.2f}")
    print_func(f"  Median: {np.median(null_distribution):.2f}")
    print_func(f"  95th percentile: {np.percentile(null_distribution, 95):.2f}")
    print_func(f"  Critical value (α={alpha}): {critical_value:.2f}")

    # Identify significant components
    sig_components = {}
    for comp_idx in range(n_components):
        comp_size = component_sizes[comp_idx]
        if comp_size >= critical_value:
            # Get nodes in this component
            nodes_in_component = np.where(labels == comp_idx)[0]

            # Get edges in this component
            comp_mask = (labels == comp_idx)
            edges_in_component = []
            for i in np.where(comp_mask)[0]:
                for j in np.where(comp_mask)[0]:
                    if i < j and np.abs(observed_f_stats[i, j]) > f_threshold:
                        edges_in_component.append((i, j, observed_f_stats[i, j]))

            sig_components[comp_idx] = {
                'size': comp_size,
                'nodes': nodes_in_component.tolist(),
                'n_nodes': len(nodes_in_component),
                'edges': edges_in_component,
                'n_edges': len(edges_in_component),
                'p_value': np.mean(null_distribution >= comp_size)
            }

    return sig_components, critical_value

def main():
    """Main NBS analysis pipeline"""
    parser = argparse.ArgumentParser(
        description='Network-Based Statistics - DATASET AWARE'
    )
    parser.add_argument('--dataset', required=True,
                       help='Dataset name (e.g., london_lsd, basel_mmm)')
    parser.add_argument('--drug_task', required=True,
                       help='Drug/experimental task name (e.g., LSD, MDMA)')
    parser.add_argument('--n_permutations', type=int, default=DEFAULT_N_PERMUTATIONS,
                       help=f'Number of permutations (default: {DEFAULT_N_PERMUTATIONS})')
    parser.add_argument('--f_threshold', type=float, default=DEFAULT_F_THRESHOLD,
                       help=f'F-statistic threshold (default: {DEFAULT_F_THRESHOLD})')
    parser.add_argument('--project_dir', help='Project directory (optional)')

    args = parser.parse_args()

    # Load dataset handler
    handler = load_dataset_handler(args.dataset, project_dir=args.project_dir)

    # Get control task
    control_task = get_control_task_name(args.dataset)

    # Setup output directory structure
    base_results_dir = handler.project_dir / 'results' / 'group_analysis'
    nbs_dir = base_results_dir / args.dataset / f'nbs_{args.drug_task}_vs_{control_task}'
    nbs_dir.mkdir(exist_ok=True, parents=True)

    # File prefix
    prefix = f"{args.dataset}_{args.drug_task}_vs_{control_task}"

    # Create markdown output file
    md_output = []

    def print_and_log(msg):
        """Print to console and log to markdown"""
        print(msg)
        md_output.append(msg)

    print_and_log("="*70)
    print_and_log(f"NETWORK-BASED STATISTICS (NBS): {args.drug_task} vs {control_task}")
    print_and_log("="*70)
    print_and_log(f"\nDataset: {handler.config['name']}")
    print_and_log(f"Replicating Luppi et al. (2024) NBS analysis")
    print_and_log(f"Comparing integrated information (Φ_R) between conditions")
    print_and_log(f"\nParameters:")
    print_and_log(f"  F-threshold: {args.f_threshold}")
    print_and_log(f"  Permutations: {args.n_permutations}")
    print_and_log(f"  Design: Paired (within-subject)")

    # -------------------------------------------------------------------------
    # Step 1: Load all subject data using analysis_loader
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 1: Loading subject data")
    print_and_log("="*70)

    try:
        control_matrices, drug_matrices, control_subjects, drug_subjects = load_drug_comparison_data(
            args.dataset,
            args.drug_task,
            results_dir=handler.project_dir / 'results' / 'phiid_matrices'
        )

        # Extract phi_r
        control_phi_r = control_matrices['phi_r']
        drug_phi_r = drug_matrices['phi_r']

    except Exception as e:
        print_and_log(f"\nERROR: Failed to load data: {e}")
        print_and_log(f"\nTroubleshooting:")
        print_and_log(f"  1. Check that ΦID matrices exist for both {control_task} and {args.drug_task}")
        print_and_log(f"  2. Verify subjects have data for both conditions")
        print_and_log(f"  3. Run: python scripts/utils/analysis_loader.py --dataset {args.dataset}")
        return

    n_subjects = len(control_subjects)
    print_and_log(f"\nSuccessfully loaded {n_subjects} subjects")
    print_and_log(f"Control subjects: {control_subjects}")
    print_and_log(f"Drug subjects: {drug_subjects}")

    print_and_log(f"\nData shapes:")
    print_and_log(f"  {control_task}: {control_phi_r.shape}")
    print_and_log(f"  {args.drug_task}: {drug_phi_r.shape}")

    # -------------------------------------------------------------------------
    # Step 2: Compute observed F-statistics
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 2: Computing observed F-statistics")
    print_and_log("="*70)

    observed_f_stats = compute_paired_f_statistic(drug_phi_r, control_phi_r)

    print_and_log(f"F-statistics computed:")
    print_and_log(f"  Shape: {observed_f_stats.shape}")
    print_and_log(f"  Mean: {observed_f_stats.mean():.4f}")
    print_and_log(f"  Max: {observed_f_stats.max():.4f}")
    print_and_log(f"  N edges > threshold ({args.f_threshold}): {(observed_f_stats > args.f_threshold).sum()}")

    # Mean differences
    mean_drug = np.mean(drug_phi_r, axis=0)
    mean_control = np.mean(control_phi_r, axis=0)
    mean_diff = mean_drug - mean_control

    print_and_log(f"\nMean Φ_R differences ({args.drug_task} - {control_task}):")
    print_and_log(f"  Mean: {mean_diff.mean():.6f}")
    print_and_log(f"  Increases: {(mean_diff > 0).sum()} edges")
    print_and_log(f"  Decreases: {(mean_diff < 0).sum()} edges")

    # Save observed statistics
    np.save(nbs_dir / f'{prefix}_observed_f_stats.npy', observed_f_stats)
    np.save(nbs_dir / f'{prefix}_mean_{args.drug_task}.npy', mean_drug)
    np.save(nbs_dir / f'{prefix}_mean_{control_task}.npy', mean_control)
    np.save(nbs_dir / f'{prefix}_mean_diff.npy', mean_diff)

    # -------------------------------------------------------------------------
    # Step 3: Run permutation test
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 3: Permutation testing")
    print_and_log("="*70)

    null_distribution = permutation_test_nbs(
        drug_phi_r, control_phi_r, args.f_threshold, args.n_permutations
    )

    np.save(nbs_dir / f'{prefix}_null_distribution.npy', null_distribution)

    # -------------------------------------------------------------------------
    # Step 4: Identify significant components
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 4: Identifying significant components")
    print_and_log("="*70)

    sig_components, critical_value = identify_significant_components(
        observed_f_stats, args.f_threshold, null_distribution, alpha=0.05, print_func=print_and_log
    )

    print_and_log(f"\nSignificant components found: {len(sig_components)}")
    for comp_idx, comp_info in sig_components.items():
        print_and_log(f"\nComponent {comp_idx}:")
        print_and_log(f"  Size: {comp_info['size']} edges (NBS statistic)")
        print_and_log(f"  N nodes: {comp_info['n_nodes']}")
        print_and_log(f"  N edges: {comp_info['n_edges']}")
        print_and_log(f"  p-value: {comp_info['p_value']:.4f}")

    # -------------------------------------------------------------------------
    # Step 5: Analyze direction of effects (increases vs decreases)
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 5: Analyzing direction of effects")
    print_and_log("="*70)

    # Separate increases and decreases
    sig_increases_mask = (observed_f_stats > args.f_threshold) & (mean_diff > 0)
    sig_decreases_mask = (observed_f_stats > args.f_threshold) & (mean_diff < 0)

    print_and_log(f"\nSignificant changes (F > {args.f_threshold}):")
    print_and_log(f"  Increases ({args.drug_task} > {control_task}): {sig_increases_mask.sum()} edges")
    print_and_log(f"  Decreases ({args.drug_task} < {control_task}): {sig_decreases_mask.sum()} edges")

    np.save(nbs_dir / f'{prefix}_sig_increases_mask.npy', sig_increases_mask)
    np.save(nbs_dir / f'{prefix}_sig_decreases_mask.npy', sig_decreases_mask)

    # -------------------------------------------------------------------------
    # Step 6: Save results
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 6: Saving results")
    print_and_log("="*70)

    # Save significant components as YAML
    # Convert numpy types to Python types for YAML serialization
    sig_components_serializable = {}
    for comp_idx, comp_info in sig_components.items():
        sig_components_serializable[int(comp_idx)] = {
            'size': int(comp_info['size']),
            'n_nodes': int(comp_info['n_nodes']),
            'n_edges': int(comp_info['n_edges']),
            'p_value': float(comp_info['p_value']),
            'nodes': [int(x) for x in comp_info['nodes']],
            'edges': [(int(i), int(j), float(f)) for i, j, f in comp_info['edges'][:100]]  # Limit for file size
        }

    with open(nbs_dir / f'{prefix}_significant_components.yaml', 'w') as f:
        yaml.dump(sig_components_serializable, f, default_flow_style=False)

    # Metadata
    metadata = {
        'dataset': args.dataset,
        'dataset_name': handler.config['name'],
        'drug_task': args.drug_task,
        'control_task': control_task,
        'n_subjects': n_subjects,
        'control_subjects': control_subjects,
        'drug_subjects': drug_subjects,
        'n_rois': N_ROIS,
        'f_threshold': float(args.f_threshold),
        'n_permutations': int(args.n_permutations),
        'alpha': 0.05,
        'critical_value': float(critical_value),
        'n_significant_components': len(sig_components),
        'n_sig_increases': int(sig_increases_mask.sum()),
        'n_sig_decreases': int(sig_decreases_mask.sum()),
        'analysis_date': pd.Timestamp.now().isoformat(),
        'output_directory': str(nbs_dir),
        'file_prefix': prefix
    }

    with open(nbs_dir / f'{prefix}_nbs_metadata.yaml', 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)

    print_and_log("✓ Saved all results")

    print_and_log("\n" + "="*70)
    print_and_log("NBS ANALYSIS COMPLETE")
    print_and_log("="*70)
    print_and_log(f"\nResults saved to: {nbs_dir}")
    print_and_log(f"File prefix: {prefix}_*")
    print_and_log(f"\nKey findings:")
    print_and_log(f"  - Significant components: {len(sig_components)}")
    print_and_log(f"  - Edges increased in {args.drug_task}: {sig_increases_mask.sum()}")
    print_and_log(f"  - Edges decreased in {args.drug_task}: {sig_decreases_mask.sum()}")
    print_and_log(f"\nNext steps:")
    print_and_log(f"  1. Run 08_analyze_regional_changes.py to identify affected regions")
    print_and_log(f"  2. Run visualization scripts")

    # Save markdown output
    md_file = nbs_dir / f'{prefix}_nbs_results.md'
    with open(md_file, 'w') as f:
        f.write('\n'.join(md_output))
    print(f"\n✓ Markdown output saved to: {md_file}")

if __name__ == '__main__':
    main()
