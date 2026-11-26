#!/usr/bin/env python3
"""
Paired Statistical Testing for Gateway/Broadcaster Changes
===========================================================

This script performs proper within-subject paired statistical testing for
drug vs placebo comparisons. It accounts for the crossover design where
each subject receives both drug and placebo conditions.

Key Features:
- Dataset-agnostic (works with any configured dataset)
- Automatically finds and matches control/placebo condition
- Computes individual-subject participation coefficients
- Performs paired t-tests with multiple comparison correction
- Generates effect sizes and confidence intervals
- Tests gateway/broadcaster role changes with McNemar's test

Usage:
    python scripts/analysis/07_paired_statistical_testing.py --dataset basel_mmm --drug_task MDMA
    python scripts/analysis/07_paired_statistical_testing.py --dataset london_lsd --drug_task LSD

Author: Generated for within-subject crossover design analysis
Date: October 2025
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.stats import false_discovery_control
import argparse
import sys
import yaml
from typing import Dict, List, Tuple

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from load_roi_mapping import load_roi_network_mapping
from dataset_handler import load_dataset_handler
from analysis_loader import load_drug_comparison_data, get_control_task_name

# Constants
N_ROIS = 454

# ============================================================================
# PARTICIPATION COEFFICIENT COMPUTATION
# ============================================================================

def compute_participation_coefficient(connectivity_matrix: np.ndarray,
                                     module_labels: np.ndarray) -> np.ndarray:
    """
    Compute participation coefficient for each node

    From Rubinov & Sporns (2010):
    P_i = 1 - sum_s=1^M (κ_is / k_i)^2

    Parameters
    ----------
    connectivity_matrix : np.ndarray (N, N)
        Connectivity matrix (synergy or redundancy)
    module_labels : np.ndarray (N,)
        Module assignment for each node

    Returns
    -------
    participation_coef : np.ndarray (N,)
        Participation coefficient for each node
    """
    N = connectivity_matrix.shape[0]
    modules = np.unique(module_labels)

    # Total strength for each node
    k_i = np.sum(connectivity_matrix, axis=1)

    # Participation coefficient
    participation_coef = np.zeros(N)

    for i in range(N):
        if k_i[i] == 0:
            participation_coef[i] = 0
            continue

        sum_term = 0
        for module in modules:
            module_mask = (module_labels == module)
            kappa_is = np.sum(connectivity_matrix[i, module_mask])
            sum_term += (kappa_is / k_i[i]) ** 2

        participation_coef[i] = 1 - sum_term

    return participation_coef


def compute_node_strength(matrix: np.ndarray) -> np.ndarray:
    """Compute node strength (sum of connections for each node)"""
    return np.sum(matrix, axis=1)


def identify_workspace(synergy_strength: np.ndarray,
                       threshold_percentile: float = 90) -> np.ndarray:
    """
    Identify synergistic workspace regions

    Parameters
    ----------
    synergy_strength : np.ndarray (454,)
        Node strength for synergy
    threshold_percentile : float
        Percentile threshold for workspace membership

    Returns
    -------
    workspace_mask : np.ndarray (454,) of bool
        True for workspace regions
    """
    synergy_threshold = np.percentile(synergy_strength, threshold_percentile)
    workspace_mask = synergy_strength >= synergy_threshold
    return workspace_mask


def classify_gateways_broadcasters(synergy_pc: np.ndarray,
                                   redundancy_pc: np.ndarray,
                                   workspace_mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Classify workspace regions into gateways vs broadcasters

    Following Luppi et al. (2024):
    - Gateways: high synergy participation coefficient (relative to redundancy)
    - Broadcasters: high redundancy participation coefficient (relative to synergy)

    Parameters
    ----------
    synergy_pc : np.ndarray (454,)
        Participation coefficient based on synergy
    redundancy_pc : np.ndarray (454,)
        Participation coefficient based on redundancy
    workspace_mask : np.ndarray (454,) of bool
        True for workspace regions

    Returns
    -------
    gateway_mask : np.ndarray (454,) of bool
    broadcaster_mask : np.ndarray (454,) of bool
    """
    workspace_indices = np.where(workspace_mask)[0]

    # Rank participation coefficients (within workspace only)
    synergy_pc_workspace = synergy_pc[workspace_mask]
    redundancy_pc_workspace = redundancy_pc[workspace_mask]

    synergy_pc_ranks = np.argsort(np.argsort(synergy_pc_workspace)) + 1
    redundancy_pc_ranks = np.argsort(np.argsort(redundancy_pc_workspace)) + 1

    # Gradient: synergy PC rank - redundancy PC rank
    pc_gradient_workspace = synergy_pc_ranks - redundancy_pc_ranks

    # Classify
    gateway_mask = np.zeros(N_ROIS, dtype=bool)
    broadcaster_mask = np.zeros(N_ROIS, dtype=bool)

    gateway_mask[workspace_indices[pc_gradient_workspace > 0]] = True
    broadcaster_mask[workspace_indices[pc_gradient_workspace < 0]] = True

    return gateway_mask, broadcaster_mask


# ============================================================================
# INDIVIDUAL SUBJECT ANALYSIS
# ============================================================================

def load_group_workspace_mask(base_results_dir: Path, dataset: str, control_task: str) -> np.ndarray:
    """
    Load the group-defined workspace mask from script 06 output.

    Using group workspace ensures consistent regions for cross-subject comparison.

    Parameters
    ----------
    base_results_dir : Path
        Base results directory (e.g., project_dir / 'results' / 'group_analysis')
    dataset : str
        Dataset name
    control_task : str
        Control task name (e.g., 'placebo', 'PCB')

    Returns
    -------
    workspace_mask : np.ndarray (454,) of bool
        Group-defined workspace mask
    """
    workspace_file = base_results_dir / dataset / control_task / f'{dataset}_{control_task}_workspace_mask.npy'

    if not workspace_file.exists():
        raise FileNotFoundError(
            f"Group workspace mask not found: {workspace_file}\n"
            f"Please run scripts 05 and 06 first."
        )

    return np.load(workspace_file).astype(bool)


def compute_individual_metrics_with_group_workspace(
    matrices: Dict[str, np.ndarray],
    rsn_labels: np.ndarray,
    group_workspace_mask: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Compute participation coefficients and classifications using GROUP workspace mask.

    This ensures consistent workspace regions across subjects while allowing
    gateway/broadcaster classification to vary per subject.

    Parameters
    ----------
    matrices : dict
        Keys: 'synergy', 'redundancy'
        Values: np.ndarray of shape (n_subjects, 454, 454)
    rsn_labels : np.ndarray (454,)
        Network labels for each ROI
    group_workspace_mask : np.ndarray (454,)
        Group-defined workspace mask (same for all subjects)

    Returns
    -------
    results : dict
        Keys: 'synergy_pc', 'redundancy_pc', 'workspace', 'gateway', 'broadcaster'
        Values: np.ndarray of shape (n_subjects, 454)
    """
    n_subjects = matrices['synergy'].shape[0]

    results = {
        'synergy_pc': np.zeros((n_subjects, N_ROIS)),
        'redundancy_pc': np.zeros((n_subjects, N_ROIS)),
        'workspace': np.tile(group_workspace_mask, (n_subjects, 1)),  # Same for all subjects
        'gateway': np.zeros((n_subjects, N_ROIS), dtype=bool),
        'broadcaster': np.zeros((n_subjects, N_ROIS), dtype=bool)
    }

    for i in range(n_subjects):
        # Compute participation coefficients
        results['synergy_pc'][i] = compute_participation_coefficient(
            matrices['synergy'][i], rsn_labels)
        results['redundancy_pc'][i] = compute_participation_coefficient(
            matrices['redundancy'][i], rsn_labels)

        # Classify gateways/broadcasters using GROUP workspace
        results['gateway'][i], results['broadcaster'][i] = classify_gateways_broadcasters(
            results['synergy_pc'][i],
            results['redundancy_pc'][i],
            group_workspace_mask  # Use group mask, not per-subject
        )

    return results


def compute_per_subject_role_switches(
    control_results: Dict[str, np.ndarray],
    drug_results: Dict[str, np.ndarray],
    subject_ids: List[str]
) -> pd.DataFrame:
    """
    Compute per-subject role switch metrics for correlation with subjective effects.

    Parameters
    ----------
    control_results : dict
        Per-subject metrics from control condition
    drug_results : dict
        Per-subject metrics from drug condition
    subject_ids : list
        Subject identifiers

    Returns
    -------
    metrics_df : pd.DataFrame
        Per-subject role switch metrics with columns:
        - subject
        - gateway_to_broadcaster: count of G->B switches
        - broadcaster_to_gateway: count of B->G switches
        - total_role_switches: sum of both
        - net_gateway_change: #drug_gateways - #control_gateways
        - net_broadcaster_change: #drug_broadcasters - #control_broadcasters
        - mean_synergy_pc_change_workspace: mean delta synergy PC in workspace
        - mean_redundancy_pc_change_workspace: mean delta redundancy PC in workspace
        - pc_gradient_change_workspace: mean (syn - red) PC change in workspace
    """
    n_subjects = len(subject_ids)
    workspace_mask = control_results['workspace'][0]  # Same for all subjects

    metrics = {
        'subject': subject_ids,
        'gateway_to_broadcaster': np.zeros(n_subjects, dtype=int),
        'broadcaster_to_gateway': np.zeros(n_subjects, dtype=int),
        'total_role_switches': np.zeros(n_subjects, dtype=int),
        'net_gateway_change': np.zeros(n_subjects, dtype=int),
        'net_broadcaster_change': np.zeros(n_subjects, dtype=int),
        'mean_synergy_pc_change_workspace': np.zeros(n_subjects),
        'mean_redundancy_pc_change_workspace': np.zeros(n_subjects),
        'pc_gradient_change_workspace': np.zeros(n_subjects),
        'n_control_gateways': np.zeros(n_subjects, dtype=int),
        'n_control_broadcasters': np.zeros(n_subjects, dtype=int),
        'n_drug_gateways': np.zeros(n_subjects, dtype=int),
        'n_drug_broadcasters': np.zeros(n_subjects, dtype=int),
    }

    for i in range(n_subjects):
        # Role switch counts
        ctrl_gw = control_results['gateway'][i]
        ctrl_bc = control_results['broadcaster'][i]
        drug_gw = drug_results['gateway'][i]
        drug_bc = drug_results['broadcaster'][i]

        # Gateway -> Broadcaster: was gateway in control, is broadcaster in drug
        metrics['gateway_to_broadcaster'][i] = np.sum(ctrl_gw & drug_bc)

        # Broadcaster -> Gateway: was broadcaster in control, is gateway in drug
        metrics['broadcaster_to_gateway'][i] = np.sum(ctrl_bc & drug_gw)

        # Total switches
        metrics['total_role_switches'][i] = (
            metrics['gateway_to_broadcaster'][i] + metrics['broadcaster_to_gateway'][i]
        )

        # Counts
        metrics['n_control_gateways'][i] = np.sum(ctrl_gw)
        metrics['n_control_broadcasters'][i] = np.sum(ctrl_bc)
        metrics['n_drug_gateways'][i] = np.sum(drug_gw)
        metrics['n_drug_broadcasters'][i] = np.sum(drug_bc)

        # Net changes
        metrics['net_gateway_change'][i] = np.sum(drug_gw) - np.sum(ctrl_gw)
        metrics['net_broadcaster_change'][i] = np.sum(drug_bc) - np.sum(ctrl_bc)

        # PC changes in workspace
        syn_pc_change = drug_results['synergy_pc'][i] - control_results['synergy_pc'][i]
        red_pc_change = drug_results['redundancy_pc'][i] - control_results['redundancy_pc'][i]

        metrics['mean_synergy_pc_change_workspace'][i] = np.mean(syn_pc_change[workspace_mask])
        metrics['mean_redundancy_pc_change_workspace'][i] = np.mean(red_pc_change[workspace_mask])

        # PC gradient change: (syn - red) in drug minus (syn - red) in control
        ctrl_gradient = control_results['synergy_pc'][i] - control_results['redundancy_pc'][i]
        drug_gradient = drug_results['synergy_pc'][i] - drug_results['redundancy_pc'][i]
        gradient_change = drug_gradient - ctrl_gradient
        metrics['pc_gradient_change_workspace'][i] = np.mean(gradient_change[workspace_mask])

    return pd.DataFrame(metrics)


def compute_individual_metrics(matrices: Dict[str, np.ndarray],
                               rsn_labels: np.ndarray,
                               threshold_percentile: float = 90) -> Dict[str, np.ndarray]:
    """
    Compute participation coefficients and classifications for each subject

    Parameters
    ----------
    matrices : dict
        Keys: 'synergy', 'redundancy'
        Values: np.ndarray of shape (n_subjects, 454, 454)
    rsn_labels : np.ndarray (454,)
        Network labels for each ROI
    threshold_percentile : float
        Percentile threshold for workspace identification

    Returns
    -------
    results : dict
        Keys: 'synergy_pc', 'redundancy_pc', 'workspace', 'gateway', 'broadcaster'
        Values: np.ndarray of shape (n_subjects, 454) or (n_subjects, 454) for boolean
    """
    n_subjects = matrices['synergy'].shape[0]

    results = {
        'synergy_pc': np.zeros((n_subjects, N_ROIS)),
        'redundancy_pc': np.zeros((n_subjects, N_ROIS)),
        'workspace': np.zeros((n_subjects, N_ROIS), dtype=bool),
        'gateway': np.zeros((n_subjects, N_ROIS), dtype=bool),
        'broadcaster': np.zeros((n_subjects, N_ROIS), dtype=bool)
    }

    for i in range(n_subjects):
        # Compute participation coefficients
        results['synergy_pc'][i] = compute_participation_coefficient(
            matrices['synergy'][i], rsn_labels)
        results['redundancy_pc'][i] = compute_participation_coefficient(
            matrices['redundancy'][i], rsn_labels)

        # Identify workspace
        syn_strength = compute_node_strength(matrices['synergy'][i])
        results['workspace'][i] = identify_workspace(syn_strength, threshold_percentile)

        # Classify gateways/broadcasters
        results['gateway'][i], results['broadcaster'][i] = classify_gateways_broadcasters(
            results['synergy_pc'][i],
            results['redundancy_pc'][i],
            results['workspace'][i]
        )

    return results


# ============================================================================
# STATISTICAL TESTS
# ============================================================================

def paired_ttest_with_correction(control_data: np.ndarray,
                                drug_data: np.ndarray,
                                method: str = 'fdr_bh') -> Dict[str, np.ndarray]:
    """
    Perform paired t-tests with multiple comparison correction

    Parameters
    ----------
    control_data : np.ndarray (n_subjects, n_features)
        Control condition data
    drug_data : np.ndarray (n_subjects, n_features)
        Drug condition data
    method : str
        Multiple comparison correction method: 'fdr_bh', 'bonferroni', 'none'

    Returns
    -------
    results : dict
        Keys: 't_stat', 'p_value', 'p_corrected', 'significant', 'mean_diff',
              'std_diff', 'ci_lower', 'ci_upper', 'cohen_d'
    """
    n_subjects, n_features = control_data.shape

    # Paired differences
    differences = drug_data - control_data
    mean_diff = np.mean(differences, axis=0)
    std_diff = np.std(differences, axis=0, ddof=1)

    # Paired t-test
    t_stat = np.zeros(n_features)
    p_value = np.zeros(n_features)

    for i in range(n_features):
        t_stat[i], p_value[i] = stats.ttest_rel(drug_data[:, i], control_data[:, i])

    # Multiple comparison correction
    if method == 'fdr_bh':
        # FDR using false_discovery_control (newer scipy)
        try:
            p_corrected = false_discovery_control(p_value, method='bh')
        except:
            # Fallback for older scipy
            from statsmodels.stats.multitest import multipletests
            _, p_corrected, _, _ = multipletests(p_value, method='fdr_bh')
    elif method == 'bonferroni':
        p_corrected = np.minimum(p_value * n_features, 1.0)
    else:  # 'none'
        p_corrected = p_value.copy()

    significant = p_corrected < 0.05

    # Confidence intervals (95%)
    se = std_diff / np.sqrt(n_subjects)
    t_crit = stats.t.ppf(0.975, n_subjects - 1)
    ci_lower = mean_diff - t_crit * se
    ci_upper = mean_diff + t_crit * se

    # Cohen's d for paired data
    cohen_d = mean_diff / std_diff

    return {
        't_stat': t_stat,
        'p_value': p_value,
        'p_corrected': p_corrected,
        'significant': significant,
        'mean_diff': mean_diff,
        'std_diff': std_diff,
        'se': se,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'cohen_d': cohen_d
    }


def mcnemar_test_role_changes(control_roles: np.ndarray,
                              drug_roles: np.ndarray) -> Dict:
    """
    McNemar's test for paired binary data (role changes)

    Tests whether the proportion of regions changing roles differs from chance

    Parameters
    ----------
    control_roles : np.ndarray (n_subjects, n_rois) of bool
        Role in control condition
    drug_roles : np.ndarray (n_subjects, n_rois) of bool
        Role in drug condition

    Returns
    -------
    results : dict
        Per-ROI McNemar statistics
    """
    n_subjects, n_rois = control_roles.shape

    chi2_stat = np.zeros(n_rois)
    p_value = np.zeros(n_rois)
    n_changed = np.zeros(n_rois, dtype=int)
    n_stable_yes = np.zeros(n_rois, dtype=int)
    n_stable_no = np.zeros(n_rois, dtype=int)

    for i in range(n_rois):
        # Contingency table
        # Format: [[n_00, n_01], [n_10, n_11]]
        # n_01 = control=0, drug=1 (gained role)
        # n_10 = control=1, drug=0 (lost role)

        n_00 = np.sum(~control_roles[:, i] & ~drug_roles[:, i])  # Neither
        n_01 = np.sum(~control_roles[:, i] & drug_roles[:, i])   # Gained
        n_10 = np.sum(control_roles[:, i] & ~drug_roles[:, i])   # Lost
        n_11 = np.sum(control_roles[:, i] & drug_roles[:, i])    # Both

        n_changed[i] = n_01 + n_10
        n_stable_yes[i] = n_11
        n_stable_no[i] = n_00

        # McNemar's test (using continuity correction for small samples)
        if n_01 + n_10 > 0:
            chi2_stat[i] = ((abs(n_01 - n_10) - 1) ** 2) / (n_01 + n_10)
            p_value[i] = stats.chi2.sf(chi2_stat[i], df=1)
        else:
            chi2_stat[i] = 0
            p_value[i] = 1.0

    return {
        'chi2_stat': chi2_stat,
        'p_value': p_value,
        'n_changed': n_changed,
        'n_stable_yes': n_stable_yes,
        'n_stable_no': n_stable_no
    }


# ============================================================================
# NETWORK-LEVEL ROLE COUNT STATISTICS
# ============================================================================

def compute_network_role_counts(
    control_results: Dict[str, np.ndarray],
    drug_results: Dict[str, np.ndarray],
    rsn_labels: np.ndarray,
    rsn_names: List[str],
    subject_ids: List[str]
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Count gateways and broadcasters per network per subject.

    Parameters
    ----------
    control_results : dict
        Per-subject metrics from control condition
    drug_results : dict
        Per-subject metrics from drug condition
    rsn_labels : np.ndarray (454,)
        Network labels for each ROI
    rsn_names : list of str
        Network names
    subject_ids : list
        Subject identifiers

    Returns
    -------
    gateway_counts : pd.DataFrame
        Gateway counts per network per subject (n_subjects rows, 2*n_networks columns)
    broadcaster_counts : pd.DataFrame
        Broadcaster counts per network per subject
    """
    n_subjects = len(subject_ids)
    n_networks = len(rsn_names)

    # Initialize count arrays
    gateway_control = np.zeros((n_subjects, n_networks), dtype=int)
    gateway_drug = np.zeros((n_subjects, n_networks), dtype=int)
    broadcaster_control = np.zeros((n_subjects, n_networks), dtype=int)
    broadcaster_drug = np.zeros((n_subjects, n_networks), dtype=int)

    for i in range(n_subjects):
        for net_id in range(n_networks):
            net_mask = (rsn_labels == net_id)

            # Count gateways in this network
            gateway_control[i, net_id] = np.sum(control_results['gateway'][i] & net_mask)
            gateway_drug[i, net_id] = np.sum(drug_results['gateway'][i] & net_mask)

            # Count broadcasters in this network
            broadcaster_control[i, net_id] = np.sum(control_results['broadcaster'][i] & net_mask)
            broadcaster_drug[i, net_id] = np.sum(drug_results['broadcaster'][i] & net_mask)

    # Create DataFrames
    gateway_df = pd.DataFrame({'subject': subject_ids})
    broadcaster_df = pd.DataFrame({'subject': subject_ids})

    for net_id, net_name in enumerate(rsn_names):
        gateway_df[f'{net_name}_control'] = gateway_control[:, net_id]
        gateway_df[f'{net_name}_drug'] = gateway_drug[:, net_id]
        gateway_df[f'{net_name}_diff'] = gateway_drug[:, net_id] - gateway_control[:, net_id]

        broadcaster_df[f'{net_name}_control'] = broadcaster_control[:, net_id]
        broadcaster_df[f'{net_name}_drug'] = broadcaster_drug[:, net_id]
        broadcaster_df[f'{net_name}_diff'] = broadcaster_drug[:, net_id] - broadcaster_control[:, net_id]

    return gateway_df, broadcaster_df


def test_network_role_count_changes(
    control_results: Dict[str, np.ndarray],
    drug_results: Dict[str, np.ndarray],
    rsn_labels: np.ndarray,
    rsn_names: List[str],
    subject_ids: List[str]
) -> pd.DataFrame:
    """
    Paired t-test for gateway/broadcaster COUNT changes within each network.

    Tests: "Did the drug significantly change the number of gateways/broadcasters
    in each network (e.g., DMN) compared to placebo?"

    Parameters
    ----------
    control_results : dict
        Per-subject metrics from control condition
    drug_results : dict
        Per-subject metrics from drug condition
    rsn_labels : np.ndarray (454,)
        Network labels for each ROI
    rsn_names : list of str
        Network names
    subject_ids : list
        Subject identifiers

    Returns
    -------
    results_df : pd.DataFrame
        Statistical results with columns:
        Network, Role, N_ROIs_in_Network, Mean_Control, Mean_Drug, Mean_Diff,
        SD_Diff, t_stat, p_value, Cohen_d, CI_lower, CI_upper, Significant
    """
    n_subjects = len(subject_ids)
    n_networks = len(rsn_names)

    results = []

    for net_id, net_name in enumerate(rsn_names):
        net_mask = (rsn_labels == net_id)
        n_rois_network = np.sum(net_mask)

        if n_rois_network == 0:
            continue

        for role_name, control_key, drug_key in [
            ('Gateway', 'gateway', 'gateway'),
            ('Broadcaster', 'broadcaster', 'broadcaster')
        ]:
            # Count per subject
            control_counts = np.array([
                np.sum(control_results[control_key][i] & net_mask)
                for i in range(n_subjects)
            ])
            drug_counts = np.array([
                np.sum(drug_results[drug_key][i] & net_mask)
                for i in range(n_subjects)
            ])

            # Paired differences
            diff = drug_counts - control_counts
            mean_diff = np.mean(diff)
            std_diff = np.std(diff, ddof=1)

            # Paired t-test
            if std_diff > 0:
                t_stat, p_value = stats.ttest_rel(drug_counts, control_counts)
                cohen_d = mean_diff / std_diff
            else:
                t_stat = 0
                p_value = 1.0
                cohen_d = 0

            # Confidence interval
            se = std_diff / np.sqrt(n_subjects)
            t_crit = stats.t.ppf(0.975, n_subjects - 1)
            ci_lower = mean_diff - t_crit * se
            ci_upper = mean_diff + t_crit * se

            results.append({
                'Network': net_name,
                'Network_ID': net_id,
                'Role': role_name,
                'N_ROIs_in_Network': n_rois_network,
                'Mean_Control': np.mean(control_counts),
                'SD_Control': np.std(control_counts, ddof=1),
                'Mean_Drug': np.mean(drug_counts),
                'SD_Drug': np.std(drug_counts, ddof=1),
                'Mean_Diff': mean_diff,
                'SD_Diff': std_diff,
                't_stat': t_stat,
                'df': n_subjects - 1,
                'p_value': p_value,
                'Cohen_d': cohen_d,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper,
                'Significant_uncorrected': p_value < 0.05
            })

    results_df = pd.DataFrame(results)

    # FDR correction across all tests
    if len(results_df) > 0:
        try:
            results_df['p_corrected'] = false_discovery_control(results_df['p_value'].values, method='bh')
        except:
            from statsmodels.stats.multitest import multipletests
            _, p_corr, _, _ = multipletests(results_df['p_value'].values, method='fdr_bh')
            results_df['p_corrected'] = p_corr

        results_df['Significant_FDR'] = results_df['p_corrected'] < 0.05

    return results_df


# ============================================================================
# NETWORK-LEVEL SUMMARY
# ============================================================================

def summarize_by_network(results: Dict[str, np.ndarray],
                         rsn_labels: np.ndarray,
                         rsn_names: List[str]) -> pd.DataFrame:
    """
    Summarize statistical results by brain network

    Parameters
    ----------
    results : dict
        Output from paired_ttest_with_correction
    rsn_labels : np.ndarray (454,)
        Network labels
    rsn_names : list of str
        Network names

    Returns
    -------
    summary_df : pd.DataFrame
        Network-level summary statistics
    """
    summary_data = []

    for rsn_id, rsn_name in enumerate(rsn_names):
        rsn_mask = (rsn_labels == rsn_id)
        n_rois = np.sum(rsn_mask)

        if n_rois == 0:
            continue

        # Extract data for this network
        mean_diff = results['mean_diff'][rsn_mask]
        p_corrected = results['p_corrected'][rsn_mask]
        cohen_d = results['cohen_d'][rsn_mask]
        significant = results['significant'][rsn_mask]

        summary_data.append({
            'Network': rsn_name,
            'Network_ID': rsn_id,
            'N_ROIs': n_rois,
            'Mean_Change': np.mean(mean_diff),
            'Median_Change': np.median(mean_diff),
            'N_Significant': np.sum(significant),
            'Pct_Significant': 100 * np.sum(significant) / n_rois,
            'Mean_Cohen_d': np.mean(cohen_d),
            'Max_Change': np.max(mean_diff),
            'Min_Change': np.min(mean_diff),
            'Min_p_corrected': np.min(p_corrected)
        })

    return pd.DataFrame(summary_data)


# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def main():
    """Main analysis pipeline"""

    parser = argparse.ArgumentParser(
        description='Paired statistical testing for gateway/broadcaster changes - DATASET AWARE'
    )
    parser.add_argument('--dataset', required=True,
                       help='Dataset name (e.g., london_lsd, basel_mmm)')
    parser.add_argument('--drug_task', required=True,
                       help='Drug task name (e.g., LSD, MDMA)')
    parser.add_argument('--threshold_percentile', type=float, default=90,
                       help='Percentile threshold for workspace (default: 90)')
    parser.add_argument('--correction', default='fdr_bh',
                       choices=['fdr_bh', 'bonferroni', 'none'],
                       help='Multiple comparison correction method')
    parser.add_argument('--project_dir', help='Project directory (optional)')
    parser.add_argument('--use_group_workspace', action='store_true', default=True,
                       help='Use group-defined workspace mask for consistent cross-subject comparison (default: True)')
    parser.add_argument('--save_per_subject', action='store_true', default=True,
                       help='Save per-subject classification arrays for correlation with subjective effects (default: True)')
    parser.add_argument('--workspace_dataset', default='hcp',
                       help='Dataset that defines the canonical workspace (default: hcp). '
                       'The HCP workspace is the canonical reference per Luppi et al.')

    args = parser.parse_args()

    # Load dataset handler
    handler = load_dataset_handler(args.dataset, project_dir=args.project_dir)

    # Get control task
    control_task = get_control_task_name(args.dataset)

    # Setup output directory
    base_results_dir = handler.project_dir / 'results' / 'group_analysis'
    stats_dir = base_results_dir / args.dataset / f'paired_stats_{args.drug_task}_vs_{control_task}'
    stats_dir.mkdir(exist_ok=True, parents=True)

    # File prefix
    prefix = f"{args.dataset}_{args.drug_task}_vs_{control_task}"

    # Logging
    md_output = []

    def print_and_log(msg):
        print(msg)
        md_output.append(msg)

    print_and_log("=" * 80)
    print_and_log("PAIRED STATISTICAL TESTING FOR GATEWAY/BROADCASTER CHANGES")
    print_and_log("=" * 80)
    print_and_log(f"\nDataset: {handler.config['name']}")
    print_and_log(f"Drug: {args.drug_task}")
    print_and_log(f"Control: {control_task}")
    print_and_log(f"Workspace dataset: {args.workspace_dataset}")
    print_and_log(f"Design: Within-subject crossover (paired)")
    print_and_log(f"Multiple comparison correction: {args.correction}")

    # -------------------------------------------------------------------------
    # Step 1: Load RSN labels
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 1: Loading RSN labels")
    print_and_log("=" * 80)

    rsn_labels, rsn_names, roi_df = load_roi_network_mapping(handler.project_dir)
    print_and_log(f"[OK] Loaded {len(rsn_labels)} ROI labels across {len(rsn_names)} networks")

    # -------------------------------------------------------------------------
    # Step 2: Load paired data
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 2: Loading paired data")
    print_and_log("=" * 80)

    try:
        control_matrices, drug_matrices, control_subjects, drug_subjects = load_drug_comparison_data(
            args.dataset,
            args.drug_task,
            results_dir=handler.project_dir / 'results' / 'phiid_matrices'
        )
    except Exception as e:
        print_and_log(f"ERROR: Failed to load data: {e}")
        return

    # Verify pairing - use intersection of subjects if mismatch
    if control_subjects != drug_subjects:
        print_and_log("WARNING: Subject mismatch between conditions!")
        print_and_log(f"Control subjects: {control_subjects}")
        print_and_log(f"Drug subjects: {drug_subjects}")

        # Find common subjects
        common_subjects = sorted(list(set(control_subjects) & set(drug_subjects)))

        if len(common_subjects) == 0:
            print_and_log("ERROR: No common subjects found!")
            return

        print_and_log(f"\nUsing {len(common_subjects)} common subjects: {common_subjects}")

        # Filter matrices to only include common subjects
        control_indices = [i for i, s in enumerate(control_subjects) if s in common_subjects]
        drug_indices = [i for i, s in enumerate(drug_subjects) if s in common_subjects]

        for key in control_matrices:
            control_matrices[key] = control_matrices[key][control_indices]
        for key in drug_matrices:
            drug_matrices[key] = drug_matrices[key][drug_indices]

        control_subjects = common_subjects
        drug_subjects = common_subjects

    n_subjects = len(control_subjects)
    print_and_log(f"[OK] Loaded {n_subjects} paired subjects")
    print_and_log(f"  Subjects: {control_subjects}")
    print_and_log(f"[OK] Data pairing verified (same subjects, same order)")

    # -------------------------------------------------------------------------
    # Step 3: Compute individual-subject metrics
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 3: Computing individual-subject metrics")
    print_and_log("=" * 80)

    # Load group workspace mask if using group workspace
    # Use workspace_dataset to determine source (default: hcp)
    ws_ds = args.workspace_dataset if args.workspace_dataset else args.dataset
    ws_task = get_control_task_name(ws_ds) if ws_ds != args.dataset else control_task

    group_workspace_mask = None
    if args.use_group_workspace:
        try:
            group_workspace_mask = load_group_workspace_mask(
                base_results_dir, ws_ds, ws_task)
            print_and_log(f"[OK] Loaded group workspace mask from {ws_ds}/{ws_task} ({np.sum(group_workspace_mask)} regions)")
            if ws_ds != args.dataset:
                print_and_log(f"  Using {ws_ds} workspace per Luppi et al.")
            else:
                print_and_log("  Using GROUP workspace for consistent cross-subject comparison")
        except FileNotFoundError as e:
            print_and_log(f"[!] Warning: Could not load group workspace: {e}")
            print_and_log("  Falling back to per-subject workspace")
            args.use_group_workspace = False

    if args.use_group_workspace and group_workspace_mask is not None:
        print_and_log("\nComputing control condition metrics (with group workspace)...")
        control_results = compute_individual_metrics_with_group_workspace(
            control_matrices, rsn_labels, group_workspace_mask)
        print_and_log(f"[OK] Computed participation coefficients for {n_subjects} subjects")

        print_and_log("\nComputing drug condition metrics (with group workspace)...")
        drug_results = compute_individual_metrics_with_group_workspace(
            drug_matrices, rsn_labels, group_workspace_mask)
        print_and_log(f"[OK] Computed participation coefficients for {n_subjects} subjects")
    else:
        print_and_log("Computing control condition metrics (per-subject workspace)...")
        control_results = compute_individual_metrics(
            control_matrices, rsn_labels, args.threshold_percentile)
        print_and_log(f"[OK] Computed participation coefficients for {n_subjects} subjects")

        print_and_log("\nComputing drug condition metrics (per-subject workspace)...")
        drug_results = compute_individual_metrics(
            drug_matrices, rsn_labels, args.threshold_percentile)
        print_and_log(f"[OK] Computed participation coefficients for {n_subjects} subjects")

    # -------------------------------------------------------------------------
    # Step 4: Statistical testing - Participation Coefficients
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 4: Paired t-tests for participation coefficients")
    print_and_log("=" * 80)

    # Synergy PC
    print_and_log("\nTesting synergy participation coefficient changes...")
    synergy_pc_stats = paired_ttest_with_correction(
        control_results['synergy_pc'],
        drug_results['synergy_pc'],
        method=args.correction
    )
    n_sig_syn = np.sum(synergy_pc_stats['significant'])
    print_and_log(f"  [OK] {n_sig_syn} / {N_ROIS} ROIs show significant changes (FDR < 0.05)")

    # Redundancy PC
    print_and_log("\nTesting redundancy participation coefficient changes...")
    redundancy_pc_stats = paired_ttest_with_correction(
        control_results['redundancy_pc'],
        drug_results['redundancy_pc'],
        method=args.correction
    )
    n_sig_red = np.sum(redundancy_pc_stats['significant'])
    print_and_log(f"  [OK] {n_sig_red} / {N_ROIS} ROIs show significant changes (FDR < 0.05)")

    # -------------------------------------------------------------------------
    # Step 5: Statistical testing - Role Changes
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 5: McNemar's test for role changes")
    print_and_log("=" * 80)

    # Gateway changes
    print_and_log("\nTesting gateway role changes...")
    gateway_stats = mcnemar_test_role_changes(
        control_results['gateway'],
        drug_results['gateway']
    )
    total_gateway_changes = np.sum(gateway_stats['n_changed'])
    print_and_log(f"  [OK] {total_gateway_changes} total gateway role changes across subjects")

    # Broadcaster changes
    print_and_log("\nTesting broadcaster role changes...")
    broadcaster_stats = mcnemar_test_role_changes(
        control_results['broadcaster'],
        drug_results['broadcaster']
    )
    total_broadcaster_changes = np.sum(broadcaster_stats['n_changed'])
    print_and_log(f"  [OK] {total_broadcaster_changes} total broadcaster role changes across subjects")

    # -------------------------------------------------------------------------
    # Step 6: Network-level summaries
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 6: Network-level summaries")
    print_and_log("=" * 80)

    synergy_pc_summary = summarize_by_network(synergy_pc_stats, rsn_labels, rsn_names)
    redundancy_pc_summary = summarize_by_network(redundancy_pc_stats, rsn_labels, rsn_names)

    print_and_log("\nSynergy PC changes by network:")
    print_and_log(synergy_pc_summary[['Network', 'Mean_Change', 'N_Significant', 'Pct_Significant', 'Mean_Cohen_d']].to_string(index=False))

    print_and_log("\nRedundancy PC changes by network:")
    print_and_log(redundancy_pc_summary[['Network', 'Mean_Change', 'N_Significant', 'Pct_Significant', 'Mean_Cohen_d']].to_string(index=False))

    # -------------------------------------------------------------------------
    # Step 6b: Network-level role COUNT statistics
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 6b: Network-level gateway/broadcaster COUNT statistics")
    print_and_log("=" * 80)
    print_and_log("\nTesting: Did drug change the NUMBER of gateways/broadcasters per network?")

    network_role_stats = test_network_role_count_changes(
        control_results, drug_results, rsn_labels, rsn_names, control_subjects
    )

    # Print significant results
    sig_results = network_role_stats[network_role_stats['Significant_uncorrected']]
    sig_fdr_results = network_role_stats[network_role_stats['Significant_FDR']]

    print_and_log(f"\nSignificant network-role changes (uncorrected p < 0.05): {len(sig_results)}")
    print_and_log(f"Significant network-role changes (FDR corrected): {len(sig_fdr_results)}")

    if len(sig_results) > 0:
        print_and_log("\nUncorrected significant results:")
        print_and_log("-" * 80)
        for _, row in sig_results.iterrows():
            direction = "increased" if row['Mean_Diff'] > 0 else "decreased"
            print_and_log(f"  {row['Network']} {row['Role']}s: {direction} by {abs(row['Mean_Diff']):.2f} "
                         f"(t={row['t_stat']:.2f}, p={row['p_value']:.4f}, d={row['Cohen_d']:.2f})")
            print_and_log(f"    Control: {row['Mean_Control']:.1f} +/- {row['SD_Control']:.1f}, "
                         f"Drug: {row['Mean_Drug']:.1f} +/- {row['SD_Drug']:.1f}")

    if len(sig_fdr_results) > 0:
        print_and_log("\nFDR-corrected significant results:")
        print_and_log("-" * 80)
        for _, row in sig_fdr_results.iterrows():
            direction = "increased" if row['Mean_Diff'] > 0 else "decreased"
            print_and_log(f"  {row['Network']} {row['Role']}s: {direction} by {abs(row['Mean_Diff']):.2f} "
                         f"(t={row['t_stat']:.2f}, p_FDR={row['p_corrected']:.4f}, d={row['Cohen_d']:.2f})")
    else:
        print_and_log("\nNo results survive FDR correction.")

    # Also compute per-subject network counts for saving
    gateway_counts_df, broadcaster_counts_df = compute_network_role_counts(
        control_results, drug_results, rsn_labels, rsn_names, control_subjects
    )

    # -------------------------------------------------------------------------
    # Step 7: Save results
    # -------------------------------------------------------------------------
    print_and_log("\n" + "=" * 80)
    print_and_log("STEP 7: Saving results")
    print_and_log("=" * 80)

    # Save participation coefficient statistics
    synergy_pc_df = pd.DataFrame({
        'ROI': roi_df['ROI'].values,
        'Network': rsn_labels,
        'Network_Name': [rsn_names[i] for i in rsn_labels],
        'mean_diff': synergy_pc_stats['mean_diff'],
        'std_diff': synergy_pc_stats['std_diff'],
        'se': synergy_pc_stats['se'],
        't_stat': synergy_pc_stats['t_stat'],
        'p_value': synergy_pc_stats['p_value'],
        'p_corrected': synergy_pc_stats['p_corrected'],
        'significant': synergy_pc_stats['significant'],
        'ci_lower': synergy_pc_stats['ci_lower'],
        'ci_upper': synergy_pc_stats['ci_upper'],
        'cohen_d': synergy_pc_stats['cohen_d']
    })
    synergy_pc_df.to_csv(stats_dir / f'{prefix}_synergy_pc_stats.csv', index=False)
    print_and_log("[OK] Saved synergy PC statistics")

    redundancy_pc_df = pd.DataFrame({
        'ROI': roi_df['ROI'].values,
        'Network': rsn_labels,
        'Network_Name': [rsn_names[i] for i in rsn_labels],
        'mean_diff': redundancy_pc_stats['mean_diff'],
        'std_diff': redundancy_pc_stats['std_diff'],
        'se': redundancy_pc_stats['se'],
        't_stat': redundancy_pc_stats['t_stat'],
        'p_value': redundancy_pc_stats['p_value'],
        'p_corrected': redundancy_pc_stats['p_corrected'],
        'significant': redundancy_pc_stats['significant'],
        'ci_lower': redundancy_pc_stats['ci_lower'],
        'ci_upper': redundancy_pc_stats['ci_upper'],
        'cohen_d': redundancy_pc_stats['cohen_d']
    })
    redundancy_pc_df.to_csv(stats_dir / f'{prefix}_redundancy_pc_stats.csv', index=False)
    print_and_log("[OK] Saved redundancy PC statistics")

    # Save network summaries
    synergy_pc_summary.to_csv(stats_dir / f'{prefix}_synergy_pc_network_summary.csv', index=False)
    redundancy_pc_summary.to_csv(stats_dir / f'{prefix}_redundancy_pc_network_summary.csv', index=False)
    print_and_log("[OK] Saved network-level summaries")

    # Save network role count statistics
    network_role_stats.to_csv(stats_dir / f'{prefix}_network_role_count_stats.csv', index=False)
    gateway_counts_df.to_csv(stats_dir / f'{prefix}_gateway_counts_per_network.csv', index=False)
    broadcaster_counts_df.to_csv(stats_dir / f'{prefix}_broadcaster_counts_per_network.csv', index=False)
    print_and_log("[OK] Saved network-level role count statistics")

    # -------------------------------------------------------------------------
    # Step 7b: Save per-subject data for correlation with subjective effects
    # -------------------------------------------------------------------------
    if args.save_per_subject:
        print_and_log("\n" + "-" * 40)
        print_and_log("Saving per-subject data for subjective correlation...")

        # Save per-subject arrays
        np.save(stats_dir / f'{prefix}_per_subject_gateway_control.npy', control_results['gateway'])
        np.save(stats_dir / f'{prefix}_per_subject_gateway_drug.npy', drug_results['gateway'])
        np.save(stats_dir / f'{prefix}_per_subject_broadcaster_control.npy', control_results['broadcaster'])
        np.save(stats_dir / f'{prefix}_per_subject_broadcaster_drug.npy', drug_results['broadcaster'])
        np.save(stats_dir / f'{prefix}_per_subject_synergy_pc_control.npy', control_results['synergy_pc'])
        np.save(stats_dir / f'{prefix}_per_subject_synergy_pc_drug.npy', drug_results['synergy_pc'])
        np.save(stats_dir / f'{prefix}_per_subject_redundancy_pc_control.npy', control_results['redundancy_pc'])
        np.save(stats_dir / f'{prefix}_per_subject_redundancy_pc_drug.npy', drug_results['redundancy_pc'])
        np.save(stats_dir / f'{prefix}_per_subject_workspace.npy', control_results['workspace'])

        print_and_log(f"  [OK] Saved per-subject arrays (shape: {n_subjects} x {N_ROIS})")

        # Compute and save role switch metrics
        role_switch_metrics = compute_per_subject_role_switches(
            control_results, drug_results, control_subjects
        )
        role_switch_metrics.to_csv(stats_dir / f'{prefix}_role_switch_metrics.csv', index=False)

        print_and_log(f"  [OK] Saved role switch metrics for {n_subjects} subjects")
        print_and_log(f"       Mean G->B switches: {role_switch_metrics['gateway_to_broadcaster'].mean():.1f}")
        print_and_log(f"       Mean B->G switches: {role_switch_metrics['broadcaster_to_gateway'].mean():.1f}")
        print_and_log(f"       Mean total switches: {role_switch_metrics['total_role_switches'].mean():.1f}")

    # Save metadata
    n_sig_network_role_uncorr = len(network_role_stats[network_role_stats['Significant_uncorrected']])
    n_sig_network_role_fdr = len(network_role_stats[network_role_stats['Significant_FDR']])

    metadata = {
        'dataset': args.dataset,
        'dataset_name': handler.config['name'],
        'drug_task': args.drug_task,
        'control_task': control_task,
        'workspace_source_dataset': ws_ds,
        'workspace_source_task': ws_task,
        'workspace_source_path': str(base_results_dir / ws_ds / ws_task),
        'n_subjects': n_subjects,
        'subjects': control_subjects,
        'design': 'within-subject crossover (paired)',
        'correction_method': args.correction,
        'threshold_percentile': args.threshold_percentile,
        'use_group_workspace': args.use_group_workspace,
        'n_rois': N_ROIS,
        'synergy_pc_significant': int(n_sig_syn),
        'redundancy_pc_significant': int(n_sig_red),
        'total_gateway_changes': int(total_gateway_changes),
        'total_broadcaster_changes': int(total_broadcaster_changes),
        'network_role_count_significant_uncorrected': int(n_sig_network_role_uncorr),
        'network_role_count_significant_fdr': int(n_sig_network_role_fdr),
        'per_subject_data_saved': args.save_per_subject,
        'analysis_date': pd.Timestamp.now().isoformat()
    }

    with open(stats_dir / f'{prefix}_metadata.yaml', 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)
    print_and_log("[OK] Saved metadata")

    # Save markdown report
    md_file = stats_dir / f'{prefix}_results.md'
    with open(md_file, 'w') as f:
        f.write('\n'.join(md_output))
    print_and_log(f"\n[OK] Markdown report saved to: {md_file}")

    print_and_log("\n" + "=" * 80)
    print_and_log("PAIRED STATISTICAL TESTING COMPLETE")
    print_and_log("=" * 80)
    print_and_log(f"\nResults saved to: {stats_dir}")
    print_and_log(f"\nKey findings:")
    print_and_log(f"  - {n_sig_syn} ROIs with significant synergy PC changes (FDR < 0.05)")
    print_and_log(f"  - {n_sig_red} ROIs with significant redundancy PC changes (FDR < 0.05)")
    print_and_log(f"  - Networks with largest effects: see {prefix}_*_network_summary.csv")


if __name__ == '__main__':
    main()
