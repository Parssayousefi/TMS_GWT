#!/usr/bin/env python3
"""
Condition-Specific Gateway/Broadcaster Analysis
================================================

Compare gateway and broadcaster classifications between drug and control conditions:
1. Load pre-computed control results from script 06 (no redundant computation)
2. Compute workspace/gateway/broadcaster for drug condition
3. Identify regions that change roles under drug conditions
4. Analyze network-specific effects

**UPDATED**: Now dataset-aware - works with any configured dataset
**REFACTORED**: Loads control results from script 06 output instead of recomputing

Requires: Run 06_identify_gateways_broadcasters.py first for the control condition
"""

import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from load_roi_mapping import load_roi_network_mapping
from analysis_loader import (
    load_group_all_measures,
    get_control_task_name,
    get_drug_task_names
)
from dataset_handler import load_dataset_handler

# Parameters
N_ROIS = 454
N_CORTICAL = 400  # Schaefer
N_SUBCORTICAL = 54  # Tian

def get_rsn_labels(base_dir):
    """Get RSN labels for all ROIs"""
    rsn_labels, rsn_names, roi_df = load_roi_network_mapping(base_dir)
    return rsn_labels, rsn_names

def compute_node_strength(matrix):
    """Compute node strength (sum of connections for each node)"""
    return np.sum(matrix, axis=1)

def identify_workspace(synergy_strength, redundancy_strength):
    """
    Identify synergistic workspace regions using gradient method

    Following Luppi et al. (2024): "regions where synergy predominates"
    means synergy rank > redundancy rank (not top percentile)

    Parameters
    ----------
    synergy_strength : array (454,)
        Node strength for synergy
    redundancy_strength : array (454,)
        Node strength for redundancy

    Returns
    -------
    workspace_mask : array (454,) of bool
        True for workspace regions (where synergy rank > redundancy rank)
    """
    # Rank regions by strength (higher strength = higher rank)
    synergy_ranks = np.argsort(np.argsort(synergy_strength)) + 1
    redundancy_ranks = np.argsort(np.argsort(redundancy_strength)) + 1

    # Gradient: synergy rank - redundancy rank
    # Positive = more synergistic than redundant (relative to other regions)
    gradient = synergy_ranks - redundancy_ranks

    # Workspace = regions where synergy predominates (gradient > 0)
    workspace_mask = gradient > 0

    return workspace_mask

def compute_participation_coefficient(connectivity_matrix, module_labels):
    """
    Compute participation coefficient for each node

    From Rubinov & Sporns (2010) and Luppi et al. (2024):
    P_i = 1 - sum_s=1^M (κ_is / k_i)^2
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

def classify_gateways_broadcasters(synergy_pc, redundancy_pc, workspace_mask):
    """
    Classify workspace regions into gateways vs broadcasters

    Returns
    -------
    gateway_mask : array (454,) of bool
    broadcaster_mask : array (454,) of bool
    pc_gradient : array (454,)
        Synergy PC rank - Redundancy PC rank (for workspace regions only)
    """
    workspace_indices = np.where(workspace_mask)[0]

    # Rank participation coefficients (within workspace only)
    synergy_pc_workspace = synergy_pc[workspace_mask]
    redundancy_pc_workspace = redundancy_pc[workspace_mask]

    synergy_pc_ranks = np.argsort(np.argsort(synergy_pc_workspace)) + 1
    redundancy_pc_ranks = np.argsort(np.argsort(redundancy_pc_workspace)) + 1

    # Gradient: synergy PC rank - redundancy PC rank
    pc_gradient_workspace = synergy_pc_ranks - redundancy_pc_ranks

    # Full gradient array (NaN for non-workspace)
    pc_gradient = np.full(N_ROIS, np.nan)
    pc_gradient[workspace_mask] = pc_gradient_workspace

    # Classify
    gateway_mask = np.zeros(N_ROIS, dtype=bool)
    broadcaster_mask = np.zeros(N_ROIS, dtype=bool)

    gateway_mask[workspace_indices[pc_gradient_workspace > 0]] = True
    broadcaster_mask[workspace_indices[pc_gradient_workspace < 0]] = True

    return gateway_mask, broadcaster_mask, pc_gradient

def compare_conditions(control_results, drug_results, rsn_labels, rsn_names):
    """
    Compare gateway/broadcaster classifications between conditions

    Returns comparison statistics and changed regions
    """
    # Identify regions that change roles
    control_gw = control_results['gateway_mask']
    control_bc = control_results['broadcaster_mask']
    drug_gw = drug_results['gateway_mask']
    drug_bc = drug_results['broadcaster_mask']

    # Changed roles
    became_gateway = (~control_gw) & drug_gw
    lost_gateway = control_gw & (~drug_gw)
    became_broadcaster = (~control_bc) & drug_bc
    lost_broadcaster = control_bc & (~drug_bc)

    # Workspace membership changes
    control_ws = control_results['workspace_mask']
    drug_ws = drug_results['workspace_mask']
    entered_workspace = (~control_ws) & drug_ws
    left_workspace = control_ws & (~drug_ws)

    # Changes in participation coefficients
    delta_synergy_pc = drug_results['synergy_pc'] - control_results['synergy_pc']
    delta_redundancy_pc = drug_results['redundancy_pc'] - control_results['redundancy_pc']

    comparison = {
        'became_gateway': became_gateway,
        'lost_gateway': lost_gateway,
        'became_broadcaster': became_broadcaster,
        'lost_broadcaster': lost_broadcaster,
        'entered_workspace': entered_workspace,
        'left_workspace': left_workspace,
        'delta_synergy_pc': delta_synergy_pc,
        'delta_redundancy_pc': delta_redundancy_pc,
        'stable_gateway': control_gw & drug_gw,
        'stable_broadcaster': control_bc & drug_bc
    }

    # Summary statistics by RSN
    summary_by_rsn = []
    for rsn_id, rsn_name in enumerate(rsn_names):
        rsn_mask = (rsn_labels == rsn_id)
        summary_by_rsn.append({
            'RSN': rsn_name,
            'control_gateways': np.sum(control_gw & rsn_mask),
            'drug_gateways': np.sum(drug_gw & rsn_mask),
            'control_broadcasters': np.sum(control_bc & rsn_mask),
            'drug_broadcasters': np.sum(drug_bc & rsn_mask),
            'became_gateway': np.sum(became_gateway & rsn_mask),
            'lost_gateway': np.sum(lost_gateway & rsn_mask),
            'became_broadcaster': np.sum(became_broadcaster & rsn_mask),
            'lost_broadcaster': np.sum(lost_broadcaster & rsn_mask),
            'entered_workspace': np.sum(entered_workspace & rsn_mask),
            'left_workspace': np.sum(left_workspace & rsn_mask),
            'mean_delta_syn_pc': np.mean(delta_synergy_pc[rsn_mask]),
            'mean_delta_red_pc': np.mean(delta_redundancy_pc[rsn_mask])
        })

    return comparison, pd.DataFrame(summary_by_rsn)

def main():
    """Main analysis pipeline"""
    parser = argparse.ArgumentParser(
        description='Condition-specific gateway/broadcaster comparison - DATASET AWARE'
    )
    parser.add_argument('--dataset', required=True,
                       help='Dataset name (e.g., london_lsd, basel_mmm)')
    parser.add_argument('--drug_task', required=True,
                       help='Drug task name (e.g., LSD, MDMA)')
    parser.add_argument('--project_dir', help='Project directory (optional)')
    parser.add_argument('--workspace_dataset', default='hcp',
                       help='Dataset that defines the canonical workspace (default: hcp). '
                       'The HCP workspace is the canonical reference per Luppi et al.')

    args = parser.parse_args()

    # Load dataset handler
    handler = load_dataset_handler(args.dataset, project_dir=args.project_dir)

    # Get control task
    control_task = get_control_task_name(args.dataset)

    # Setup output directory: results/group_analysis/{dataset}/comparison_{drug}_vs_{control}/
    base_results_dir = handler.project_dir / 'results' / 'group_analysis'
    comparison_dir = base_results_dir / args.dataset / f'comparison_{args.drug_task}_vs_{control_task}'
    comparison_dir.mkdir(exist_ok=True, parents=True)

    md_output = []

    def print_and_log(msg):
        """Print to console and log to markdown"""
        print(msg)
        md_output.append(msg)

    print_and_log("="*70)
    print_and_log("CONDITION-SPECIFIC GATEWAY/BROADCASTER ANALYSIS")
    print_and_log("="*70)
    print_and_log(f"\nDataset: {handler.config['name']}")
    print_and_log(f"Comparing {args.drug_task} vs {control_task}")

    # -------------------------------------------------------------------------
    # Step 1: Get RSN labels
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 1: Loading RSN labels")
    print_and_log("="*70)

    rsn_labels, rsn_names = get_rsn_labels(handler.project_dir)
    print_and_log("[OK] Loaded RSN assignments")

    # -------------------------------------------------------------------------
    # Step 2: Load pre-computed control results from script 06
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 2: Loading pre-computed control results from script 06")
    print_and_log("="*70)

    # Determine workspace source: use workspace_dataset if different from analysis dataset
    ws_dataset = args.workspace_dataset if args.workspace_dataset else args.dataset
    ws_task = get_control_task_name(ws_dataset) if ws_dataset != args.dataset else control_task
    ws_dir = base_results_dir / ws_dataset / ws_task
    ws_prefix = f"{ws_dataset}_{ws_task}"

    # Control results (gateway/broadcaster/PC) come from the analysis dataset's control condition
    control_dir = base_results_dir / args.dataset / control_task
    control_prefix = f"{args.dataset}_{control_task}"

    # Check workspace mask from workspace dataset
    ws_mask_file = ws_dir / f'{ws_prefix}_workspace_mask.npy'
    if not ws_mask_file.exists():
        print_and_log(f"ERROR: Workspace mask not found: {ws_mask_file}")
        print_and_log(f"\nPlease run script 05 for {ws_dataset} first:")
        print_and_log(f"  python scripts/analysis/05_group_average_workspace.py --dataset {ws_dataset}")
        return

    # Check that script 06 has been run for the analysis dataset's control condition
    required_files = [
        f'{control_prefix}_gateway_mask.npy',
        f'{control_prefix}_broadcaster_mask.npy',
        f'{control_prefix}_synergy_participation_coef.npy',
        f'{control_prefix}_redundancy_participation_coef.npy',
        f'{control_prefix}_group_avg_synergy.npy',
        f'{control_prefix}_group_avg_redundancy.npy'
    ]

    missing_files = [f for f in required_files if not (control_dir / f).exists()]
    if missing_files:
        print_and_log(f"ERROR: Missing control condition files from script 06:")
        for f in missing_files:
            print_and_log(f"  - {f}")
        print_and_log(f"\nPlease run script 06 first:")
        print_and_log(f"  python scripts/analysis/06_identify_gateways_broadcasters.py --dataset {args.dataset}")
        return

    # Load workspace mask from workspace dataset (e.g., HCP)
    control_workspace = np.load(ws_mask_file)

    # Load gateway/broadcaster/PC results from analysis dataset's control condition
    try:
        control_gateway_mask = np.load(control_dir / f'{control_prefix}_gateway_mask.npy')
        control_broadcaster_mask = np.load(control_dir / f'{control_prefix}_broadcaster_mask.npy')
        control_syn_pc = np.load(control_dir / f'{control_prefix}_synergy_participation_coef.npy')
        control_red_pc = np.load(control_dir / f'{control_prefix}_redundancy_participation_coef.npy')
        control_synergy = np.load(control_dir / f'{control_prefix}_group_avg_synergy.npy')
        control_redundancy = np.load(control_dir / f'{control_prefix}_group_avg_redundancy.npy')
        control_pc_gradient = np.load(control_dir / f'{control_prefix}_pc_gradient.npy')
    except Exception as e:
        print_and_log(f"ERROR: Failed to load control results: {e}")
        return

    if ws_dataset != args.dataset:
        print_and_log(f"[OK] Loaded workspace mask from {ws_dataset}/{ws_task} ({control_workspace.sum()} regions)")
        print_and_log(f"[OK] Loaded control gateway/broadcaster from {args.dataset}/{control_task}")
    else:
        print_and_log(f"[OK] Loaded pre-computed {control_task} results from script 06")
    print_and_log(f"  Workspace: {control_workspace.sum()} regions")
    print_and_log(f"  Gateways: {control_gateway_mask.sum()}")
    print_and_log(f"  Broadcasters: {control_broadcaster_mask.sum()}")

    control_results = {
        'synergy': control_synergy,
        'redundancy': control_redundancy,
        'synergy_pc': control_syn_pc,
        'redundancy_pc': control_red_pc,
        'workspace_mask': control_workspace,
        'gateway_mask': control_gateway_mask,
        'broadcaster_mask': control_broadcaster_mask,
        'pc_gradient': control_pc_gradient
    }

    # -------------------------------------------------------------------------
    # Step 3: Load drug condition data and compute from scratch
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log(f"STEP 3: Loading {args.drug_task} data")
    print_and_log("="*70)

    try:
        drug_matrices, drug_subjects = load_group_all_measures(
            args.dataset,
            args.drug_task,
            results_dir=handler.project_dir / 'results' / 'phiid_matrices'
        )
    except Exception as e:
        print_and_log(f"ERROR: Failed to load drug data: {e}")
        print_and_log("\nTroubleshooting:")
        print_and_log("  1. Check that ΦID matrices have been computed for the drug condition")
        print_and_log("  2. Verify task name in dataset configuration")
        print_and_log(f"  3. Run: python scripts/utils/analysis_loader.py --dataset {args.dataset}")
        return

    print_and_log(f"[OK] Loaded {len(drug_subjects)} {args.drug_task} subjects")

    # -------------------------------------------------------------------------
    # Step 4: Compute workspace/gateway/broadcaster for drug condition
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log(f"STEP 4: Computing {args.drug_task} workspace and classification")
    print_and_log("="*70)

    drug_synergy = np.mean(drug_matrices['synergy'], axis=0)
    drug_redundancy = np.mean(drug_matrices['redundancy'], axis=0)

    drug_syn_strength = compute_node_strength(drug_synergy)
    drug_red_strength = compute_node_strength(drug_redundancy)
    drug_workspace = identify_workspace(drug_syn_strength, drug_red_strength)

    drug_syn_pc = compute_participation_coefficient(drug_synergy, rsn_labels)
    drug_red_pc = compute_participation_coefficient(drug_redundancy, rsn_labels)

    drug_gateway_mask, drug_broadcaster_mask, drug_pc_gradient = classify_gateways_broadcasters(
        drug_syn_pc, drug_red_pc, drug_workspace
    )

    print_and_log(f"{args.drug_task} results:")
    print_and_log(f"  Workspace: {drug_workspace.sum()} regions")
    print_and_log(f"  Gateways: {drug_gateway_mask.sum()}")
    print_and_log(f"  Broadcasters: {drug_broadcaster_mask.sum()}")

    drug_results = {
        'synergy': drug_synergy,
        'redundancy': drug_redundancy,
        'synergy_pc': drug_syn_pc,
        'redundancy_pc': drug_red_pc,
        'workspace_mask': drug_workspace,
        'gateway_mask': drug_gateway_mask,
        'broadcaster_mask': drug_broadcaster_mask,
        'pc_gradient': drug_pc_gradient
    }

    # -------------------------------------------------------------------------
    # Step 5: Compare conditions
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 5: Comparing conditions")
    print_and_log("="*70)

    comparison, summary_df = compare_conditions(control_results, drug_results, rsn_labels, rsn_names)

    print_and_log("\nWorkspace changes:")
    print_and_log(f"  Entered workspace: {comparison['entered_workspace'].sum()} regions")
    print_and_log(f"  Left workspace: {comparison['left_workspace'].sum()} regions")

    print_and_log("\nGateway changes:")
    print_and_log(f"  Became gateway: {comparison['became_gateway'].sum()} regions")
    print_and_log(f"  Lost gateway status: {comparison['lost_gateway'].sum()} regions")
    print_and_log(f"  Stable gateways: {comparison['stable_gateway'].sum()} regions")

    print_and_log("\nBroadcaster changes:")
    print_and_log(f"  Became broadcaster: {comparison['became_broadcaster'].sum()} regions")
    print_and_log(f"  Lost broadcaster status: {comparison['lost_broadcaster'].sum()} regions")
    print_and_log(f"  Stable broadcasters: {comparison['stable_broadcaster'].sum()} regions")

    print_and_log("\nSummary by RSN:")
    print_and_log(summary_df.to_string())

    # -------------------------------------------------------------------------
    # Step 6: Save results with structured naming
    # -------------------------------------------------------------------------
    print_and_log("\n" + "="*70)
    print_and_log("STEP 6: Saving results")
    print_and_log("="*70)

    # File naming: dataset_drugtask_vs_controltask_measure.npy
    prefix = f"{args.dataset}_{args.drug_task}_vs_{control_task}"

    # Save control results
    np.save(comparison_dir / f'{prefix}_control_gateway_mask.npy', control_gateway_mask)
    np.save(comparison_dir / f'{prefix}_control_broadcaster_mask.npy', control_broadcaster_mask)
    np.save(comparison_dir / f'{prefix}_control_workspace_mask.npy', control_workspace)
    np.save(comparison_dir / f'{prefix}_control_synergy_pc.npy', control_syn_pc)
    np.save(comparison_dir / f'{prefix}_control_redundancy_pc.npy', control_red_pc)

    # Save drug results
    np.save(comparison_dir / f'{prefix}_drug_gateway_mask.npy', drug_gateway_mask)
    np.save(comparison_dir / f'{prefix}_drug_broadcaster_mask.npy', drug_broadcaster_mask)
    np.save(comparison_dir / f'{prefix}_drug_workspace_mask.npy', drug_workspace)
    np.save(comparison_dir / f'{prefix}_drug_synergy_pc.npy', drug_syn_pc)
    np.save(comparison_dir / f'{prefix}_drug_redundancy_pc.npy', drug_red_pc)

    # Save comparison
    np.save(comparison_dir / f'{prefix}_became_gateway.npy', comparison['became_gateway'])
    np.save(comparison_dir / f'{prefix}_lost_gateway.npy', comparison['lost_gateway'])
    np.save(comparison_dir / f'{prefix}_became_broadcaster.npy', comparison['became_broadcaster'])
    np.save(comparison_dir / f'{prefix}_lost_broadcaster.npy', comparison['lost_broadcaster'])
    np.save(comparison_dir / f'{prefix}_delta_synergy_pc.npy', comparison['delta_synergy_pc'])
    np.save(comparison_dir / f'{prefix}_delta_redundancy_pc.npy', comparison['delta_redundancy_pc'])

    print_and_log("[OK] Saved condition-specific masks")

    # Save summary
    summary_df.to_csv(comparison_dir / f'{prefix}_summary.csv', index=False)
    print_and_log("[OK] Saved summary CSV")

    # Create detailed region-level comparison
    region_df = pd.DataFrame({
        'ROI': np.arange(1, N_ROIS + 1),  # ROIs 1-454 to match atlas numbering
        'RSN': rsn_labels,
        'RSN_name': [rsn_names[label] for label in rsn_labels],
        'control_workspace': control_workspace,
        'drug_workspace': drug_workspace,
        'control_gateway': control_gateway_mask,
        'drug_gateway': drug_gateway_mask,
        'control_broadcaster': control_broadcaster_mask,
        'drug_broadcaster': drug_broadcaster_mask,
        'entered_workspace': comparison['entered_workspace'],
        'left_workspace': comparison['left_workspace'],
        'became_gateway': comparison['became_gateway'],
        'lost_gateway': comparison['lost_gateway'],
        'became_broadcaster': comparison['became_broadcaster'],
        'lost_broadcaster': comparison['lost_broadcaster'],
        'control_synergy_pc': control_syn_pc,
        'drug_synergy_pc': drug_syn_pc,
        'delta_synergy_pc': comparison['delta_synergy_pc'],
        'control_redundancy_pc': control_red_pc,
        'drug_redundancy_pc': drug_red_pc,
        'delta_redundancy_pc': comparison['delta_redundancy_pc']
    })
    region_df.to_csv(comparison_dir / f'{prefix}_region_level.csv', index=False)
    print_and_log("[OK] Saved region-level comparison")

    # Load control metadata from script 06 to get subject info
    control_metadata_file = control_dir / f'{control_prefix}_gateway_broadcaster_metadata.yaml'
    if control_metadata_file.exists():
        with open(control_metadata_file, 'r') as f:
            control_metadata = yaml.safe_load(f)
        n_control_subjects = control_metadata.get('n_workspace', 'unknown')  # workspace metadata has subject count
    else:
        # Try workspace metadata
        workspace_metadata_file = control_dir / f'{control_prefix}_workspace_metadata.yaml'
        if workspace_metadata_file.exists():
            with open(workspace_metadata_file, 'r') as f:
                control_metadata = yaml.safe_load(f)
            n_control_subjects = control_metadata.get('n_subjects', 'unknown')
        else:
            n_control_subjects = 'unknown (metadata not found)'

    # Metadata
    metadata = {
        'dataset': args.dataset,
        'dataset_name': handler.config['name'],
        'control_task': control_task,
        'drug_task': args.drug_task,
        'workspace_source_dataset': ws_dataset,
        'workspace_source_task': ws_task,
        'workspace_source_path': str(ws_dir),
        'n_control_subjects': n_control_subjects,
        'n_drug_subjects': len(drug_subjects),
        'drug_subjects': drug_subjects,
        'n_rois': N_ROIS,
        'control_workspace': int(control_workspace.sum()),
        'drug_workspace': int(drug_workspace.sum()),
        'control_gateways': int(control_gateway_mask.sum()),
        'drug_gateways': int(drug_gateway_mask.sum()),
        'control_broadcasters': int(control_broadcaster_mask.sum()),
        'drug_broadcasters': int(drug_broadcaster_mask.sum()),
        'became_gateway': int(comparison['became_gateway'].sum()),
        'lost_gateway': int(comparison['lost_gateway'].sum()),
        'became_broadcaster': int(comparison['became_broadcaster'].sum()),
        'lost_broadcaster': int(comparison['lost_broadcaster'].sum()),
        'analysis_date': pd.Timestamp.now().isoformat(),
        'output_directory': str(comparison_dir),
        'file_prefix': prefix,
        'control_results_source': str(control_dir),
        'notes': f'Condition comparison: {args.drug_task} vs {control_task}. '
                 f'Workspace mask from {ws_dataset}/{ws_task}. Control gateway/broadcaster from script 06.'
    }

    with open(comparison_dir / f'{prefix}_metadata.yaml', 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)
    print_and_log("[OK] Saved metadata")

    print_and_log("\n" + "="*70)
    print_and_log("CONDITION COMPARISON COMPLETE")
    print_and_log("="*70)
    print_and_log(f"\nResults saved to: {comparison_dir}")
    print_and_log(f"File prefix: {prefix}_*")

    print_and_log("\nKey findings:")
    print_and_log(f"  - {args.drug_task} workspace: {drug_workspace.sum()} vs {control_task}: {control_workspace.sum()} regions")
    print_and_log(f"  - {args.drug_task} gateways: {drug_gateway_mask.sum()} vs {control_task}: {control_gateway_mask.sum()}")
    print_and_log(f"  - {args.drug_task} broadcasters: {drug_broadcaster_mask.sum()} vs {control_task}: {control_broadcaster_mask.sum()}")
    print_and_log(f"  - Changed gateway status: {comparison['became_gateway'].sum() + comparison['lost_gateway'].sum()} regions")
    print_and_log(f"  - Changed broadcaster status: {comparison['became_broadcaster'].sum() + comparison['lost_broadcaster'].sum()} regions")

    # Save markdown output
    md_file = comparison_dir / f'{prefix}_results.md'
    with open(md_file, 'w') as f:
        f.write('\n'.join(md_output))
    print(f"\n[OK] Markdown output saved to: {md_file}")

if __name__ == '__main__':
    main()
