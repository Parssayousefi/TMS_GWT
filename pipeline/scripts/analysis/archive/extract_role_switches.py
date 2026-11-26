#!/usr/bin/env python3
"""
Extract Workspace Regions That Switch Between Gateway and Broadcaster Roles
============================================================================

This script identifies workspace regions that change their role (gateway ↔ broadcaster)
under drug conditions compared to baseline (placebo/control).

Usage:
    python scripts/analysis/extract_role_switches.py --dataset london_lsd --drug_task LSD
    python scripts/analysis/extract_role_switches.py --dataset basel_mmm --drug_task MDMA

Output:
    - CSV file with all role-switching regions
    - Summary by network
    - Visualization-ready tables
"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
import sys
import argparse

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from load_roi_mapping import load_roi_network_mapping
from dataset_handler import load_dataset_handler
from analysis_loader import get_control_task_name

def load_role_switch_data(dataset, drug_task, control_task, base_results_dir):
    """
    Load role-switching data from comparison results

    Returns
    -------
    dict with keys:
        - became_gateway: array of ROI indices that became gateways
        - became_broadcaster: array of ROI indices that became broadcasters
        - control_workspace: baseline workspace mask
        - drug_workspace: drug workspace mask
        - control_gateway: baseline gateway mask
        - control_broadcaster: baseline broadcaster mask
        - drug_gateway: drug gateway mask
        - drug_broadcaster: drug broadcaster mask
    """
    comparison_dir = base_results_dir / dataset / f'comparison_{drug_task}_vs_{control_task}'
    prefix = f'{dataset}_{drug_task}_vs_{control_task}'

    # Load all relevant masks
    data = {
        'became_gateway': np.load(comparison_dir / f'{prefix}_became_gateway.npy'),
        'became_broadcaster': np.load(comparison_dir / f'{prefix}_became_broadcaster.npy'),
        'lost_gateway': np.load(comparison_dir / f'{prefix}_lost_gateway.npy'),
        'lost_broadcaster': np.load(comparison_dir / f'{prefix}_lost_broadcaster.npy'),
        'control_workspace': np.load(comparison_dir / f'{prefix}_control_workspace_mask.npy'),
        'drug_workspace': np.load(comparison_dir / f'{prefix}_drug_workspace_mask.npy'),
        'control_gateway': np.load(comparison_dir / f'{prefix}_control_gateway_mask.npy'),
        'control_broadcaster': np.load(comparison_dir / f'{prefix}_control_broadcaster_mask.npy'),
        'drug_gateway': np.load(comparison_dir / f'{prefix}_drug_gateway_mask.npy'),
        'drug_broadcaster': np.load(comparison_dir / f'{prefix}_drug_broadcaster_mask.npy'),
    }

    # Load metadata for summary
    with open(comparison_dir / f'{prefix}_metadata.yaml') as f:
        metadata = yaml.safe_load(f)

    return data, metadata

def identify_role_switches(data):
    """
    Identify workspace regions that switched roles

    Returns
    -------
    dict with:
        - broadcaster_to_gateway: ROIs that were broadcasters, became gateways
        - gateway_to_broadcaster: ROIs that were gateways, became broadcasters
    """
    # Broadcaster → Gateway: Was broadcaster in control, became gateway in drug
    broadcaster_to_gateway = data['control_broadcaster'] & data['became_gateway']

    # Gateway → Broadcaster: Was gateway in control, became broadcaster in drug
    gateway_to_broadcaster = data['control_gateway'] & data['became_broadcaster']

    # Get ROI indices
    b2g_rois = np.where(broadcaster_to_gateway)[0]
    g2b_rois = np.where(gateway_to_broadcaster)[0]

    return {
        'broadcaster_to_gateway': b2g_rois,
        'gateway_to_broadcaster': g2b_rois,
        'b2g_mask': broadcaster_to_gateway,
        'g2b_mask': gateway_to_broadcaster
    }

def create_role_switch_table(switches, rsn_labels, rsn_names, dataset, drug, control):
    """
    Create detailed table of role-switching regions
    """
    all_switches = []

    # Broadcaster → Gateway switches
    for roi in switches['broadcaster_to_gateway']:
        all_switches.append({
            'ROI': int(roi) + 1,  # Convert array index (0-453) to ROI label (1-454)
            'Network_ID': int(rsn_labels[roi]),
            'Network': rsn_names[rsn_labels[roi]],
            'Baseline_Role': 'Broadcaster',
            'Drug_Role': 'Gateway',
            'Switch_Direction': 'Broadcaster → Gateway',
            'Interpretation': 'Shifted from broadcasting to gathering information'
        })

    # Gateway → Broadcaster switches
    for roi in switches['gateway_to_broadcaster']:
        all_switches.append({
            'ROI': int(roi) + 1,  # Convert array index (0-453) to ROI label (1-454)
            'Network_ID': int(rsn_labels[roi]),
            'Network': rsn_names[rsn_labels[roi]],
            'Baseline_Role': 'Gateway',
            'Drug_Role': 'Broadcaster',
            'Switch_Direction': 'Gateway → Broadcaster',
            'Interpretation': 'Shifted from gathering to broadcasting information'
        })

    df = pd.DataFrame(all_switches)

    # Sort by network and ROI
    if len(df) > 0:
        df = df.sort_values(['Network_ID', 'ROI']).reset_index(drop=True)

    return df

def summarize_by_network(switches_df):
    """
    Summarize role switches by network
    """
    if len(switches_df) == 0:
        return pd.DataFrame()

    # Count by network and direction
    summary = switches_df.groupby(['Network', 'Switch_Direction']).size().reset_index(name='Count')

    # Pivot for easier reading
    pivot = summary.pivot(index='Network', columns='Switch_Direction', values='Count').fillna(0).astype(int)

    # Add total
    pivot['Total_Switches'] = pivot.sum(axis=1)

    # Sort by total switches
    pivot = pivot.sort_values('Total_Switches', ascending=False)

    return pivot

def main():
    parser = argparse.ArgumentParser(
        description='Extract workspace regions that switch between gateway and broadcaster roles'
    )
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--drug_task', required=True, help='Drug task name')
    parser.add_argument('--project_dir', help='Project directory (optional)')
    parser.add_argument('--output_dir', help='Output directory (optional)')

    args = parser.parse_args()

    # Load dataset handler
    handler = load_dataset_handler(args.dataset, project_dir=args.project_dir)
    control_task = get_control_task_name(args.dataset)

    print("="*70)
    print("WORKSPACE REGION ROLE SWITCHES")
    print("="*70)
    print(f"\nDataset: {handler.config['name']}")
    print(f"Comparison: {args.drug_task} vs {control_task}")
    print(f"\nIdentifying workspace regions that switched between gateway and broadcaster roles...")

    # Setup paths
    base_results_dir = handler.project_dir / 'results' / 'group_analysis'

    # Load role switch data
    print(f"\n{'='*70}")
    print("STEP 1: Loading role switch data")
    print("="*70)

    try:
        data, metadata = load_role_switch_data(
            args.dataset, args.drug_task, control_task, base_results_dir
        )
        print(f"✓ Loaded comparison data")
        print(f"  Control workspace: {data['control_workspace'].sum()} regions")
        print(f"  Drug workspace: {data['drug_workspace'].sum()} regions")
        print(f"  Control gateways: {data['control_gateway'].sum()}")
        print(f"  Control broadcasters: {data['control_broadcaster'].sum()}")
        print(f"  Drug gateways: {data['drug_gateway'].sum()}")
        print(f"  Drug broadcasters: {data['drug_broadcaster'].sum()}")
    except FileNotFoundError as e:
        print(f"\n✗ ERROR: Comparison files not found")
        print(f"  {e}")
        print(f"\nYou need to run the comparison script first:")
        print(f"  python scripts/analysis/06b_compare_conditions_gateways.py \\")
        print(f"      --dataset {args.dataset} --drug_task {args.drug_task}")
        return

    # Load ROI network mapping
    print(f"\n{'='*70}")
    print("STEP 2: Loading ROI network mapping")
    print("="*70)

    rsn_labels, rsn_names, roi_df = load_roi_network_mapping(handler.project_dir)
    print(f"✓ Loaded {len(rsn_labels)} ROI labels across {len(rsn_names)} networks")

    # Identify role switches
    print(f"\n{'='*70}")
    print("STEP 3: Identifying role switches")
    print("="*70)

    switches = identify_role_switches(data)

    n_b2g = len(switches['broadcaster_to_gateway'])
    n_g2b = len(switches['gateway_to_broadcaster'])

    print(f"\n✓ Role switches identified:")
    print(f"  Broadcaster → Gateway: {n_b2g} regions")
    print(f"  Gateway → Broadcaster: {n_g2b} regions")
    print(f"  Total switches: {n_b2g + n_g2b} regions")

    if n_b2g + n_g2b == 0:
        print(f"\n✗ No role switches detected!")
        print(f"  This could mean:")
        print(f"    1. Drug doesn't change workspace architecture")
        print(f"    2. Changes are subtle (PC changes without role switches)")
        print(f"    3. Effects are on non-workspace regions")
        return

    # Create detailed table
    print(f"\n{'='*70}")
    print("STEP 4: Creating detailed tables")
    print("="*70)

    switches_df = create_role_switch_table(
        switches, rsn_labels, rsn_names,
        args.dataset, args.drug_task, control_task
    )

    print(f"✓ Created detailed table with {len(switches_df)} entries")

    # Summarize by network
    network_summary = summarize_by_network(switches_df)

    print(f"\n{'='*70}")
    print("NETWORK SUMMARY")
    print("="*70)
    print(network_summary.to_string())

    # Save results
    print(f"\n{'='*70}")
    print("STEP 5: Saving results")
    print("="*70)

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = base_results_dir / args.dataset / f'role_switches_{args.drug_task}_vs_{control_task}'

    output_dir.mkdir(exist_ok=True, parents=True)

    prefix = f'{args.dataset}_{args.drug_task}_vs_{control_task}'

    # Save detailed table
    switches_csv = output_dir / f'{prefix}_role_switches.csv'
    switches_df.to_csv(switches_csv, index=False)
    print(f"✓ Saved detailed table: {switches_csv}")

    # Save network summary
    if len(network_summary) > 0:
        summary_csv = output_dir / f'{prefix}_role_switches_by_network.csv'
        network_summary.to_csv(summary_csv)
        print(f"✓ Saved network summary: {summary_csv}")

    # Save as numpy arrays (for visualization scripts)
    np.save(output_dir / f'{prefix}_broadcaster_to_gateway_rois.npy',
            switches['broadcaster_to_gateway'])
    np.save(output_dir / f'{prefix}_gateway_to_broadcaster_rois.npy',
            switches['gateway_to_broadcaster'])
    print(f"✓ Saved ROI arrays for visualization")

    # Save metadata
    output_metadata = {
        'dataset': args.dataset,
        'dataset_name': handler.config['name'],
        'drug_task': args.drug_task,
        'control_task': control_task,
        'n_broadcaster_to_gateway': int(n_b2g),
        'n_gateway_to_broadcaster': int(n_g2b),
        'total_switches': int(n_b2g + n_g2b),
        'analysis_date': pd.Timestamp.now().isoformat(),
        'output_directory': str(output_dir)
    }

    with open(output_dir / f'{prefix}_role_switches_metadata.yaml', 'w') as f:
        yaml.dump(output_metadata, f, default_flow_style=False)
    print(f"✓ Saved metadata")

    # Print specific examples
    print(f"\n{'='*70}")
    print("EXAMPLES OF ROLE-SWITCHING REGIONS")
    print("="*70)

    if n_b2g > 0:
        print(f"\n📡 Broadcaster → Gateway (top 5):")
        top_b2g = switches_df[switches_df['Switch_Direction'] == 'Broadcaster → Gateway'].head(5)
        for _, row in top_b2g.iterrows():
            print(f"  ROI {row['ROI']:3d} ({row['Network']:15s}): {row['Interpretation']}")

    if n_g2b > 0:
        print(f"\n🔁 Gateway → Broadcaster (top 5):")
        top_g2b = switches_df[switches_df['Switch_Direction'] == 'Gateway → Broadcaster'].head(5)
        for _, row in top_g2b.iterrows():
            print(f"  ROI {row['ROI']:3d} ({row['Network']:15s}): {row['Interpretation']}")

    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved to: {output_dir}")
    print(f"\nKey files:")
    print(f"  - {prefix}_role_switches.csv               (detailed table)")
    print(f"  - {prefix}_role_switches_by_network.csv    (network summary)")
    print(f"  - {prefix}_role_switches_metadata.yaml     (metadata)")

    print(f"\n💡 Interpretation:")
    if n_b2g > n_g2b:
        print(f"  More regions shifted TO gateway role → {args.drug_task} enhances information gathering")
    elif n_g2b > n_b2g:
        print(f"  More regions shifted TO broadcaster role → {args.drug_task} enhances information distribution")
    else:
        print(f"  Balanced shifts → {args.drug_task} causes workspace reorganization")

if __name__ == '__main__':
    main()
