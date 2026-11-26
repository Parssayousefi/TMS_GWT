#!/usr/bin/env python3
"""
Extract STATISTICALLY SIGNIFICANT Role Switches
================================================

This script identifies workspace regions that:
1. Switch roles (gateway ↔ broadcaster) based on group-average PC rankings
2. Show SIGNIFICANT participation coefficient changes (p_FDR < 0.05)

This combines descriptive role classification with inferential statistics.

Usage:
    python scripts/analysis/extract_significant_role_switches.py --dataset london_lsd --drug_task LSD
"""

import numpy as np
import pandas as pd
import yaml
from pathlib import Path
import sys
import argparse

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from dataset_handler import load_dataset_handler
from analysis_loader import get_control_task_name

def main():
    parser = argparse.ArgumentParser(
        description='Extract statistically significant role switches'
    )
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--drug_task', required=True, help='Drug task name')
    parser.add_argument('--project_dir', help='Project directory (optional)')
    parser.add_argument('--alpha', type=float, default=0.05, help='Significance threshold (default: 0.05)')

    args = parser.parse_args()

    # Load dataset handler
    handler = load_dataset_handler(args.dataset, project_dir=args.project_dir)
    control_task = get_control_task_name(args.dataset)

    print("="*70)
    print("STATISTICALLY SIGNIFICANT ROLE SWITCHES")
    print("="*70)
    print(f"\nDataset: {handler.config['name']}")
    print(f"Comparison: {args.drug_task} vs {control_task}")
    print(f"Significance threshold: p_FDR < {args.alpha}")

    # Setup paths
    base_results_dir = handler.project_dir / 'results' / 'group_analysis'
    role_switches_dir = base_results_dir / args.dataset / f'role_switches_{args.drug_task}_vs_{control_task}'
    paired_stats_dir = base_results_dir / args.dataset / f'paired_stats_{args.drug_task}_vs_{control_task}'

    # Load role switches (descriptive)
    print(f"\n{'='*70}")
    print("STEP 1: Loading descriptive role switches")
    print("="*70)

    try:
        switches_csv = role_switches_dir / f'{args.dataset}_{args.drug_task}_vs_{control_task}_role_switches.csv'
        switches_df = pd.read_csv(switches_csv)
        print(f"✓ Loaded {len(switches_df)} descriptive role switches")
    except FileNotFoundError:
        print(f"✗ ERROR: Role switches file not found")
        print(f"\nRun first: python scripts/analysis/extract_role_switches.py --dataset {args.dataset} --drug_task {args.drug_task}")
        return

    if len(switches_df) == 0:
        print("No role switches detected")
        return

    # Load paired statistics
    print(f"\n{'='*70}")
    print("STEP 2: Loading paired statistics")
    print("="*70)

    try:
        synergy_stats = pd.read_csv(paired_stats_dir / f'{args.dataset}_{args.drug_task}_vs_{control_task}_synergy_pc_stats.csv')
        redundancy_stats = pd.read_csv(paired_stats_dir / f'{args.dataset}_{args.drug_task}_vs_{control_task}_redundancy_pc_stats.csv')
        print(f"✓ Loaded synergy PC statistics ({len(synergy_stats)} ROIs)")
        print(f"✓ Loaded redundancy PC statistics ({len(redundancy_stats)} ROIs)")
    except FileNotFoundError:
        print(f"✗ ERROR: Paired statistics not found")
        print(f"\nRun first: python scripts/analysis/07_paired_statistical_testing.py --dataset {args.dataset} --drug_task {args.drug_task}")
        return

    # Identify significant switches
    print(f"\n{'='*70}")
    print("STEP 3: Identifying statistically significant switches")
    print("="*70)

    significant_switches = []

    for _, row in switches_df.iterrows():
        roi = row['ROI']
        direction = row['Switch_Direction']

        # For Broadcaster → Gateway: Check if synergy PC significantly increased
        if 'Gateway' in direction.split('→')[1]:
            matching_rows = synergy_stats[synergy_stats['ROI'] == roi]
            if len(matching_rows) == 0:
                print(f"WARNING: ROI {roi} not found in synergy PC statistics. Skipping.")
                continue
            stats_row = matching_rows.iloc[0]
            measure = 'synergy_pc'
            pc_type = 'Synergy'

        # For Gateway → Broadcaster: Check if redundancy PC significantly increased
        else:
            matching_rows = redundancy_stats[redundancy_stats['ROI'] == roi]
            if len(matching_rows) == 0:
                print(f"WARNING: ROI {roi} not found in redundancy PC statistics. Skipping.")
                continue
            stats_row = matching_rows.iloc[0]
            measure = 'redundancy_pc'
            pc_type = 'Redundancy'

        # Check significance
        is_significant = stats_row['p_corrected'] < args.alpha

        switch_info = {
            'ROI': int(roi),
            'Network': row['Network'],
            'Baseline_Role': row['Baseline_Role'],
            'Drug_Role': row['Drug_Role'],
            'Switch_Direction': direction,
            'Mean_PC_Change': stats_row['mean_diff'],
            'Cohen_d': stats_row['cohen_d'],
            'p_value': stats_row['p_value'],
            'p_corrected': stats_row['p_corrected'],
            'Significant': is_significant,
            'CI_lower': stats_row['ci_lower'],
            'CI_upper': stats_row['ci_upper'],
            'PC_Type': pc_type
        }

        significant_switches.append(switch_info)

    sig_df = pd.DataFrame(significant_switches)

    # Summary
    n_significant = sig_df['Significant'].sum()
    n_total = len(sig_df)

    print(f"\nDescriptive role switches: {n_total}")
    print(f"Statistically significant: {n_significant} ({100*n_significant/n_total:.1f}%)")
    print(f"Non-significant (descriptive only): {n_total - n_significant}")

    if n_significant > 0:
        print(f"\n{'='*70}")
        print("STATISTICALLY SIGNIFICANT SWITCHES")
        print("="*70)

        sig_only = sig_df[sig_df['Significant']]
        for _, row in sig_only.iterrows():
            print(f"\n✓ ROI {row['ROI']} ({row['Network']}): {row['Switch_Direction']}")
            print(f"  Mean Δ{row['PC_Type']} PC: {row['Mean_PC_Change']:+.6f}")
            print(f"  Cohen's d: {row['Cohen_d']:+.3f}")
            print(f"  p_FDR: {row['p_corrected']:.4f}")
            print(f"  95% CI: [{row['CI_lower']:.6f}, {row['CI_upper']:.6f}]")

    # Non-significant
    print(f"\n{'='*70}")
    print("NON-SIGNIFICANT SWITCHES (Descriptive Only)")
    print("="*70)

    non_sig = sig_df[~sig_df['Significant']]
    if len(non_sig) > 0:
        print(f"\nThese {len(non_sig)} regions show role switches in group-average rankings")
        print(f"but do NOT show significant PC changes at the individual level:")

        for _, row in non_sig.iterrows():
            print(f"\n✗ ROI {row['ROI']} ({row['Network']}): {row['Switch_Direction']}")
            print(f"  Mean Δ{row['PC_Type']} PC: {row['Mean_PC_Change']:+.6f}")
            print(f"  Cohen's d: {row['Cohen_d']:+.3f}")
            print(f"  p_FDR: {row['p_corrected']:.4f} (n.s.)")
    else:
        print("\nAll switches are statistically significant!")

    # Save results
    print(f"\n{'='*70}")
    print("STEP 4: Saving results")
    print("="*70)

    output_dir = base_results_dir / args.dataset / f'significant_role_switches_{args.drug_task}_vs_{control_task}'
    output_dir.mkdir(exist_ok=True, parents=True)

    prefix = f'{args.dataset}_{args.drug_task}_vs_{control_task}'

    # Save full table
    sig_df.to_csv(output_dir / f'{prefix}_role_switches_with_statistics.csv', index=False)
    print(f"✓ Saved full table with statistics")

    # Save significant only
    if n_significant > 0:
        sig_only_df = sig_df[sig_df['Significant']]
        sig_only_df.to_csv(output_dir / f'{prefix}_significant_role_switches_only.csv', index=False)
        print(f"✓ Saved significant switches only")

    # Save metadata
    metadata = {
        'dataset': args.dataset,
        'drug_task': args.drug_task,
        'control_task': control_task,
        'alpha': args.alpha,
        'n_total_switches': int(n_total),
        'n_significant_switches': int(n_significant),
        'n_non_significant': int(n_total - n_significant),
        'pct_significant': float(100 * n_significant / n_total) if n_total > 0 else 0,
        'analysis_date': pd.Timestamp.now().isoformat()
    }

    with open(output_dir / f'{prefix}_metadata.yaml', 'w') as f:
        yaml.dump(metadata, f, default_flow_style=False)
    print(f"✓ Saved metadata")

    # Summary
    print(f"\n{'='*70}")
    print("INTERPRETATION & RECOMMENDATIONS")
    print("="*70)

    if n_significant == 0:
        print(f"\n⚠️  NO STATISTICALLY SIGNIFICANT ROLE SWITCHES DETECTED")
        print(f"\nWhat this means:")
        print(f"  - {n_total} regions show role switches in group-average rankings")
        print(f"  - BUT these switches are NOT reliable at the individual subject level")
        print(f"  - Likely due to small effect sizes and/or high inter-subject variability")
        print(f"\nRecommendations:")
        print(f"  1. Report as EXPLORATORY: 'Descriptive analysis suggests...'")
        print(f"  2. Focus on NETWORK-LEVEL PC changes (from paired stats network summary)")
        print(f"  3. Do NOT claim 'X regions significantly switched roles'")
        print(f"  4. Consider: Larger sample size may reveal significant switches")

    elif n_significant < n_total:
        print(f"\n⚠️  PARTIAL SIGNIFICANCE: {n_significant}/{n_total} switches are significant")
        print(f"\nRecommendations:")
        print(f"  1. Report ONLY the {n_significant} significant switch(es)")
        print(f"  2. Mention {n_total - n_significant} additional descriptive switches in supplement")
        print(f"  3. Use language: '{n_significant} workspace regions showed significant role switches'")

    else:
        print(f"\n✓ ALL ROLE SWITCHES ARE STATISTICALLY SIGNIFICANT")
        print(f"\nRecommendations:")
        print(f"  1. Strong finding! All {n_total} switches are individually reliable")
        print(f"  2. Report with confidence: '{n_total} workspace regions significantly switched roles'")
        print(f"  3. Highlight in main results (not exploratory)")

    print(f"\n{'='*70}")
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved to: {output_dir}")

if __name__ == '__main__':
    main()
