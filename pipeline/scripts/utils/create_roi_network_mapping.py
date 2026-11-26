#!/usr/bin/env python3
"""
Create ROI to Network Name Mapping CSV
=======================================

This script creates a CSV file that maps each of the 454 ROIs to their
corresponding network name, using nilearn's Schaefer 2018 atlas.

Output: roi_to_network_mapping.csv with columns:
    - ROI: ROI index (0-454)
        * 0 = Background (not used in analysis)
        * 1-400 = Schaefer 400 cortical ROIs
        * 401-454 = Tian 54 subcortical ROIs
    - Network_ID: Network ID (-1 for background, 0-7 for networks)
    - Network_Name: Full network name
    - Atlas: Which atlas the ROI comes from (Schaefer400 or Tian54)
    - Hemisphere: LH, RH, or Bilateral (for subcortical)
    - ROI_Label: Full label from the atlas
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent.parent
OUTPUT_FILE = BASE_DIR / 'atlases' / 'roi_to_network_mapping.csv'
OUTPUT_FILE.parent.mkdir(exist_ok=True, parents=True)

# Parameters
N_ROIS = 455  # 1 background + 454 brain regions
N_CORTICAL = 400  # Schaefer
N_SUBCORTICAL = 54  # Tian

def get_roi_network_mapping():
    """
    Create detailed mapping of each ROI to its network

    Returns
    -------
    pd.DataFrame
        Mapping with columns: ROI, Network_ID, Network_Name, Atlas, Hemisphere, ROI_Label
    """

    # Network ID to name mapping
    network_names = {
        0: 'Visual',
        1: 'Somatomotor',
        2: 'DorsAttn',
        3: 'SalVentAttn',
        4: 'Limbic',
        5: 'Cont',
        6: 'Default',
        7: 'Subcortical'
    }

    # Initialize lists for DataFrame
    roi_indices = []
    network_ids = []
    network_names_list = []
    atlas_list = []
    hemisphere_list = []
    roi_labels = []

    # -------------------------------------------------------------------------
    # Part 0: Background ROI - index 0
    # -------------------------------------------------------------------------
    print("Adding background ROI at index 0...")
    roi_indices.append(0)
    network_ids.append(-1)
    network_names_list.append('Background')
    atlas_list.append('Background')
    hemisphere_list.append('Unknown')
    roi_labels.append('Background')
    print("✓ Added background ROI")

    # -------------------------------------------------------------------------
    # Part 1: Cortical ROIs (Schaefer 400) - indices 1-400
    # -------------------------------------------------------------------------
    print("\nLoading Schaefer 400 parcellation from nilearn...")

    try:
        from nilearn import datasets

        # Fetch Schaefer 400 with 7 networks
        atlas = datasets.fetch_atlas_schaefer_2018(n_rois=400, yeo_networks=7, verbose=0)
        labels = atlas['labels']

        # Network abbreviation to ID mapping
        network_map = {
            'Vis': 0,        # Visual
            'SomMot': 1,     # Somatomotor
            'DorsAttn': 2,   # Dorsal Attention
            'SalVentAttn': 3,  # Salience/Ventral Attention
            'Limbic': 4,     # Limbic
            'Cont': 5,       # Control (Frontoparietal)
            'Default': 6     # Default Mode
        }

        print(f"✓ Loaded {len(labels)} labels from nilearn")

        # Check if first label is background and skip it
        first_label = labels[0].decode() if isinstance(labels[0], bytes) else labels[0]
        if first_label == 'Background':
            print("  Skipping background label from nilearn")
            labels = labels[1:]  # Skip background

        print(f"✓ Processing {len(labels)} cortical ROI labels")

        # Process each cortical ROI (indices 1-400, not 0-399)
        for i, label in enumerate(labels):
            # Decode label if bytes
            label_str = label.decode() if isinstance(label, bytes) else label

            # Parse label: format is "7Networks_LH_Vis_1"
            parts = label_str.split('_')

            roi_indices.append(i + 1)  # ROI 1-400, not 0-399
            roi_labels.append(label_str)
            atlas_list.append('Schaefer400')

            # Extract hemisphere (LH or RH)
            if len(parts) >= 2:
                hemisphere_list.append(parts[1])  # 'LH' or 'RH'
            else:
                hemisphere_list.append('Unknown')

            # Extract network abbreviation and map to ID
            if len(parts) >= 3:
                network_abbr = parts[2]
                if network_abbr in network_map:
                    net_id = network_map[network_abbr]
                    network_ids.append(net_id)
                    network_names_list.append(network_names[net_id])
                else:
                    print(f"WARNING: Unknown network '{network_abbr}' for ROI {i+1}")
                    network_ids.append(-1)
                    network_names_list.append('Unknown')
            else:
                print(f"WARNING: Cannot parse label '{label_str}' for ROI {i+1}")
                network_ids.append(-1)
                network_names_list.append('Unknown')

    except ImportError as e:
        print("ERROR: nilearn library not found")
        print("Please install nilearn: pip install nilearn")
        raise RuntimeError("Cannot proceed without nilearn") from e

    # -------------------------------------------------------------------------
    # Part 2: Subcortical ROIs (Tian 54) - indices 401-454
    # -------------------------------------------------------------------------
    print(f"\nAdding {N_SUBCORTICAL} subcortical ROIs (Tian atlas)...")

    # Tian Subcortex S4 3T labels (from official atlas repository)
    # Source: https://github.com/yetianmed/subcortex/blob/master/Group-Parcellation/3T/Subcortex-Only/Tian_Subcortex_S4_3T_label.txt
    # Scale 4 has 54 ROIs (27 per hemisphere)
    tian_labels = [
        'HIP-head-m1-rh', 'HIP-head-m2-rh', 'THA-VAip-rh', 'THA-VAia-rh',
        'HIP-head-l-rh', 'HIP-body-rh', 'HIP-tail-rh',
        'THA-VPm-rh', 'THA-VPl-rh', 'THA-VAs-rh',
        'THA-DAm-rh', 'THA-DAl-rh',
        'PUT-VA-rh', 'PUT-DA-rh', 'PUT-VP-rh', 'PUT-DP-rh',
        'CAU-VA-rh', 'CAU-DA-rh', 'CAU-body-rh', 'CAU-tail-rh',
        'lAMY-rh', 'mAMY-rh', 'THA-DP-rh',
        'NAc-shell-rh', 'NAc-core-rh', 'pGP-rh', 'aGP-rh',
        'HIP-head-m1-lh', 'HIP-head-m2-lh', 'THA-VAip-lh', 'THA-VAia-lh',
        'HIP-head-l-lh', 'HIP-body-lh', 'HIP-tail-lh',
        'THA-VPm-lh', 'THA-VPl-lh', 'THA-VAs-lh',
        'THA-DAm-lh', 'THA-DAl-lh',
        'PUT-VA-lh', 'PUT-DA-lh', 'PUT-VP-lh', 'PUT-DP-lh',
        'CAU-VA-lh', 'CAU-DA-lh', 'CAU-body-lh', 'CAU-tail-lh',
        'lAMY-lh', 'mAMY-lh', 'THA-DP-lh',
        'NAc-shell-lh', 'NAc-core-lh', 'pGP-lh', 'aGP-lh'
    ]

    # Verify we have exactly 54 labels
    if len(tian_labels) != N_SUBCORTICAL:
        print(f"WARNING: Expected {N_SUBCORTICAL} Tian labels but got {len(tian_labels)}")

    for i in range(N_SUBCORTICAL):
        roi_idx = N_CORTICAL + i + 1  # ROI 401-454, not 400-453
        roi_indices.append(roi_idx)
        network_ids.append(7)  # Subcortical network ID
        network_names_list.append('Subcortical')
        atlas_list.append('Tian54')

        # Extract hemisphere from label (ends with -rh or -lh)
        label = tian_labels[i] if i < len(tian_labels) else f'Unknown_{i+1}'
        if label.endswith('-rh'):
            hemisphere_list.append('RH')
        elif label.endswith('-lh'):
            hemisphere_list.append('LH')
        else:
            hemisphere_list.append('Bilateral')

        roi_labels.append(label)

    print(f"✓ Added subcortical ROIs (indices {N_CORTICAL + 1}-{N_ROIS-1})")
    print(f"  Structures: Hippocampus, Thalamus, Putamen, Caudate, Amygdala, NAc, GP")

    # -------------------------------------------------------------------------
    # Create DataFrame
    # -------------------------------------------------------------------------
    df = pd.DataFrame({
        'ROI': roi_indices,
        'Network_ID': network_ids,
        'Network_Name': network_names_list,
        'Atlas': atlas_list,
        'Hemisphere': hemisphere_list,
        'ROI_Label': roi_labels
    })

    return df

def main():
    """Main execution"""

    print("="*70)
    print("CREATE ROI TO NETWORK MAPPING")
    print("="*70)
    print()

    # Get mapping
    df = get_roi_network_mapping()

    # Print summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nTotal ROIs: {len(df)}")
    print(f"\nROIs per network:")
    for net_id in sorted(df['Network_ID'].unique()):
        net_name = df[df['Network_ID'] == net_id]['Network_Name'].iloc[0]
        count = (df['Network_ID'] == net_id).sum()
        print(f"  {net_id}: {net_name:15s} - {count:3d} ROIs")

    print(f"\nROIs per atlas:")
    for atlas in df['Atlas'].unique():
        count = (df['Atlas'] == atlas).sum()
        print(f"  {atlas:15s} - {count:3d} ROIs")

    # Save to CSV
    print("\n" + "="*70)
    print("SAVING")
    print("="*70)

    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Saved to: {OUTPUT_FILE}")

    # Show first and last few rows
    print("\nFirst 10 rows:")
    print(df.head(10).to_string(index=False))

    print("\nLast 10 rows:")
    print(df.tail(10).to_string(index=False))

    print("\n" + "="*70)
    print("COMPLETE")
    print("="*70)
    print(f"\nMapping file created: {OUTPUT_FILE}")
    print(f"Total ROIs: {len(df)}")
    print("\nYou can now use this CSV to:")
    print("  - Look up which network any ROI belongs to")
    print("  - Filter ROIs by network")
    print("  - Map ROI indices to readable network names")

if __name__ == '__main__':
    main()
