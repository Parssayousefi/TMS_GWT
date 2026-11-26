#!/usr/bin/env python3
"""
Utility function to load ROI to network mapping from CSV
========================================================

This module provides a centralized function to load the ROI-to-network
mapping that all scripts should use instead of fetching from nilearn directly.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def load_roi_network_mapping(base_dir=None):
    """
    Load ROI to network mapping from CSV file

    Parameters
    ----------
    base_dir : Path or str, optional
        Base directory of the project. If None, attempts to find it automatically.

    Returns
    -------
    rsn_labels : np.ndarray (454,)
        Network ID for each ROI (0-7)
    rsn_names : list of str
        Network names ['Visual', 'Somatomotor', 'DorsAttn', 'SalVentAttn',
                       'Limbic', 'Cont', 'Default', 'Subcortical']
    roi_df : pd.DataFrame
        Full mapping DataFrame with columns: ROI, Network_ID, Network_Name,
        Atlas, Hemisphere, ROI_Label
    """

    # Try to find the CSV file
    if base_dir is None:
        # Try to locate based on this file's location
        current_file = Path(__file__).resolve()
        base_dir = current_file.parent.parent.parent
    else:
        base_dir = Path(base_dir)

    csv_path = base_dir / 'atlases' / 'roi_to_network_mapping.csv'

    if not csv_path.exists():
        raise FileNotFoundError(
            f"ROI mapping CSV not found at: {csv_path}\n"
            f"Please run: python scripts/utils/create_roi_network_mapping.py"
        )

    # Load CSV
    roi_df = pd.read_csv(csv_path)

    # FILTER OUT BACKGROUND (ROI 0, Network_ID = -1)
    # The CSV contains ROI 0 (background) which is not used in analysis
    # The PhiID matrices only contain actual brain ROIs (1-454), not background
    roi_df = roi_df[roi_df['ROI'] != 0].reset_index(drop=True)

    # Verify we have 454 ROIs
    n_rois = len(roi_df)
    if n_rois != 454:
        print(f"WARNING: Expected 454 ROIs, found {n_rois}")

    # Extract network labels array (now indexed 0-453 for the 454 brain ROIs)
    rsn_labels = roi_df['Network_ID'].values

    # Get unique network names in order of ID
    network_id_to_name = dict(zip(roi_df['Network_ID'], roi_df['Network_Name']))
    unique_ids = sorted(network_id_to_name.keys())
    rsn_names = [network_id_to_name[i] for i in unique_ids]

    return rsn_labels, rsn_names, roi_df

def get_rsn_names():
    """
    Get network names from CSV in proper order

    Returns
    -------
    list of str
        Network names ordered by Network_ID (0-7)
        Dynamically loaded from the CSV file
    """
    _, rsn_names, _ = load_roi_network_mapping()
    return rsn_names

if __name__ == '__main__':
    # Test the function
    rsn_labels, rsn_names, roi_df = load_roi_network_mapping()

    print("="*70)
    print("ROI NETWORK MAPPING LOADED")
    print("="*70)
    print(f"\nTotal ROIs: {len(rsn_labels)}")
    print(f"\nNetwork names: {rsn_names}")
    print(f"\nROIs per network:")
    for i, name in enumerate(rsn_names):
        count = np.sum(rsn_labels == i)
        print(f"  {i}: {name:15s} - {count} ROIs")

    print(f"\nFirst 5 ROIs:")
    print(roi_df.head())

    print(f"\nLast 5 ROIs:")
    print(roi_df.tail())
