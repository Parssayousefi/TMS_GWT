#!/usr/bin/env python3
"""
Analysis Data Loader - Dataset-Aware Loading for Group-Level Analyses

This module provides functions to load ΦID results across multiple subjects
and conditions, handling different dataset structures automatically.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import sys

sys.path.insert(0, str(Path(__file__).parent))
from dataset_handler import load_dataset_handler
from qc_loader import QCExclusions


def load_phiid_matrices_for_scan(
    dataset: str,
    subject: str,
    task: str,
    run: int,
    session: Optional[str] = None,
    results_dir: Optional[Path] = None,
    synergy_suffix: str = 'synergy'
) -> Optional[Dict[str, np.ndarray]]:
    """
    Load ΦID matrices for a single scan

    Parameters:
    -----------
    dataset : str
        Dataset name
    subject : str
        Subject ID (without sub- prefix)
    task : str
        Task name
    run : int
        Run number
    session : str, optional
        Session ID (without ses- prefix)
    results_dir : Path, optional
        Results directory (defaults to project_dir/results/phiid_matrices)
    synergy_suffix : str
        File suffix for synergy matrices (default: 'synergy').
        Use 'synergy_total' to load total source synergy (Eq. 5) files.

    Returns:
    --------
    dict or None
        Dictionary with keys: 'synergy', 'redundancy', 'phi', 'phi_r'
        Returns None if files not found
    """
    # Load dataset handler to get project_dir
    handler = load_dataset_handler(dataset)

    if results_dir is None:
        results_dir = handler.project_dir / 'results' / 'phiid_matrices'

    # Map logical measure names to file suffixes
    measure_suffixes = {
        'synergy': synergy_suffix,
        'redundancy': 'redundancy',
        'phi': 'phi',
        'phi_r': 'phi_r'
    }

    # Use dataset handler to generate standardized filenames with dataset prefix
    files = {}
    for measure, suffix in measure_suffixes.items():
        filename = handler.get_output_filename(
            subject=subject,
            task=task,
            run=run,
            session=session,
            suffix=suffix,
            extension='npy',
            include_dataset=True
        )
        files[measure] = results_dir / filename

    # Check all files exist
    if not all(f.exists() for f in files.values()):
        return None

    # Load all matrices
    matrices = {key: np.load(str(path)) for key, path in files.items()}

    return matrices


def discover_processed_scans(
    dataset: str,
    task_filter: Optional[str] = None,
    results_dir: Optional[Path] = None,
    synergy_suffix: str = 'synergy'
) -> List[Dict]:
    """
    Discover all processed scans for a dataset

    Parameters:
    -----------
    dataset : str
        Dataset name
    task_filter : str, optional
        Logical condition/task name (e.g., 'PCB', 'LSD', 'placebo', 'ketamine').
        This is the KEY from the dataset config's 'tasks' section, NOT the
        BIDS task name. For session-based datasets where condition is encoded
        in the session name, filtering matches the session field.
    results_dir : Path, optional
        Results directory
    synergy_suffix : str
        File suffix for synergy matrices (default: 'synergy').

    Returns:
    --------
    list of dict
        Each dict contains: subject, session, task, run, files_exist
    """
    handler = load_dataset_handler(dataset)

    if results_dir is None:
        results_dir = handler.project_dir / 'results' / 'phiid_matrices'

    # Optional QC gate: when the dataset config sets
    # `qc.honor_exclusions_in_analysis: true`, skip scans flagged as
    # EXCLUDE in the merged QCExclusions store (motion + ART). Default is
    # off so all other datasets retain their existing behaviour.
    honor_qc = handler.config.get('qc', {}).get(
        'honor_exclusions_in_analysis', False
    )
    qc = QCExclusions(dataset) if honor_qc else None

    # Resolve task_filter: determine whether it maps to a BIDS task name
    # or a session name (for datasets where condition = session)
    filter_by_session = None
    filter_by_bids_task = None

    if task_filter:
        tasks_config = handler.config.get('tasks', {})
        has_sessions = handler.config.get('bids_structure', {}).get('has_sessions', False)

        if task_filter in tasks_config:
            # task_filter is a logical condition name (config key)
            task_info = tasks_config[task_filter]
            bids_task_name = task_info.get('bids_task_name', task_filter) if isinstance(task_info, dict) else task_filter

            if has_sessions and bids_task_name != task_filter:
                # Condition is encoded in session name, not BIDS task
                # e.g., config key='PCB', bids_task_name='rest', session='PCB'
                filter_by_session = task_filter
                filter_by_bids_task = bids_task_name
            else:
                # Condition IS the BIDS task name
                filter_by_bids_task = task_filter
        else:
            # Not a config key — treat as literal BIDS task name
            filter_by_bids_task = task_filter

    # Get all scans from dataset
    all_scans = handler.discover_all_scans()

    # Check which scans have been processed
    processed_scans = []

    n_qc_dropped = 0
    for scan in all_scans:
        # Apply filters
        if filter_by_bids_task and scan['task'] != filter_by_bids_task:
            continue
        if filter_by_session and scan.get('session') != filter_by_session:
            continue

        # QC gate (config-opt-in via qc.honor_exclusions_in_analysis)
        if qc is not None and qc.is_excluded(
            scan['subject'], scan.get('session'), scan['task']
        ):
            n_qc_dropped += 1
            print(f"  QC drop: sub-{scan['subject']} "
                  f"task-{scan['task']}"
                  f"{' ses-' + scan['session'] if scan.get('session') else ''} "
                  f"(status={qc.get_status(scan['subject'], scan.get('session'), scan['task'])})")
            continue

        # Check if ΦID matrices exist for this scan
        matrices = load_phiid_matrices_for_scan(
            dataset,
            scan['subject'],
            scan['task'],
            scan['run'],
            scan.get('session'),
            results_dir,
            synergy_suffix=synergy_suffix
        )

        if matrices is not None:
            scan_info = {
                'subject': scan['subject'],
                'session': scan.get('session'),
                'task': scan['task'],
                'run': scan['run'],
                'files_exist': True
            }
            processed_scans.append(scan_info)

    if qc is not None and n_qc_dropped > 0:
        print(f"  QC summary: {n_qc_dropped} scan(s) dropped, "
              f"{len(processed_scans)} retained for {dataset}"
              f"{' task=' + task_filter if task_filter else ''}")

    return processed_scans


def load_group_matrices(
    dataset: str,
    task: str,
    measure: str = 'synergy',
    results_dir: Optional[Path] = None,
    synergy_suffix: str = 'synergy'
) -> Tuple[np.ndarray, List[str]]:
    """
    Load matrices for all subjects in a dataset for a specific task

    Parameters:
    -----------
    dataset : str
        Dataset name
    task : str
        Task name
    measure : str
        Which measure to load: 'synergy', 'redundancy', 'phi', or 'phi_r'
    results_dir : Path, optional
        Results directory
    synergy_suffix : str
        File suffix for synergy matrices (default: 'synergy').

    Returns:
    --------
    matrices : np.ndarray
        Shape (n_subjects, n_rois, n_rois)
    subjects : list of str
        Subject IDs that were successfully loaded
    """
    # Discover processed scans
    processed_scans = discover_processed_scans(dataset, task_filter=task, results_dir=results_dir, synergy_suffix=synergy_suffix)

    if len(processed_scans) == 0:
        raise ValueError(f"No processed scans found for dataset={dataset}, task={task}")

    # Load matrices for each subject
    matrices = []
    subjects = []

    for scan in processed_scans:
        mats = load_phiid_matrices_for_scan(
            dataset,
            scan['subject'],
            scan['task'],
            scan['run'],
            scan.get('session'),
            results_dir,
            synergy_suffix=synergy_suffix
        )

        if mats and measure in mats:
            matrices.append(mats[measure])
            subjects.append(scan['subject'])

    if len(matrices) == 0:
        raise ValueError(f"No matrices loaded for dataset={dataset}, task={task}, measure={measure}")

    matrices = np.array(matrices)

    return matrices, subjects


def load_group_all_measures(
    dataset: str,
    task: str,
    results_dir: Optional[Path] = None,
    synergy_suffix: str = 'synergy'
) -> Tuple[Dict[str, np.ndarray], List[str]]:
    """
    Load all ΦID measures for all subjects in a dataset for a specific task

    Parameters:
    -----------
    dataset : str
        Dataset name
    task : str
        Task name
    results_dir : Path, optional
        Results directory
    synergy_suffix : str
        File suffix for synergy matrices (default: 'synergy').
        Use 'synergy_total' to load total source synergy (Eq. 5) files.

    Returns:
    --------
    all_matrices : dict
        Keys: 'synergy', 'redundancy', 'phi', 'phi_r'
        Values: np.ndarray of shape (n_subjects, n_rois, n_rois)
    subjects : list of str
        Subject IDs that were successfully loaded
    """
    # Discover processed scans
    processed_scans = discover_processed_scans(dataset, task_filter=task, results_dir=results_dir, synergy_suffix=synergy_suffix)

    if len(processed_scans) == 0:
        raise ValueError(f"No processed scans found for dataset={dataset}, task={task}")

    # Initialize storage
    all_matrices = {
        'synergy': [],
        'redundancy': [],
        'phi': [],
        'phi_r': []
    }
    subjects = []

    # Load matrices for each subject
    for scan in processed_scans:
        mats = load_phiid_matrices_for_scan(
            dataset,
            scan['subject'],
            scan['task'],
            scan['run'],
            scan.get('session'),
            results_dir,
            synergy_suffix=synergy_suffix
        )

        if mats:
            for measure in all_matrices.keys():
                if measure in mats:
                    all_matrices[measure].append(mats[measure])

            subjects.append(scan['subject'])

    if len(subjects) == 0:
        raise ValueError(f"No matrices loaded for dataset={dataset}, task={task}")

    # Convert to arrays
    for measure in all_matrices.keys():
        all_matrices[measure] = np.array(all_matrices[measure])

    return all_matrices, subjects


def get_control_task_name(dataset: str) -> str:
    """
    Get the control/placebo task name for a dataset

    Parameters:
    -----------
    dataset : str
        Dataset name

    Returns:
    --------
    str
        Name of control/placebo task
    """
    handler = load_dataset_handler(dataset)

    # Find control task
    for task_name, task_info in handler.config['tasks'].items():
        if task_info.get('is_control', False):
            return task_name

    raise ValueError(f"No control task found for dataset {dataset}")


def get_drug_task_names(dataset: str) -> List[str]:
    """
    Get the drug/experimental task names for a dataset

    Parameters:
    -----------
    dataset : str
        Dataset name

    Returns:
    --------
    list of str
        Names of drug/experimental tasks
    """
    handler = load_dataset_handler(dataset)

    drug_tasks = []
    for task_name, task_info in handler.config['tasks'].items():
        if not task_info.get('is_control', False):
            drug_tasks.append(task_name)

    return drug_tasks


# Convenience functions for specific comparisons

def load_baseline_group_data(dataset: str, results_dir: Optional[Path] = None):
    """
    Load group-level data for baseline/control condition

    Parameters:
    -----------
    dataset : str
        Dataset name
    results_dir : Path, optional
        Results directory

    Returns:
    --------
    all_matrices : dict
        Keys: 'synergy', 'redundancy', 'phi', 'phi_r'
        Values: np.ndarray of shape (n_subjects, n_rois, n_rois)
    subjects : list of str
        Subject IDs
    control_task : str
        Name of control task
    """
    control_task = get_control_task_name(dataset)
    all_matrices, subjects = load_group_all_measures(dataset, control_task, results_dir)

    return all_matrices, subjects, control_task


def load_drug_comparison_data(
    dataset: str,
    drug_task: str,
    results_dir: Optional[Path] = None
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], List[str], List[str]]:
    """
    Load group-level data for control vs drug comparison

    IMPORTANT: This function filters to MATCHED SUBJECTS ONLY.
    Only subjects with data in BOTH conditions are returned.

    Parameters:
    -----------
    dataset : str
        Dataset name
    drug_task : str
        Drug task name
    results_dir : Path, optional
        Results directory

    Returns:
    --------
    control_matrices : dict
        Control condition matrices (matched subjects only)
    drug_matrices : dict
        Drug condition matrices (matched subjects only)
    matched_subjects : list of str
        Subjects present in BOTH conditions
    matched_subjects : list of str
        Same as above (for backward compatibility)
    """
    control_task = get_control_task_name(dataset)

    control_matrices, control_subjects = load_group_all_measures(
        dataset, control_task, results_dir
    )
    drug_matrices, drug_subjects = load_group_all_measures(
        dataset, drug_task, results_dir
    )

    # Find matched subjects (present in both conditions)
    control_set = set(control_subjects)
    drug_set = set(drug_subjects)
    matched_subjects = sorted(control_set & drug_set)

    # Filter to matched subjects only
    if len(matched_subjects) < len(control_subjects) or len(matched_subjects) < len(drug_subjects):
        print(f"WARNING: Subject mismatch detected!")
        print(f"  Control subjects: {len(control_subjects)}")
        print(f"  Drug subjects: {len(drug_subjects)}")
        print(f"  Matched subjects: {len(matched_subjects)}")

        # Identify mismatches
        control_only = control_set - drug_set
        drug_only = drug_set - control_set
        if control_only:
            print(f"  Subjects with {control_task} but not {drug_task}: {sorted(control_only)}")
        if drug_only:
            print(f"  Subjects with {drug_task} but not {control_task}: {sorted(drug_only)}")
        print(f"  Filtering to {len(matched_subjects)} matched subjects")

    # Get indices for matched subjects
    control_indices = [control_subjects.index(s) for s in matched_subjects]
    drug_indices = [drug_subjects.index(s) for s in matched_subjects]

    # Filter all matrices
    filtered_control = {}
    filtered_drug = {}

    for key in control_matrices.keys():
        filtered_control[key] = control_matrices[key][control_indices]
        filtered_drug[key] = drug_matrices[key][drug_indices]

    return filtered_control, filtered_drug, matched_subjects, matched_subjects


# Testing/CLI
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Test analysis data loader')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--task', help='Task name (optional)')
    args = parser.parse_args()

    print(f"Testing analysis loader for dataset: {args.dataset}")

    # Find control task
    try:
        control_task = get_control_task_name(args.dataset)
        print(f"Control task: {control_task}")

        drug_tasks = get_drug_task_names(args.dataset)
        print(f"Drug tasks: {drug_tasks}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Discover processed scans
    task_to_check = args.task if args.task else control_task
    print(f"\nDiscovering processed scans for task: {task_to_check}")

    processed = discover_processed_scans(args.dataset, task_filter=task_to_check)
    print(f"Found {len(processed)} processed scans")

    if len(processed) > 0:
        print("\nFirst 5 scans:")
        for scan in processed[:5]:
            session_str = f"ses-{scan['session']}" if scan['session'] else "no-session"
            print(f"  sub-{scan['subject']} {session_str} task-{scan['task']} run-{scan['run']}")

        # Try loading group data
        print(f"\nLoading group matrices for task: {task_to_check}")
        try:
            all_mats, subjects = load_group_all_measures(args.dataset, task_to_check)
            print(f"Loaded {len(subjects)} subjects")
            print(f"Synergy shape: {all_mats['synergy'].shape}")
            print(f"Redundancy shape: {all_mats['redundancy'].shape}")
        except Exception as e:
            print(f"Error loading group data: {e}")
