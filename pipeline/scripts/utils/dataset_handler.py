#!/usr/bin/env python3
"""
Dataset Handler - Universal data discovery and path resolution
Handles multiple dataset structures transparently

This utility abstracts away dataset-specific details and provides a uniform
interface for accessing fMRIPrep data across different BIDS structures.

Also integrates with QC exclusion system to automatically filter scans
based on motion quality criteria.
"""

import os
import re
import yaml
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from bids import BIDSLayout
from bids.layout import BIDSLayoutIndexer
import warnings

# Import QC loader if available
try:
    from qc_loader import QCExclusions, load_qc_exclusions
    HAS_QC_LOADER = True
except ImportError:
    HAS_QC_LOADER = False


class DatasetHandler:
    """
    Handles dataset-specific operations including data discovery and path resolution
    """

    def __init__(self, dataset_name: str, project_dir: Optional[str] = None):
        """
        Initialize dataset handler

        Parameters:
        -----------
        dataset_name : str
            Name of dataset (e.g., 'london_lsd', 'basel_mmm')
        project_dir : str, optional
            Path to project directory. If None, infers from script location
        """
        self.dataset_name = dataset_name

        # Find project directory
        if project_dir is None:
            # Assume we're in scripts/utils/ and go up two levels
            project_dir = Path(__file__).parent.parent.parent
        self.project_dir = Path(project_dir)

        # Load configurations
        self._load_configs()

        # Initialize BIDS layout (lazy loading)
        self._layout = None

    def _load_configs(self):
        """Load dataset and global configurations"""
        # Load main dataset registry
        registry_path = self.project_dir / 'configs' / 'dataset_configs.yaml'
        with open(registry_path, 'r') as f:
            registry = yaml.safe_load(f)

        # Store global configs
        self.global_preprocessing = registry.get('global_preprocessing', {})
        self.global_analysis = registry.get('global_analysis', {})

        # Get dataset config file path
        if self.dataset_name not in registry['datasets']:
            available = ', '.join(registry['datasets'].keys())
            raise ValueError(
                f"Dataset '{self.dataset_name}' not found. "
                f"Available datasets: {available}"
            )

        dataset_config_file = registry['datasets'][self.dataset_name]['config_file']
        dataset_config_path = self.project_dir / dataset_config_file

        # Load dataset-specific config
        with open(dataset_config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        # Merge preprocessing configs (dataset-specific overrides global)
        self.preprocessing = {**self.global_preprocessing}
        if 'preprocessing' in self.config:
            self.preprocessing.update(self.config['preprocessing'])

        print(f"Loaded dataset config: {self.config['name']}")

    def _get_env_paths(self) -> dict:
        """Get paths for current environment (local vs hpc)"""
        paths = self.config.get('paths', {})

        # Check if paths use local/hpc structure
        if 'local' in paths or 'hpc' in paths:
            # Auto-detect environment based on OS
            import platform
            if platform.system() == 'Windows':
                env = 'local'
            else:
                # Check if we're on HPC by looking for typical HPC paths
                env = 'hpc' if Path('/home/s3648540').exists() else 'local'

            return paths.get(env, paths.get('local', {}))
        else:
            # Legacy format - paths directly in config
            return paths

    @property
    def fmriprep_dir(self) -> str:
        """Get fMRIPrep directory for current environment"""
        env_paths = self._get_env_paths()
        return env_paths.get('fmriprep_dir', self.config['paths'].get('fmriprep_dir', ''))

    @property
    def bids_root(self) -> str:
        """Get BIDS root directory for current environment"""
        env_paths = self._get_env_paths()
        return env_paths.get('bids_root', self.config['paths'].get('bids_root', ''))

    @property
    def layout(self) -> BIDSLayout:
        """Lazy-load BIDS layout for fMRIPrep derivatives"""
        if self._layout is None:
            fmriprep_dir = self.fmriprep_dir
            print(f"Initializing BIDS layout for: {fmriprep_dir}")
            # Don't use derivatives=True since we're pointing directly at the derivative folder
            # derivatives=True expects a raw BIDS root and looks for /derivatives inside
            # Use indexer to skip directories with symlinks that break on Windows
            # Pass pre-compiled regex patterns to avoid PyBIDS _regexfy bug on Windows
            # (_regexfy generates forward-slash patterns that don't match Windows backslash paths)
            indexer = BIDSLayoutIndexer(ignore=[
                re.compile(r'sourcedata'),
                re.compile(r'freesurfer'),
                re.compile(r'\.git'),
            ])
            self._layout = BIDSLayout(
                fmriprep_dir,
                validate=False,
                derivatives=False,
                indexer=indexer
            )
        return self._layout

    def get_subjects(self) -> List[str]:
        """
        Get list of subjects to process

        Returns:
        --------
        list of str
            Subject IDs (without 'sub-' prefix)
        """
        subjects_config = self.config['subjects']

        # Auto-discover from BIDS if enabled
        if subjects_config.get('auto_discover', False):
            subjects = self.layout.get_subjects()
        else:
            # Use explicitly listed subjects
            subjects = subjects_config.get('include', [])

        # Apply exclusions
        exclude = subjects_config.get('exclude', [])
        subjects = [s for s in subjects if s not in exclude]

        print(f"Found {len(subjects)} subjects")
        return sorted(subjects)

    def get_sessions(self, subject: str) -> List[Optional[str]]:
        """
        Get list of sessions for a subject

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)

        Returns:
        --------
        list
            Session IDs (without 'ses-' prefix), or [None] if no sessions
        """
        has_sessions = self.config['bids_structure']['has_sessions']

        if not has_sessions:
            return [None]

        # Get sessions from BIDS
        sessions = self.layout.get_sessions(subject=subject)
        return sorted(sessions) if sessions else [None]

    def get_tasks(self, subject: str, session: Optional[str] = None) -> List[str]:
        """
        Get list of tasks for a subject/session

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        session : str, optional
            Session ID (without 'ses-' prefix)

        Returns:
        --------
        list of str
            Task names (without 'task-' prefix)
        """
        # Get tasks from BIDS
        kwargs = {'subject': subject, 'suffix': 'bold'}
        if session is not None:
            kwargs['session'] = session

        files = self.layout.get(**kwargs)

        # Extract unique tasks
        tasks = set()
        for f in files:
            entities = self.layout.parse_file_entities(f)
            if 'task' in entities:
                tasks.add(entities['task'])

        # Filter to configured tasks (using bids_task_name from config)
        valid_bids_tasks = set()
        for task_key, task_info in self.config['tasks'].items():
            if isinstance(task_info, dict) and 'bids_task_name' in task_info:
                valid_bids_tasks.add(task_info['bids_task_name'])
            else:
                valid_bids_tasks.add(task_key)  # Fallback to key
        tasks = [t for t in tasks if t in valid_bids_tasks]

        return sorted(tasks)

    def get_runs(self, subject: str, task: str, session: Optional[str] = None) -> List[int]:
        """
        Get list of runs for a subject/task/session

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        task : str
            Task name (without 'task-' prefix)
        session : str, optional
            Session ID (without 'ses-' prefix)

        Returns:
        --------
        list of int
            Run numbers (returns [1] if files exist but have no run entity)
        """
        # Get runs from BIDS
        kwargs = {'subject': subject, 'task': task, 'suffix': 'bold'}
        if session is not None:
            kwargs['session'] = session

        files = self.layout.get(**kwargs)

        # Extract run numbers
        runs = set()
        for f in files:
            entities = self.layout.parse_file_entities(f)
            if 'run' in entities:
                runs.add(entities['run'])

        # If no runs found but files exist, treat as implicit run 1
        if not runs and len(files) > 0:
            return [1]

        return sorted(runs)

    def select_run(self, subject: str, task: str, session: Optional[str] = None) -> Optional[int]:
        """
        Select run based on dataset's run selection strategy

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        task : str
            Task name (without 'task-' prefix)
        session : str, optional
            Session ID (without 'ses-' prefix)

        Returns:
        --------
        int or None
            Selected run number, or None if no runs found
        """
        runs = self.get_runs(subject, task, session)

        if not runs:
            return None

        strategy = self.config['run_selection']['strategy']

        if strategy == 'all':
            # Return all runs (caller will handle multiple)
            return runs
        elif strategy == 'last':
            return max(runs)
        elif strategy == 'first':
            return min(runs)
        elif strategy == 'specific':
            specific_run = self.config['run_selection']['specific_run']
            if specific_run in runs:
                return specific_run
            else:
                warnings.warn(
                    f"Specific run {specific_run} not found for "
                    f"sub-{subject} task-{task}. Using last available run."
                )
                return max(runs)
        else:
            raise ValueError(f"Unknown run selection strategy: {strategy}")

    def get_functional_file(
        self,
        subject: str,
        task: str,
        run: int,
        session: Optional[str] = None,
        space: str = 'MNI152NLin2009cAsym',
        desc: str = 'preproc'
    ) -> Optional[str]:
        """
        Get path to functional (BOLD) file

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        task : str
            Task name (without 'task-' prefix)
        run : int
            Run number (if 1 and no run entity exists, will search without run)
        session : str, optional
            Session ID (without 'ses-' prefix)
        space : str
            Space (default: MNI152NLin2009cAsym)
        desc : str
            Description (default: preproc)

        Returns:
        --------
        str or None
            Path to BOLD file, or None if not found
        """
        kwargs = {
            'subject': subject,
            'task': task,
            'run': run,
            'suffix': 'bold',
            'space': space,
            'desc': desc,
            'extension': 'nii.gz',
            'return_type': 'file'
        }
        if session is not None:
            kwargs['session'] = session

        files = self.layout.get(**kwargs)

        # If no files found and run=1, try without run entity (implicit run)
        if not files and run == 1:
            kwargs_no_run = kwargs.copy()
            del kwargs_no_run['run']
            files = self.layout.get(**kwargs_no_run)

        if not files:
            return None

        return files[0]  # Should be unique

    def get_confounds_file(
        self,
        subject: str,
        task: str,
        run: int,
        session: Optional[str] = None
    ) -> Optional[str]:
        """
        Get path to confounds file

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        task : str
            Task name (without 'task-' prefix)
        run : int
            Run number (if 1 and no run entity exists, will search without run)
        session : str, optional
            Session ID (without 'ses-' prefix)

        Returns:
        --------
        str or None
            Path to confounds file, or None if not found
        """
        kwargs = {
            'subject': subject,
            'task': task,
            'run': run,
            'desc': 'confounds',
            'suffix': 'timeseries',
            'extension': 'tsv',
            'return_type': 'file'
        }
        if session is not None:
            kwargs['session'] = session

        files = self.layout.get(**kwargs)

        # If no files found and run=1, try without run entity (implicit run)
        if not files and run == 1:
            kwargs_no_run = kwargs.copy()
            del kwargs_no_run['run']
            files = self.layout.get(**kwargs_no_run)

        if not files:
            return None

        return files[0]

    def get_output_filename(
        self,
        subject: str,
        task: str,
        run: int,
        session: Optional[str] = None,
        suffix: str = 'timeseries',
        extension: str = 'npy',
        include_dataset: bool = True
    ) -> str:
        """
        Generate standardized output filename with dataset prefix

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        task : str
            Task name (without 'task-' prefix)
        run : int
            Run number
        session : str, optional
            Session ID (without 'ses-' prefix)
        suffix : str
            Suffix (e.g., 'timeseries', 'phiid', 'synergy')
        extension : str
            File extension (without dot)
        include_dataset : bool
            Whether to include dataset name as prefix (default: True)

        Returns:
        --------
        str
            Standardized filename with dataset prefix
            Format: dataset-{dataset}_sub-{subject}_[ses-{session}_]task-{task}_run-{run}_{suffix}.{extension}
        """
        parts = []

        # Add dataset prefix
        if include_dataset:
            parts.append(f'dataset-{self.dataset_name}')

        # Add subject
        parts.append(f'sub-{subject}')

        # Add session if present
        if session is not None:
            parts.append(f'ses-{session}')

        # Add task, run, and suffix
        parts.extend([
            f'task-{task}',
            f'run-{run}',
            suffix
        ])

        filename = '_'.join(parts) + f'.{extension}'
        return filename

    def discover_all_scans(self) -> List[Dict]:
        """
        Discover all available scans in the dataset

        Returns:
        --------
        list of dict
            Each dict contains: subject, session, task, run, func_file, confounds_file
        """
        scans = []

        subjects = self.get_subjects()

        for subject in subjects:
            sessions = self.get_sessions(subject)

            for session in sessions:
                tasks = self.get_tasks(subject, session)

                for task in tasks:
                    # Get selected run(s)
                    selected_run = self.select_run(subject, task, session)

                    if selected_run is None:
                        continue

                    # Handle 'all' strategy returning multiple runs
                    if isinstance(selected_run, list):
                        runs_to_process = selected_run
                    else:
                        runs_to_process = [selected_run]

                    for run in runs_to_process:
                        func_file = self.get_functional_file(subject, task, run, session)
                        confounds_file = self.get_confounds_file(subject, task, run, session)

                        if func_file and confounds_file:
                            # Normalise run to plain int. Pybids returns PaddedInt
                            # for zero-padded BIDS entities (e.g. cambridge_propofol's
                            # run-01) which stringifies as "01"; downstream filename
                            # builders use str(run) and so produce run-01 names,
                            # mismatching the run-1 convention used elsewhere.
                            scan_info = {
                                'subject': subject,
                                'session': session,
                                'task': task,
                                'run': int(run),
                                'func_file': func_file,
                                'confounds_file': confounds_file
                            }
                            scans.append(scan_info)

        return scans

    def get_task_info(self, task: str) -> Dict:
        """
        Get information about a task

        Parameters:
        -----------
        task : str
            Task name (without 'task-' prefix)

        Returns:
        --------
        dict
            Task information including description, is_control flag
        """
        return self.config['tasks'].get(task, {})

    # =========================================================================
    # QC Exclusion Integration
    # =========================================================================

    def load_qc_exclusions(self) -> Optional['QCExclusions']:
        """
        Load QC exclusions for this dataset.

        Returns:
        --------
        QCExclusions or None
            QC exclusions object, or None if not available
        """
        if not HAS_QC_LOADER:
            warnings.warn(
                "QC loader not available. Install qc_loader.py in scripts/utils/."
            )
            return None

        try:
            return load_qc_exclusions(self.dataset_name, self.project_dir)
        except Exception as e:
            warnings.warn(f"Could not load QC exclusions: {e}")
            return None

    def is_scan_excluded(
        self,
        subject: str,
        session: Optional[str] = None,
        task: Optional[str] = None,
        qc: Optional['QCExclusions'] = None
    ) -> bool:
        """
        Check if a scan should be excluded based on QC criteria.

        Parameters:
        -----------
        subject : str
            Subject ID (without 'sub-' prefix)
        session : str, optional
            Session ID (without 'ses-' prefix)
        task : str, optional
            BIDS task name. Required for datasets that encode conditions in the
            task name (has_sessions=false); if omitted, QC collapses to worst-case.
        qc : QCExclusions, optional
            Pre-loaded QC exclusions. If None, loads automatically.

        Returns:
        --------
        bool
            True if scan should be excluded
        """
        if qc is None:
            qc = self.load_qc_exclusions()

        if qc is None:
            # No QC data available, include by default
            return False

        return qc.is_excluded(subject, session, task)

    def get_scan_qc_status(
        self,
        subject: str,
        session: Optional[str] = None,
        task: Optional[str] = None,
        qc: Optional['QCExclusions'] = None
    ) -> str:
        """
        Get QC status for a scan.

        Returns: 'INCLUDE', 'EXCLUDE', 'TRUNCATE', or 'UNKNOWN'
        """
        if qc is None:
            qc = self.load_qc_exclusions()

        if qc is None:
            return 'UNKNOWN'

        return qc.get_status(subject, session, task)

    def get_truncation_info(
        self,
        subject: str,
        session: Optional[str] = None,
        task: Optional[str] = None,
        qc: Optional['QCExclusions'] = None
    ) -> Optional[Dict]:
        """
        Get truncation info for a scan if it's salvageable.

        Returns dict with start_volume, end_volume, n_volumes, duration_minutes
        or None if not applicable.
        """
        if qc is None:
            qc = self.load_qc_exclusions()

        if qc is None:
            return None

        truncation = qc.get_truncation(subject, session, task)
        if truncation is None:
            return None

        return {
            'start_volume': truncation.start_volume,
            'end_volume': truncation.end_volume,
            'n_volumes': truncation.n_volumes,
            'duration_minutes': truncation.duration_minutes,
            'reason': truncation.reason
        }

    def get_included_scans(self, respect_qc: bool = True) -> List[Dict]:
        """
        Discover all scans that pass QC criteria.

        This is the main method for getting scans to process.
        It combines discover_all_scans() with QC filtering.

        Parameters:
        -----------
        respect_qc : bool
            If True, filter out QC-excluded scans (default: True)

        Returns:
        --------
        list of dict
            Each dict contains: subject, session, task, run, func_file,
            confounds_file, qc_status, truncation_info (if applicable)
        """
        # Get all available scans
        all_scans = self.discover_all_scans()

        if not respect_qc:
            return all_scans

        # Load QC exclusions
        qc = self.load_qc_exclusions()

        if qc is None:
            print("Note: No QC exclusions loaded. Including all scans.")
            return all_scans

        # Filter scans based on QC
        included_scans = []
        excluded_count = 0
        truncate_count = 0

        for scan in all_scans:
            subject = scan['subject']
            session = scan['session']
            task = scan.get('task')

            status = qc.get_status(subject, session, task)
            scan['qc_status'] = status

            if qc.is_excluded(subject, session, task):
                excluded_count += 1
                continue

            # Add truncation info if applicable
            if status == 'TRUNCATE':
                truncation = qc.get_truncation(subject, session, task)
                if truncation:
                    scan['truncation_info'] = {
                        'start_volume': truncation.start_volume,
                        'end_volume': truncation.end_volume,
                        'n_volumes': truncation.n_volumes
                    }
                    truncate_count += 1

            included_scans.append(scan)

        print(f"QC filtering: {len(included_scans)} included, "
              f"{excluded_count} excluded, {truncate_count} to truncate")

        return included_scans

    def get_qc_summary(self) -> Dict:
        """
        Get summary of QC exclusions for this dataset.

        Returns:
        --------
        dict
            Summary statistics including counts and rates
        """
        qc = self.load_qc_exclusions()

        if qc is None:
            return {'error': 'No QC data available'}

        return qc.summary()

    def __repr__(self):
        """String representation"""
        return (
            f"DatasetHandler(dataset='{self.dataset_name}', "
            f"subjects={len(self.get_subjects())}, "
            f"has_sessions={self.config['bids_structure']['has_sessions']})"
        )


def load_dataset_handler(dataset_name: str, project_dir: Optional[str] = None) -> DatasetHandler:
    """
    Convenience function to load a dataset handler

    Parameters:
    -----------
    dataset_name : str
        Name of dataset (e.g., 'london_lsd', 'basel_mmm')
    project_dir : str, optional
        Path to project directory

    Returns:
    --------
    DatasetHandler
        Initialized dataset handler
    """
    return DatasetHandler(dataset_name, project_dir)


# Datasets explicitly excluded from the "active" set used by --all flags in the
# CONN preprocessing pipeline. Update this set if scope changes.
# missouri_psilo activated 2026-08-26 (rollout completed through the uniform
# CONN chain; see data/missouri_psilo/PREPROCESSING_NOTES.md).
_INACTIVE_DATASETS = frozenset({'basel_mmm'})


def list_active_datasets(project_dir: Optional[str] = None) -> List[str]:
    """
    List dataset names registered in configs/dataset_configs.yaml, excluding
    those marked inactive in `_INACTIVE_DATASETS`.

    Used by `--all` flags in the CONN preprocessing scripts to iterate every
    dataset in the project's current rollout scope (all registered datasets
    except Basel MMM, which is not part of this study).

    Parameters
    ----------
    project_dir : str, optional
        Project directory. Defaults to two levels up from this file.

    Returns
    -------
    list of str
        Sorted list of active dataset names.
    """
    if project_dir is None:
        project_dir = Path(__file__).parent.parent.parent
    else:
        project_dir = Path(project_dir)

    registry_path = project_dir / 'configs' / 'dataset_configs.yaml'
    with open(registry_path, 'r') as f:
        registry = yaml.safe_load(f)

    all_datasets = set(registry.get('datasets', {}).keys())
    active = sorted(all_datasets - _INACTIVE_DATASETS)
    return active


# CLI for testing/debugging
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Test dataset handler')
    parser.add_argument('--dataset', required=True, help='Dataset name')
    parser.add_argument('--discover', action='store_true', help='Discover all scans')
    parser.add_argument('--included', action='store_true', help='Show QC-included scans only')
    parser.add_argument('--qc-summary', action='store_true', help='Show QC summary')
    parser.add_argument('--subject', help='Test specific subject')
    args = parser.parse_args()

    # Load handler
    handler = load_dataset_handler(args.dataset)
    print(handler)
    print()

    if args.qc_summary:
        print("QC Summary:")
        print("-" * 40)
        summary = handler.get_qc_summary()
        if 'error' in summary:
            print(f"  {summary['error']}")
        else:
            for key, value in summary.items():
                print(f"  {key}: {value}")
        print()

    if args.discover:
        print("Discovering all scans...")
        scans = handler.discover_all_scans()
        print(f"\nFound {len(scans)} scans:")
        for i, scan in enumerate(scans[:10], 1):  # Show first 10
            session_str = f"ses-{scan['session']}" if scan['session'] else "no-session"
            print(f"  {i}. sub-{scan['subject']} {session_str} "
                  f"task-{scan['task']} run-{scan['run']}")
        if len(scans) > 10:
            print(f"  ... and {len(scans) - 10} more")

    if args.included:
        print("Discovering QC-included scans...")
        scans = handler.get_included_scans()
        print(f"\nFound {len(scans)} included scans:")
        for i, scan in enumerate(scans[:15], 1):  # Show first 15
            session_str = f"ses-{scan['session']}" if scan['session'] else "no-session"
            qc_status = scan.get('qc_status', 'UNKNOWN')
            trunc_info = ""
            if 'truncation_info' in scan:
                t = scan['truncation_info']
                trunc_info = f" [TRUNCATE: vol {t['start_volume']}-{t['end_volume']}]"
            print(f"  {i}. sub-{scan['subject']} {session_str} "
                  f"task-{scan['task']} run-{scan['run']} ({qc_status}){trunc_info}")
        if len(scans) > 15:
            print(f"  ... and {len(scans) - 15} more")

    if args.subject:
        print(f"\nTesting subject: {args.subject}")
        sessions = handler.get_sessions(args.subject)
        print(f"  Sessions: {sessions}")

        for session in sessions:
            # Check QC status
            qc_status = handler.get_scan_qc_status(args.subject, session)
            excluded = handler.is_scan_excluded(args.subject, session)
            truncation = handler.get_truncation_info(args.subject, session)

            status_str = f"QC: {qc_status}"
            if excluded:
                status_str += " [EXCLUDED]"
            if truncation:
                status_str += f" [Truncate vol {truncation['start_volume']}-{truncation['end_volume']}]"

            print(f"  Session {session}: {status_str}")

            tasks = handler.get_tasks(args.subject, session)
            print(f"    Tasks: {tasks}")

            for task in tasks:
                runs = handler.get_runs(args.subject, task, session)
                selected_run = handler.select_run(args.subject, task, session)
                print(f"      Task {task}: runs={runs}, selected={selected_run}")
