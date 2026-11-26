#!/usr/bin/env python3
"""
fMRIPrep Configuration Loader
=============================
Loads dataset-specific configurations and generates fMRIPrep parameters.
Used by the dataset-agnostic fMRIPrep PowerShell script.

Usage:
    python fmriprep_config_loader.py --dataset berlin_ketamine --env local
    python fmriprep_config_loader.py --dataset berlin_ketamine --env local --list-subjects
    python fmriprep_config_loader.py --dataset berlin_ketamine --env local --output json
"""

import os
import sys
import yaml
import json
import argparse
import platform
from pathlib import Path
from typing import Dict, List, Optional


def get_project_dir() -> Path:
    """Get project directory from script location."""
    # Script is in scripts/preprocessing/, go up two levels
    return Path(__file__).parent.parent.parent


def auto_detect_environment() -> str:
    """
    Auto-detect environment based on OS platform.

    Returns:
    --------
    str: 'local' on Windows, 'hpc' on Linux
    """
    system = platform.system()
    if system == 'Windows':
        return 'local'
    elif system == 'Linux':
        return 'hpc'
    else:
        # Default to local for macOS or unknown
        return 'local'


def validate_path_for_environment(path: str, env: str) -> bool:
    """
    Check if a path looks appropriate for the given environment.

    Returns True if path looks valid, False if it looks wrong.
    """
    if not path:
        return False

    # Check for obvious mismatches
    is_linux_path = path.startswith('/') and not path.startswith('//')
    is_windows_path = len(path) > 1 and path[1] == ':' or path.startswith('//')

    if env == 'local' and is_linux_path and not is_windows_path:
        return False
    if env == 'hpc' and is_windows_path:
        return False

    return True


def load_configs(dataset_name: str, project_dir: Optional[Path] = None) -> Dict:
    """
    Load dataset and global configurations.

    Returns dict with all necessary fMRIPrep parameters.
    """
    if project_dir is None:
        project_dir = get_project_dir()

    # Load main registry
    registry_path = project_dir / 'configs' / 'dataset_configs.yaml'
    with open(registry_path, 'r') as f:
        registry = yaml.safe_load(f)

    # Check dataset exists
    if dataset_name not in registry['datasets']:
        available = ', '.join(registry['datasets'].keys())
        raise ValueError(f"Dataset '{dataset_name}' not found. Available: {available}")

    # Load dataset-specific config
    dataset_config_file = registry['datasets'][dataset_name]['config_file']
    dataset_config_path = project_dir / dataset_config_file

    with open(dataset_config_path, 'r') as f:
        dataset_config = yaml.safe_load(f)

    return {
        'registry': registry,
        'dataset': dataset_config
    }


def get_paths(configs: Dict, env: str = 'local') -> Dict:
    """
    Get resolved paths for the specified environment.

    Parameters:
    -----------
    configs : dict
        Loaded configurations
    env : str
        Environment: 'local' or 'hpc'

    Returns:
    --------
    dict with bids_root, fmriprep_dir, license_file, work_dir, log_dir
    """
    registry = configs['registry']
    dataset = configs['dataset']
    dataset_name = dataset.get('dataset_id', 'unknown')

    # Get global paths
    global_paths = registry.get('global_paths', {}).get(env, {})

    # Get dataset paths
    dataset_paths = dataset.get('paths', {})

    # Try environment-specific paths first
    if env in dataset_paths:
        bids_root = dataset_paths[env].get('bids_root', '')
        fmriprep_dir = dataset_paths[env].get('fmriprep_dir', '')
    else:
        # Fall back to legacy paths, but warn if they look wrong
        bids_root = dataset_paths.get('bids_root', '')
        fmriprep_dir = dataset_paths.get('fmriprep_dir', '')

        # Validate paths match expected environment
        if not validate_path_for_environment(bids_root, env):
            print(f"WARNING: Dataset '{dataset_name}' does not have '{env}' paths configured.", file=sys.stderr)
            print(f"         The legacy paths appear to be for a different environment.", file=sys.stderr)
            print(f"         Please add a '{env}' section to configs/datasets/{dataset_name}.yaml", file=sys.stderr)
            print(f"         Current bids_root: {bids_root}", file=sys.stderr)

    # Resolve templateflow_dir with fallback to sibling of work_dir
    templateflow_dir = global_paths.get('templateflow_dir', '')
    if not templateflow_dir:
        work_dir = global_paths.get('work_dir', '')
        if work_dir:
            templateflow_dir = str(Path(work_dir).parent / 'templateflow')

    return {
        'bids_root': bids_root,
        'fmriprep_dir': fmriprep_dir,
        'license_file': global_paths.get('license_file', ''),
        'work_dir': global_paths.get('work_dir', ''),
        'log_dir': global_paths.get('log_dir', ''),
        'project_dir': global_paths.get('project_dir', str(get_project_dir())),
        'templateflow_dir': templateflow_dir
    }


def get_fmriprep_params(configs: Dict, resource_preset: str = 'medium') -> Dict:
    """
    Get fMRIPrep parameters from config.

    Parameters:
    -----------
    configs : dict
        Loaded configurations
    resource_preset : str
        Resource preset: 'low', 'medium', 'high', 'max'

    Returns:
    --------
    dict with version, output_spaces, nprocs, mem_mb, etc.
    """
    fmriprep_config = configs['registry'].get('fmriprep', {})

    # Get resource settings
    presets = fmriprep_config.get('resource_presets', {})
    resources = presets.get(resource_preset, presets.get('medium', {}))

    # Get dataset-specific extra args (e.g., --skull-strip-t1w skip)
    dataset_extra_args = configs['dataset'].get('fmriprep_extra_args', '')

    return {
        'version': fmriprep_config.get('version', '23.2.1'),
        'output_spaces': fmriprep_config.get('output_spaces', 'MNI152NLin2009cAsym:res-2'),
        'skip_bids_validation': fmriprep_config.get('skip_bids_validation', True),
        'stop_on_first_crash': fmriprep_config.get('stop_on_first_crash', True),
        'nprocs': resources.get('nprocs', 8),
        'mem_mb': resources.get('mem_mb', 32768),
        'extra_args': dataset_extra_args
    }


def discover_subjects(configs: Dict, bids_root: str) -> List[str]:
    """
    Discover subjects from BIDS directory or config.

    Returns list of subject IDs (without 'sub-' prefix).
    """
    dataset = configs['dataset']
    subjects_config = dataset.get('subjects', {})

    # Check if we should auto-discover
    if subjects_config.get('auto_discover', False):
        # Scan BIDS directory for subjects
        bids_path = Path(bids_root)
        if bids_path.exists():
            subjects = []
            for sub_dir in sorted(bids_path.glob('sub-*')):
                if sub_dir.is_dir():
                    sub_id = sub_dir.name.replace('sub-', '')
                    subjects.append(sub_id)
        else:
            subjects = []
    else:
        # Use explicitly listed subjects
        subjects = subjects_config.get('include', [])

    # Apply exclusions
    exclude = subjects_config.get('exclude', [])
    subjects = [s for s in subjects if s not in exclude]

    return sorted(subjects)


def check_processed_subjects(fmriprep_dir: str, subjects: List[str]) -> Dict[str, List[str]]:
    """
    Check which subjects have already been processed.

    A subject is considered fully processed only if it has:
    1. An HTML report (sub-XXXXXX.html)
    2. At least one preprocessed BOLD file (*desc-preproc_bold.nii.gz)
       OR at least one confounds file (*desc-confounds_timeseries.tsv)

    Returns dict with 'processed' and 'pending' lists.
    """
    fmriprep_path = Path(fmriprep_dir)

    processed = []
    pending = []

    for sub_id in subjects:
        # Check for HTML report
        report_path = fmriprep_path / f'sub-{sub_id}.html'
        if not report_path.exists():
            pending.append(sub_id)
            continue

        # Check for actual functional outputs (not just anat)
        # Use rglob to handle both session and no-session layouts
        sub_dir = fmriprep_path / f'sub-{sub_id}'
        has_bold = False
        if sub_dir.exists():
            preproc_bolds = list(sub_dir.rglob('*desc-preproc_bold.nii.gz'))
            confounds = list(sub_dir.rglob('*desc-confounds_timeseries.tsv'))
            has_bold = len(preproc_bolds) > 0 or len(confounds) > 0

        if has_bold:
            processed.append(sub_id)
        else:
            pending.append(sub_id)

    return {
        'processed': processed,
        'pending': pending
    }


def main():
    parser = argparse.ArgumentParser(
        description='Load fMRIPrep configuration for a dataset'
    )
    parser.add_argument(
        '--dataset', '-d',
        required=True,
        help='Dataset name (e.g., berlin_ketamine, london_lsd)'
    )
    parser.add_argument(
        '--env', '-e',
        choices=['auto', 'local', 'hpc'],
        default='auto',
        help='Environment (auto, local, or hpc). Auto detects based on OS.'
    )
    parser.add_argument(
        '--resources', '-r',
        choices=['low', 'medium', 'high', 'max'],
        default='medium',
        help='Resource preset'
    )
    parser.add_argument(
        '--list-subjects',
        action='store_true',
        help='List all subjects'
    )
    parser.add_argument(
        '--check-processed',
        action='store_true',
        help='Check which subjects are already processed'
    )
    parser.add_argument(
        '--output', '-o',
        choices=['text', 'json', 'powershell'],
        default='text',
        help='Output format'
    )
    parser.add_argument(
        '--project-dir',
        type=str,
        default=None,
        help='Override project directory'
    )

    args = parser.parse_args()

    # Resolve environment
    env = args.env
    if env == 'auto':
        env = auto_detect_environment()

    # Load configs
    project_dir = Path(args.project_dir) if args.project_dir else None
    configs = load_configs(args.dataset, project_dir)

    # Get paths and params
    paths = get_paths(configs, env)
    fmriprep_params = get_fmriprep_params(configs, args.resources)

    # Discover subjects
    subjects = discover_subjects(configs, paths['bids_root'])

    # Build output
    output = {
        'dataset': args.dataset,
        'dataset_name': configs['dataset'].get('name', args.dataset),
        'environment': env,
        'resource_preset': args.resources,
        'paths': paths,
        'fmriprep': fmriprep_params,
        'subjects': subjects,
        'subject_count': len(subjects)
    }

    # Check processed if requested
    if args.check_processed:
        status = check_processed_subjects(paths['fmriprep_dir'], subjects)
        output['processed'] = status['processed']
        output['pending'] = status['pending']
        output['processed_count'] = len(status['processed'])
        output['pending_count'] = len(status['pending'])

    # Format output
    if args.output == 'json':
        print(json.dumps(output, indent=2))

    elif args.output == 'powershell':
        # Output PowerShell variable assignments
        print(f"$DATASET = '{args.dataset}'")
        print(f"$DATASET_NAME = '{output['dataset_name']}'")
        print(f"$BIDS_ROOT = '{paths['bids_root']}'")
        print(f"$OUTPUT_DIR = '{paths['fmriprep_dir']}'")
        print(f"$LICENSE_FILE = '{paths['license_file']}'")
        print(f"$WORK_DIR = '{paths['work_dir']}'")
        print(f"$LOG_DIR = '{paths['log_dir']}'")
        print(f"$TEMPLATEFLOW_DIR = '{paths['templateflow_dir']}'")
        print(f"$FMRIPREP_VERSION = '{fmriprep_params['version']}'")
        print(f"$OUTPUT_SPACES = '{fmriprep_params['output_spaces']}'")
        print(f"$NPROCS = {fmriprep_params['nprocs']}")
        print(f"$MEM_MB = {fmriprep_params['mem_mb']}")
        print(f"$EXTRA_ARGS = '{fmriprep_params['extra_args']}'")
        print(f"$SUBJECTS = '{' '.join(subjects)}'")
        if args.check_processed:
            print(f"$PROCESSED = '{' '.join(output['processed'])}'")
            print(f"$PENDING = '{' '.join(output['pending'])}'")

    else:  # text
        print(f"Dataset: {output['dataset_name']} ({args.dataset})")
        print(f"Environment: {env}" + (f" (auto-detected)" if args.env == 'auto' else ""))
        print(f"Resource Preset: {args.resources}")
        print()
        print("Paths:")
        print(f"  BIDS Root: {paths['bids_root']}")
        print(f"  fMRIPrep Dir: {paths['fmriprep_dir']}")
        print(f"  License: {paths['license_file']}")
        print(f"  Work Dir: {paths['work_dir']}")
        print(f"  Log Dir: {paths['log_dir']}")
        print(f"  TemplateFlow: {paths['templateflow_dir']}")
        print()
        print("fMRIPrep Settings:")
        print(f"  Version: {fmriprep_params['version']}")
        print(f"  Output Spaces: {fmriprep_params['output_spaces']}")
        print(f"  Processors: {fmriprep_params['nprocs']}")
        print(f"  Memory: {fmriprep_params['mem_mb']} MB ({fmriprep_params['mem_mb'] // 1024} GB)")
        if fmriprep_params.get('extra_args'):
            print(f"  Extra Args: {fmriprep_params['extra_args']}")
        print()
        print(f"Subjects ({len(subjects)}):")
        if args.list_subjects:
            for sub in subjects:
                print(f"  - sub-{sub}")
        else:
            print(f"  {', '.join(subjects[:10])}", end='')
            if len(subjects) > 10:
                print(f" ... and {len(subjects) - 10} more")
            else:
                print()

        if args.check_processed:
            print()
            print(f"Processed: {output['processed_count']}")
            print(f"Pending: {output['pending_count']}")
            if output['pending']:
                print(f"  Pending subjects: {', '.join(output['pending'])}")


if __name__ == '__main__':
    main()
