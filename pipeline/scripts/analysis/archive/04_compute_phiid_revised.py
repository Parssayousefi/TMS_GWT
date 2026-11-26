#!/usr/bin/env python3
"""
Compute Integrated Information Decomposition (PhiID) measures
Using the Luppi et al. implementation from:
https://github.com/Imperial-MIND-lab/integrated-info-decomp

Supports two modes:
1. Single scan mode: --dataset X --subject Y --task Z --run N
2. Dataset batch mode: --dataset X [--subject Y] [--session S] [--task T]

UPDATED: Correct formulas per Luppi et al. (2024) Eqs 4-7
"""

import os
import sys
import yaml
import numpy as np
import argparse
from pathlib import Path
from typing import List, Dict, Optional
import time

# Add utils and tools to path
SCRIPT_DIR = Path(__file__).parent
UTILS_DIR = SCRIPT_DIR.parent / 'utils'
PROJECT_DIR = SCRIPT_DIR.parent.parent
TOOLS_DIR = PROJECT_DIR / 'tools' / 'integrated-info-decomp'
sys.path.insert(0, str(UTILS_DIR))
sys.path.insert(0, str(TOOLS_DIR))
from dataset_handler import load_dataset_handler

# Import phyid library
from phyid.calculate import calc_PhiID
from phyid.utils import PhiID_atoms_abbr

def compute_phiid_luppi_method(timeseries, tau=1, method='gaussian', redundancy='MMI'):
    """
    Compute PhiID using Luppi et al.'s exact implementation via phyid library

    Parameters:
    -----------
    timeseries : np.ndarray
        Shape (n_timepoints, n_rois)
    tau : int
        Time delay in TRs (default=1 as per paper)
    method : str
        'gaussian' for Gaussian MMI-PID (as used in the paper)
    redundancy : str
        'MMI' for minimum mutual information (default, as per paper)

    Returns:
    --------
    dict with keys: redundancy, synergy, phi, phi_r matrices

    Notes:
    ------
    Implements Luppi et al. (2024) using PhiID atoms from calc_PhiID().

    Per the Luppi supplementary code README (docs/luppi_supps/software/
    MMI_PhiID/README.md):
        "the atom that corresponds to (persistent) synergy is termed sts,
         and the atom that corresponds to (persistent) redundancy is
         termed rtr."

    Measures computed:
        Red(X,Y)   = atoms['rtr']  (persistent redundancy / double-redundancy)
        Syn(X,Y)   = atoms['sts']  (persistent synergy / double-synergy)
        Phi(X,Y)   = TDMI - I(X_past;X_fut) - I(Y_past;Y_fut)  (Eq. 6)
        Phi_R(X,Y) = Phi + Red                                   (Eq. 7)
    """

    n_timepoints, n_rois = timeseries.shape

    # Initialize output matrices
    redundancy_matrix = np.zeros((n_rois, n_rois))
    synergy_matrix = np.zeros((n_rois, n_rois))
    synergy_total_matrix = np.zeros((n_rois, n_rois))  # Eq. 5 total synergy (for Luppi comparison)
    phi_matrix = np.zeros((n_rois, n_rois))
    phi_r_matrix = np.zeros((n_rois, n_rois))

    # Process each ROI pair
    for i in range(n_rois):
        if i % 50 == 0:
            print(f"Processing ROI {i}/{n_rois}")

        for j in range(i+1, n_rois):
            # Extract timeseries for this pair
            src_ts = timeseries[:, i]
            trg_ts = timeseries[:, j]

            try:
                # Compute PhiID using phyid library
                atoms_res, calc_res = calc_PhiID(src_ts, trg_ts, tau,
                                                 kind=method, redundancy=redundancy)

                # Use intermediate MI values from calc_res (NOT raw atoms)
                I_res = calc_res['I_res']

                # Persistent redundancy: rtr atom (double-redundancy)
                # Filter non-finite values before averaging (as in Luppi MATLAB code)
                rtr_vals = atoms_res['rtr']
                rtr = np.mean(rtr_vals[np.isfinite(rtr_vals)]) if np.any(np.isfinite(rtr_vals)) else 0.0

                # Persistent synergy: sts atom (double-synergy)
                # Per Luppi supplementary code: synergy_mat(row,col) = atoms.sts
                sts_vals = atoms_res['sts']
                syn = np.mean(sts_vals[np.isfinite(sts_vals)]) if np.any(np.isfinite(sts_vals)) else 0.0

                # Helper to safely average with isfinite guard
                def _finite_mean(arr):
                    finite = arr[np.isfinite(arr)]
                    return np.mean(finite) if len(finite) > 0 else 0.0

                # Total synergy (all source-synergy atoms: str + stx + sty + sts)
                # Matches MATLAB pipeline: atoms.str + atoms.stx + atoms.sty + atoms.sts
                syn_total = (_finite_mean(atoms_res['str'])
                             + _finite_mean(atoms_res['stx'])
                             + _finite_mean(atoms_res['sty'])
                             + syn)  # syn = finite_mean of sts, already computed above

                # Eq. 6: Phi = TDMI - I(X_past; X_fut) - I(Y_past; Y_fut)
                phi = (_finite_mean(I_res['I_xytab'])
                       - _finite_mean(I_res['I_xta'])
                       - _finite_mean(I_res['I_ytb']))

                # Eq. 7: Phi_R = Phi + Red(X,Y)
                phi_r = phi + rtr

                # Store results (symmetric matrix)
                redundancy_matrix[i, j] = redundancy_matrix[j, i] = rtr
                synergy_matrix[i, j] = synergy_matrix[j, i] = syn
                synergy_total_matrix[i, j] = synergy_total_matrix[j, i] = syn_total
                phi_matrix[i, j] = phi_matrix[j, i] = phi
                phi_r_matrix[i, j] = phi_r_matrix[j, i] = phi_r

            except Exception as e:
                print(f"Warning: Failed to compute PhiID for ROI pair ({i}, {j}): {e}")
                # Leave as zeros
                continue

    return {
        'redundancy': redundancy_matrix,
        'synergy': synergy_matrix,
        'synergy_total': synergy_total_matrix,
        'phi': phi_matrix,
        'phi_r': phi_r_matrix
    }

def discover_deconv_files(
    timeseries_dir: Path,
    dataset: Optional[str] = None,
    subject: Optional[str] = None,
    session: Optional[str] = None,
    task: Optional[str] = None
) -> List[Path]:
    """
    Discover HRF-deconvolved timeseries files ready for PhiID computation.

    Returns list of *_timeseries_deconv.npy files matching the filters.
    """
    all_files = list(timeseries_dir.glob("*_timeseries_deconv.npy"))

    filtered = []
    for f in all_files:
        name = f.stem
        if dataset and f"dataset-{dataset}_" not in name:
            continue
        if subject and f"sub-{subject}_" not in name:
            continue
        if session and f"ses-{session}_" not in name:
            continue
        if task and f"task-{task}_" not in name:
            continue
        filtered.append(f)

    return sorted(filtered)


def parse_scan_info(filepath: Path) -> Dict:
    """Extract subject/session/task/run from a deconv timeseries filename."""
    import re
    name = filepath.stem  # e.g. dataset-london_lsd_sub-01_ses-LSD_task-rest_run-3_timeseries_deconv

    info = {}
    m = re.search(r'dataset-([^_]+(?:_[^_]+)*?)_sub-', name)
    if m:
        info['dataset'] = m.group(1)

    m = re.search(r'sub-([^_]+)', name)
    if m:
        info['subject'] = m.group(1)

    m = re.search(r'ses-([^_]+)', name)
    if m:
        info['session'] = m.group(1)

    m = re.search(r'task-([^_]+)', name)
    if m:
        info['task'] = m.group(1)

    m = re.search(r'run-(\d+)', name)
    if m:
        info['run'] = int(m.group(1))

    return info


def check_phiid_exists(output_dir: Path, scan_info: Dict, dataset: str) -> bool:
    """Check if PhiID output already exists for a scan (checks for phi_r.npy)."""
    # Build expected output filename pattern
    parts = [f"dataset-{dataset}"]
    parts.append(f"sub-{scan_info['subject']}")
    if scan_info.get('session'):
        parts.append(f"ses-{scan_info['session']}")
    if scan_info.get('task'):
        parts.append(f"task-{scan_info['task']}")
    if scan_info.get('run'):
        parts.append(f"run-{scan_info['run']}")
    parts.append("phi_r.npy")

    expected_file = output_dir / "_".join(parts)
    return expected_file.exists()


def process_single_scan(
    deconv_file: Path,
    output_dir: Path,
    dataset: str,
    tau: int = 1,
    method: str = 'gaussian',
    overwrite: bool = False,
    verbose: bool = True
) -> Dict:
    """
    Process a single deconvolved timeseries file through PhiID computation.

    Returns dict with status info.
    """
    result = {
        'input_path': str(deconv_file),
        'status': 'unknown',
        'message': ''
    }

    scan_info = parse_scan_info(deconv_file)

    # Check if already computed
    if not overwrite and check_phiid_exists(output_dir, scan_info, dataset):
        result['status'] = 'skipped'
        result['message'] = 'Already computed (use --overwrite to reprocess)'
        if verbose:
            print(f"  Skipping {deconv_file.name} - already computed")
        return result

    try:
        if verbose:
            print(f"  Loading: {deconv_file.name}")

        timeseries = np.load(deconv_file)
        if verbose:
            print(f"  Shape: {timeseries.shape}")

        # No z-scoring: calc_PhiID (like Luppi's PhiIDFull.m) normalizes
        # to unit variance internally. Matches MATLAB pipeline exactly.

        # Compute PhiID
        if verbose:
            print(f"  Computing PhiID (method={method}, tau={tau})...")

        start_time = time.time()
        phiid_results = compute_phiid_luppi_method(
            timeseries, tau=tau, method=method
        )
        elapsed = time.time() - start_time

        if verbose:
            print(f"  Completed in {elapsed/60:.1f} minutes")

        # Build output filename prefix
        parts = [f"dataset-{dataset}"]
        parts.append(f"sub-{scan_info['subject']}")
        if scan_info.get('session'):
            parts.append(f"ses-{scan_info['session']}")
        if scan_info.get('task'):
            parts.append(f"task-{scan_info['task']}")
        if scan_info.get('run'):
            parts.append(f"run-{scan_info['run']}")
        prefix = "_".join(parts)

        # Save measure matrices
        for measure, matrix in phiid_results.items():
            out_file = output_dir / f"{prefix}_{measure}.npy"
            np.save(out_file, matrix)
            if verbose:
                print(f"  Saved: {out_file.name}")

        # Save summary YAML
        summary = {
            'dataset': dataset,
            'subject': scan_info.get('subject'),
            'session': scan_info.get('session'),
            'task': scan_info.get('task'),
            'run': scan_info.get('run'),
            'tau': tau,
            'method': method,
            'n_timepoints': int(timeseries.shape[0]),
            'n_rois': int(timeseries.shape[1]),
            'elapsed_minutes': round(elapsed / 60, 2),
            'mean_redundancy': float(np.mean(phiid_results['redundancy'])),
            'mean_synergy': float(np.mean(phiid_results['synergy'])),
            'mean_synergy_total': float(np.mean(phiid_results['synergy_total'])),
            'mean_phi': float(np.mean(phiid_results['phi'])),
            'mean_phi_r': float(np.mean(phiid_results['phi_r'])),
            'negative_phi_percentage': float(np.mean(phiid_results['phi'] < 0) * 100),
            'positive_phi_r_check': bool(np.all(phiid_results['phi_r'] >= 0))
        }

        summary_file = output_dir / f"{prefix}_summary.yaml"
        with open(summary_file, 'w') as f:
            yaml.dump(summary, f)

        if verbose:
            print(f"  Mean Phi_R = {summary['mean_phi_r']:.4f}")
            print(f"  Phi_R non-negative check: {summary['positive_phi_r_check']}")

        result['status'] = 'success'
        result['message'] = f"Completed in {elapsed/60:.1f} min"
        result['summary'] = summary

    except Exception as e:
        result['status'] = 'error'
        result['message'] = str(e)
        if verbose:
            print(f"  ERROR: {e}")

    return result


def print_summary(results: List[Dict]):
    """Print summary of batch processing results."""
    success = sum(1 for r in results if r['status'] == 'success')
    skipped = sum(1 for r in results if r['status'] == 'skipped')
    errors = sum(1 for r in results if r['status'] == 'error')

    print("\n" + "=" * 60)
    print("PHIID COMPUTATION SUMMARY")
    print("=" * 60)
    print(f"  Successful: {success}")
    print(f"  Skipped (already done): {skipped}")
    print(f"  Errors: {errors}")
    print(f"  Total: {len(results)}")

    if success > 0:
        summaries = [r['summary'] for r in results if r.get('summary')]
        if summaries:
            total_time = sum(s['elapsed_minutes'] for s in summaries)
            print(f"\n  Total compute time: {total_time:.1f} minutes")
            phi_r_checks = [s['positive_phi_r_check'] for s in summaries]
            print(f"  Phi_R non-negative check passed: {sum(phi_r_checks)}/{len(phi_r_checks)}")

    if errors > 0:
        print("\nErrors:")
        for r in results:
            if r['status'] == 'error':
                print(f"  - {Path(r['input_path']).name}: {r['message']}")


def main():
    parser = argparse.ArgumentParser(
        description='Compute PhiID measures using Luppi et al. (2024) method',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process entire dataset (batch mode)
  python 04_compute_phiid_revised.py --dataset london_lsd

  # Process specific subject
  python 04_compute_phiid_revised.py --dataset berlin_ketamine --subject 14

  # Process specific task/session
  python 04_compute_phiid_revised.py --dataset london_lsd --task LSD
  python 04_compute_phiid_revised.py --dataset berlin_ketamine --session ketamine

  # Force reprocessing
  python 04_compute_phiid_revised.py --dataset london_lsd --overwrite

  # Single scan mode (explicit subject/task/run)
  python 04_compute_phiid_revised.py --dataset london_lsd --subject 01 --task LSD --run 3
        """
    )

    # Required
    parser.add_argument('--dataset', required=True,
                        help='Dataset name (e.g., london_lsd, berlin_ketamine)')

    # Filters (optional - if subject+task+run all given, runs single scan mode)
    parser.add_argument('--subject', help='Subject ID (without sub- prefix)')
    parser.add_argument('--session', help='Session ID (without ses- prefix)')
    parser.add_argument('--task', help='Task name')
    parser.add_argument('--run', type=int, help='Run number')

    # Output
    parser.add_argument('--output_dir',
                        help='Output directory (default: results/phiid_matrices)')

    # Processing options
    parser.add_argument('--tau', type=int, default=1,
                        help='Time delay in TRs (default: 1)')
    parser.add_argument('--method', default='gaussian',
                        help='PID method (default: gaussian)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output files')
    parser.add_argument('--project_dir',
                        help='Project directory (auto-detected if not specified)')

    args = parser.parse_args()

    print("=" * 60)
    print("PhiID Computation (Luppi et al. 2024)")
    print("Equations 4-7: Red, Syn, Phi, Phi_R")
    print("=" * 60)

    # Resolve project directory
    project_dir = Path(args.project_dir) if args.project_dir else SCRIPT_DIR.parent.parent
    paths_config = project_dir / 'configs' / 'paths.yaml'
    with open(paths_config, 'r') as f:
        paths = yaml.safe_load(f)

    # Directories
    timeseries_dir = Path(paths['output']['timeseries'])
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path(paths['output']['phiid_matrices'])
    output_dir.mkdir(parents=True, exist_ok=True)

    if not timeseries_dir.exists():
        print(f"ERROR: Timeseries directory not found: {timeseries_dir}")
        sys.exit(1)

    # Discover deconvolved files
    print(f"\nDiscovering deconvolved timeseries for dataset: {args.dataset}")
    deconv_files = discover_deconv_files(
        timeseries_dir,
        dataset=args.dataset,
        subject=args.subject,
        session=args.session,
        task=args.task
    )

    if not deconv_files:
        print("No deconvolved timeseries files found matching criteria.")
        print(f"  Dataset: {args.dataset}")
        print(f"  Subject filter: {args.subject}")
        print(f"  Session filter: {args.session}")
        print(f"  Task filter: {args.task}")
        print(f"  Directory: {timeseries_dir}")
        sys.exit(1)

    print(f"Found {len(deconv_files)} deconvolved timeseries files")
    for f in deconv_files:
        print(f"  - {f.name}")

    # Process each scan
    results = []
    for i, deconv_file in enumerate(deconv_files, 1):
        print(f"\n[{i}/{len(deconv_files)}] Processing {deconv_file.name}")
        result = process_single_scan(
            deconv_file=deconv_file,
            output_dir=output_dir,
            dataset=args.dataset,
            tau=args.tau,
            method=args.method,
            overwrite=args.overwrite,
            verbose=True
        )
        results.append(result)

    # Print summary
    print_summary(results)


if __name__ == '__main__':
    main()