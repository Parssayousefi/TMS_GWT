#!/usr/bin/env python3
"""
Parallel PhiID computation wrapper.

Launches independent Python subprocesses for each scan, avoiding Windows
ProcessPoolExecutor issues with native library conflicts. Each subprocess
runs 04_compute_phiid_revised.py with single-subject filters.

Usage:
  # Run HCP first (defines workspace for downstream analyses)
  python scripts/analysis/04b_compute_phiid_parallel.py --dataset hcp --max_workers 16

  # Run all datasets (skip-existing will skip already-computed scans)
  python scripts/analysis/04b_compute_phiid_parallel.py --max_workers 16

  # Run specific dataset with overwrite
  python scripts/analysis/04b_compute_phiid_parallel.py --dataset london_lsd --max_workers 16 --overwrite
"""

import os
import sys
import time
import argparse
import subprocess
import importlib.util
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Resolve paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
PHIID_SCRIPT = SCRIPT_DIR / "04_compute_phiid_revised.py"


def _load_phiid_module():
    """Load the PhiID module for discovery functions only (main process)."""
    spec = importlib.util.spec_from_file_location(
        "compute_phiid",
        PHIID_SCRIPT
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_single_scan(args):
    """
    Run PhiID for a single scan as an independent subprocess.

    Uses ThreadPoolExecutor to manage subprocess concurrency - each thread
    launches an independent Python process via subprocess.run().
    """
    deconv_file, output_dir, dataset, tau, method, overwrite = args
    filename = Path(deconv_file).name

    # Build command to call original script with filters matching this specific file
    cmd = [
        sys.executable,
        str(PHIID_SCRIPT),
        "--dataset", dataset,
        "--tau", str(tau),
        "--method", method,
        "--output_dir", str(output_dir),
        "--project_dir", str(PROJECT_DIR),
    ]
    if overwrite:
        cmd.append("--overwrite")

    # Extract subject/session/task/run from filename to filter to this exact scan
    import re
    name = Path(deconv_file).stem

    m = re.search(r'sub-([^_]+)', name)
    if m:
        cmd.extend(["--subject", m.group(1)])

    m = re.search(r'ses-([^_]+)', name)
    if m:
        cmd.extend(["--session", m.group(1)])

    m = re.search(r'task-([^_]+)', name)
    if m:
        cmd.extend(["--task", m.group(1)])

    m = re.search(r'run-(\d+)', name)
    if m:
        cmd.extend(["--run", m.group(1)])

    # Set thread limits for subprocess
    env = os.environ.copy()
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['OMP_NUM_THREADS'] = '1'
    env['MKL_NUM_THREADS'] = '1'
    env['NUMEXPR_NUM_THREADS'] = '1'

    start = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            timeout=1800  # 30 min timeout per scan
        )
        elapsed = time.time() - start

        if result.returncode == 0:
            return {
                'status': 'success',
                'filename': filename,
                'elapsed_minutes': round(elapsed / 60, 1),
                'message': ''
            }
        else:
            # Check if it was skipped (already computed)
            if 'already computed' in result.stdout.lower() or 'Skipping' in result.stdout:
                return {
                    'status': 'skipped',
                    'filename': filename,
                    'elapsed_minutes': 0,
                    'message': 'Already computed'
                }
            return {
                'status': 'error',
                'filename': filename,
                'elapsed_minutes': round(elapsed / 60, 1),
                'message': result.stderr[-500:] if result.stderr else result.stdout[-500:]
            }
    except subprocess.TimeoutExpired:
        return {
            'status': 'error',
            'filename': filename,
            'elapsed_minutes': 30,
            'message': 'Timeout (>30 min)'
        }
    except Exception as e:
        return {
            'status': 'error',
            'filename': filename,
            'elapsed_minutes': round((time.time() - start) / 60, 1),
            'message': str(e)
        }


def main():
    parser = argparse.ArgumentParser(
        description='Parallel PhiID computation via independent subprocesses',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run HCP first (workspace baseline)
  python 04b_compute_phiid_parallel.py --dataset hcp --max_workers 16

  # Run all datasets
  python 04b_compute_phiid_parallel.py --max_workers 16

  # Filter by dataset and condition
  python 04b_compute_phiid_parallel.py --dataset london_lsd --session LSD --max_workers 16
        """
    )

    parser.add_argument('--dataset', help='Dataset name filter (omit for all datasets)')
    parser.add_argument('--subject', help='Subject ID filter')
    parser.add_argument('--session', help='Session filter')
    parser.add_argument('--task', help='Task filter')
    parser.add_argument('--max_workers', type=int, default=16,
                        help='Number of parallel subprocesses (default: 16)')
    parser.add_argument('--tau', type=int, default=1,
                        help='Time delay in TRs (default: 1)')
    parser.add_argument('--method', default='gaussian',
                        help='PID method (default: gaussian)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output files')
    parser.add_argument('--output_dir', help='Output directory (default: from paths.yaml)')

    args = parser.parse_args()

    print("=" * 60)
    print("Parallel PhiID Computation (subprocess mode)")
    print(f"Workers: {args.max_workers}")
    print("=" * 60)

    # Load module in main process for discovery only
    mod = _load_phiid_module()

    import yaml
    paths_config = PROJECT_DIR / 'configs' / 'paths.yaml'
    with open(paths_config, 'r') as f:
        paths = yaml.safe_load(f)

    timeseries_dir = Path(paths['output']['timeseries'])
    output_dir = Path(args.output_dir) if args.output_dir else Path(paths['output']['phiid_matrices'])
    output_dir.mkdir(parents=True, exist_ok=True)

    if not timeseries_dir.exists():
        print(f"ERROR: Timeseries directory not found: {timeseries_dir}")
        sys.exit(1)

    # Discover deconvolved files
    print(f"\nDiscovering deconvolved timeseries...")
    if args.dataset:
        print(f"  Dataset filter: {args.dataset}")

    deconv_files = mod.discover_deconv_files(
        timeseries_dir,
        dataset=args.dataset,
        subject=args.subject,
        session=args.session,
        task=args.task
    )

    if not deconv_files:
        print("No deconvolved timeseries files found matching criteria.")
        sys.exit(1)

    # Pre-filter already-computed scans
    to_process = []
    skipped = 0
    for f in deconv_files:
        scan_info = mod.parse_scan_info(f)
        ds = scan_info.get('dataset', args.dataset or 'unknown')
        if not args.overwrite and mod.check_phiid_exists(output_dir, scan_info, ds):
            skipped += 1
        else:
            to_process.append((f, ds))

    print(f"\nFound {len(deconv_files)} deconvolved files")
    print(f"  Already computed (skipping): {skipped}")
    print(f"  To process: {len(to_process)}")

    if not to_process:
        print("\nAll scans already computed. Use --overwrite to reprocess.")
        return

    # Build task arguments
    task_args = [
        (str(f), str(output_dir), ds, args.tau, args.method, args.overwrite)
        for f, ds in to_process
    ]

    # Process using ThreadPoolExecutor managing subprocesses
    # (threads manage subprocess.run calls; actual computation is in independent processes)
    print(f"\nStarting {len(to_process)} scans with {args.max_workers} parallel subprocesses...")
    start_time = time.time()
    results = []
    completed = 0

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        future_to_idx = {
            executor.submit(run_single_scan, ta): i for i, ta in enumerate(task_args)
        }

        for future in as_completed(future_to_idx):
            completed += 1
            try:
                result = future.result()
                results.append(result)

                status = result['status']
                fname = result['filename']
                elapsed = result.get('elapsed_minutes', 0)

                if status == 'success':
                    print(f"  [{completed}/{len(to_process)}] OK: {fname} ({elapsed} min)")
                elif status == 'skipped':
                    print(f"  [{completed}/{len(to_process)}] SKIP: {fname}")
                else:
                    print(f"  [{completed}/{len(to_process)}] ERROR: {fname}: {result['message'][:100]}")

            except Exception as e:
                completed_result = {'status': 'error', 'filename': '?', 'message': str(e)}
                results.append(completed_result)
                print(f"  [{completed}/{len(to_process)}] EXCEPTION: {e}")

    total_elapsed = time.time() - start_time

    # Summary
    success = sum(1 for r in results if r['status'] == 'success')
    error_count = sum(1 for r in results if r['status'] == 'error')
    skip_count = skipped + sum(1 for r in results if r['status'] == 'skipped')

    print("\n" + "=" * 60)
    print("PARALLEL PHIID COMPUTATION SUMMARY")
    print("=" * 60)
    print(f"  Successful:  {success}")
    print(f"  Skipped:     {skip_count}")
    print(f"  Errors:      {error_count}")
    print(f"  Total scans: {len(deconv_files)}")
    print(f"  Wall clock:  {total_elapsed/60:.1f} minutes")
    print(f"  Workers:     {args.max_workers}")

    if success > 0:
        scan_times = [r['elapsed_minutes'] for r in results if r['status'] == 'success']
        total_cpu = sum(scan_times)
        print(f"\n  Total CPU time:  {total_cpu:.1f} minutes")
        print(f"  Avg per scan:    {total_cpu/len(scan_times):.1f} minutes")
        print(f"  Speedup:         {total_cpu/max(total_elapsed/60, 0.01):.1f}x")

    if error_count > 0:
        print("\nErrors:")
        for r in results:
            if r['status'] == 'error':
                print(f"  - {r['filename']}: {r['message'][:200]}")

    # Verify output
    import glob
    phi_r_count = len(glob.glob(str(output_dir / "*_phi_r.npy")))
    print(f"\nTotal phi_r.npy files in output: {phi_r_count}")


if __name__ == '__main__':
    main()
