#!/usr/bin/env python3
"""
Stage 1a of the CONN-based extraction pipeline.

Dataset-agnostic generalisation of scripts/test_against_luppi/prepare_conn_inputs.py.
For each scan in a dataset (one BOLD = one scan, identified by subject + optional
session + task + run), prepare the minimum inputs CONN V25.b needs:
  - Trim BOLD per dataset config + QC truncation rules.
  - Decompress trimmed BOLD and T1w to plain .nii in a per-scan scratch dir.
  - Extract a 6-column motion text file from fMRIPrep confounds, sliced to match
    the trimmed BOLD.
  - Write a single manifest.json covering all valid scans for the dataset.

The MATLAB stage (01b_conn_denoise.m) consumes the manifest. The parcellate
stage (01c_parcellate_conn_output.py) writes drop-in production timeseries
files matching the existing format from 01_extract_timeseries_luppi_method.py.

Usage:
    python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset hcp
    python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset hcp --subjects 100307
    python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset basel_lam --subjects 001
    python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset london_lsd --tasks LSD
"""

from pathlib import Path
import argparse
import gzip
import importlib.util
import json
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import nibabel as nib
import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_DIR / "scripts" / "utils"))
sys.path.insert(0, str(PROJECT_DIR / "scripts" / "preprocessing"))

from dataset_handler import load_dataset_handler, list_active_datasets  # noqa: E402
from qc_loader import QCExclusions  # noqa: E402

# Reuse production helpers via importlib so trim semantics stay identical to
# 01_extract_timeseries_luppi_method.py.
_prod_path = PROJECT_DIR / "scripts" / "archive" / "01_extract_timeseries_luppi_method.py"
_spec = importlib.util.spec_from_file_location("_prod_extract", str(_prod_path))
_prod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_prod)
remove_boundary_volumes = _prod.remove_boundary_volumes
select_volume_range = _prod.select_volume_range


CONN_INPUT_ROOT = PROJECT_DIR / "data" / "conn_input"
MOTION_COLS = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]


def decompress_to_nii(src_path: Path, dst_path: Path) -> None:
    """Decompress a .nii.gz to a plain .nii in dst_path, overwriting if present.

    Always overwrites so re-runs with different trim settings produce the
    correct output rather than a stale file from a prior invocation.
    """
    with gzip.open(str(src_path), "rb") as f_in:
        with open(str(dst_path), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out, length=64 * 1024 * 1024)


def find_t1w_file(handler, subject: str, session: str = None) -> Path:
    """Locate fMRIPrep T1w preproc output in MNI152NLin2009cAsym res-2 space.

    Tries subject-level anat first (most common), then session-level anat as
    fallback for datasets where anat was acquired per session. CONN runs its own
    segmentation on this T1 to derive WM/CSF masks for aCompCor.
    """
    fmriprep_dir = Path(handler.config["paths"]["local"]["fmriprep_dir"])

    candidates = [
        fmriprep_dir / f"sub-{subject}" / "anat" /
        f"sub-{subject}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz",
    ]
    if session is not None:
        candidates.append(
            fmriprep_dir / f"sub-{subject}" / f"ses-{session}" / "anat" /
            f"sub-{subject}_ses-{session}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz"
        )

    for cand in candidates:
        if cand.exists():
            return cand

    # Fallback: some datasets preserve extra BIDS entities from the source T1
    # (e.g. cambridge_propofol T1 source was sub-XX_acq-MPRAGE_run-01_T1w.nii.gz,
    # which fMRIPrep propagates into the MNI-space output filename). Glob for
    # the canonical suffix, scoped to the subject's anat dirs.
    glob_suffix = "*_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz"
    glob_dirs = [fmriprep_dir / f"sub-{subject}" / "anat"]
    if session is not None:
        glob_dirs.append(fmriprep_dir / f"sub-{subject}" / f"ses-{session}" / "anat")
    # Last resort: any other session's anat dir. Needed when anat exists in only
    # one session and fMRIPrep therefore emits session-level anat outputs (e.g.
    # monash_psilo: bids_filter.json restricts anat to ses-baseline, so
    # ses-psilocybin scans must borrow the ses-baseline T1). Sorted for
    # determinism; earlier candidates always win.
    glob_dirs.extend(sorted(
        (fmriprep_dir / f"sub-{subject}").glob("ses-*/anat")
    ))
    for d in glob_dirs:
        if d.is_dir():
            matches = sorted(d.glob(glob_suffix))
            if matches:
                return matches[0]

    raise FileNotFoundError(
        f"sub-{subject}" + (f" ses-{session}" if session else "") +
        f": T1w not found. Tried literal: {[str(c) for c in candidates]}; "
        f"globbed: {[str(d / glob_suffix) for d in glob_dirs]}"
    )


def make_scan_id(subject: str, session, task: str, run) -> str:
    """Build a unique scan identifier used as the per-scan input directory name.

    For datasets without sessions and a single (task, run) combo per subject
    (e.g. HCP), this collapses to ``sub-{subject}``. For multi-scan-per-subject
    datasets it includes session/task/run to avoid file collisions.
    """
    parts = [f"sub-{subject}"]
    if session is not None:
        parts.append(f"ses-{session}")
    parts.append(f"task-{task}")
    parts.append(f"run-{run}")
    return "_".join(parts)


def prepare_scan(
    scan: dict,
    handler,
    qc: QCExclusions,
    auto_truncate: bool,
    dataset_input_root: Path,
):
    """Process one scan end-to-end. Returns the manifest entry for stage 1b,
    or None if the scan is QC-skipped (EXCLUDE, or TRUNCATE with auto_truncate=False)."""
    subject = scan["subject"]
    session = scan["session"]  # may be None for HCP-style datasets
    task = scan["task"]
    run = int(scan["run"])
    func_file = Path(scan["func_file"])
    confounds_file = Path(scan["confounds_file"])

    scan_id = make_scan_id(subject, session, task, run)

    print(f"\n{'=' * 60}")
    print(scan_id)
    print(f"{'=' * 60}")

    # ---- QC gate (mirrors production main loop) ----
    # Pass task explicitly so datasets that encode conditions in the task name
    # (has_sessions=false, e.g. cambridge_propofol) get per-scan QC decisions
    # rather than a worst-case collapse over all tasks for the subject.
    qc_status = qc.get_status(subject, session, task)
    truncation = None
    if qc_status == "EXCLUDE":
        print(f"[QC SKIP] {scan_id}: EXCLUDED — {qc.get_exclusion_reason(subject, session, task)}")
        return None
    elif qc_status == "TRUNCATE":
        if auto_truncate:
            truncation = qc.get_truncation(subject, session, task)
            if truncation is not None:
                print(
                    f"[QC TRUNCATE] {scan_id}: volumes "
                    f"{truncation.start_volume}–{truncation.end_volume} "
                    f"({truncation.n_volumes} vol) — {truncation.reason}"
                )
        else:
            print(f"[QC SKIP] {scan_id}: TRUNCATE but auto_truncate=False — {qc.get_exclusion_reason(subject, session, task)}")
            return None

    # ---- Resolve T1 ----
    t1_src = find_t1w_file(handler, subject, session)

    preprocessing = handler.preprocessing
    volumes_to_remove = preprocessing.get("volumes_to_remove", 5)
    volumes_to_remove_last = preprocessing.get("volumes_to_remove_last", 0)
    tr = float(preprocessing.get("tr"))

    # ---- Per-scan input directory ----
    sub_dir = dataset_input_root / scan_id
    sub_dir.mkdir(parents=True, exist_ok=True)
    print(f"  output dir:  {sub_dir}")
    print(f"  func source: {func_file}")
    if truncation is not None:
        print(f"  trim:        QC truncation to volumes {truncation.start_volume}–{truncation.end_volume}")
    else:
        print(f"  trim:        first {volumes_to_remove} + last {volumes_to_remove_last} volumes")
    print(f"  TR:          {tr}s")

    # ---- 1. Trim BOLD ----
    print("Step 1/4: trimming BOLD...")
    if truncation is not None:
        trimmed_temp_path, n_volumes, actual_start, actual_end = select_volume_range(
            str(func_file),
            start_volume=truncation.start_volume,
            end_volume=truncation.end_volume,
            volumes_to_remove_first=volumes_to_remove,
        )
        confounds_start = actual_start
        confounds_end = actual_end
    else:
        trimmed_temp_path, n_volumes = remove_boundary_volumes(
            str(func_file),
            n_first=volumes_to_remove,
            n_last=volumes_to_remove_last,
        )
        confounds_start = volumes_to_remove
        confounds_end = None
    trimmed_temp_path = Path(trimmed_temp_path)
    print(f"  trimmed BOLD has {n_volumes} timepoints")

    bold_dst = sub_dir / "bold_trimmed.nii"
    if trimmed_temp_path.suffix == ".gz" or str(trimmed_temp_path).endswith(".nii.gz"):
        print(f"  decompressing trimmed BOLD -> {bold_dst}")
        decompress_to_nii(trimmed_temp_path, bold_dst)
    else:
        print(f"  copying trimmed BOLD -> {bold_dst}")
        shutil.copy2(str(trimmed_temp_path), str(bold_dst))

    bold_img = nib.load(str(bold_dst))
    print(f"  bold_trimmed.nii shape: {bold_img.shape}")

    try:
        trimmed_temp_path.unlink()
    except (PermissionError, FileNotFoundError):
        pass

    # ---- 2. Decompress structural T1 ----
    print("Step 2/4: decompressing structural T1...")
    t1_dst = sub_dir / "T1w.nii"
    decompress_to_nii(t1_src, t1_dst)
    print(f"  T1w: {t1_dst.name} ({t1_dst.stat().st_size / 1e6:.1f} MB)")

    # ---- 3. Extract & trim motion parameters ----
    # Only 6 motion params are bridged here; 5 WM + 5 CSF aCompCor are re-computed
    # by CONN in 01b from the T1 structural_segment (matches Luppi 2022/2024).
    print("Step 3/4: extracting motion parameters...")
    confounds_df = pd.read_csv(confounds_file, sep="\t")
    missing_cols = [c for c in MOTION_COLS if c not in confounds_df.columns]
    if missing_cols:
        raise RuntimeError(f"{scan_id}: missing motion columns {missing_cols}")

    motion_df = confounds_df[MOTION_COLS]
    if confounds_end is not None:
        motion_df = motion_df.iloc[confounds_start:confounds_end].fillna(0)
    else:
        end_idx = (
            len(motion_df) - volumes_to_remove_last
            if volumes_to_remove_last > 0
            else len(motion_df)
        )
        motion_df = motion_df.iloc[confounds_start:end_idx].fillna(0)

    if len(motion_df) != n_volumes:
        raise RuntimeError(
            f"{scan_id}: motion length {len(motion_df)} != trimmed BOLD length {n_volumes}"
        )

    motion_dst = sub_dir / "motion.txt"
    np.savetxt(str(motion_dst), motion_df.values, fmt="%.8f", delimiter=" ")
    print(f"  motion.txt: {len(motion_df)} rows x 6 columns")

    # ---- 4. Build manifest entry ----
    print("Step 4/4: building manifest entry")
    entry = {
        "dataset":      handler.dataset_name,
        "scan_id":      scan_id,
        "subject":      subject,
        "session":      session,
        "task":         task,
        "run":          run,
        "tr":           tr,
        "n_volumes":    int(n_volumes),
        "bold_trimmed": str(bold_dst).replace("\\", "/"),
        "T1w":          str(t1_dst).replace("\\", "/"),
        "motion_txt":   str(motion_dst).replace("\\", "/"),
    }
    print(f"  done: {scan_id}")
    return entry


def _manifest_entry_from_existing(scan_id: str, handler, dataset_input_root: Path) -> dict:
    """Reconstruct a manifest entry for a scan whose inputs already exist on disk."""
    sub_dir = dataset_input_root / scan_id
    bold_path = sub_dir / "bold_trimmed.nii"
    t1_path = sub_dir / "T1w.nii"
    motion_path = sub_dir / "motion.txt"

    parts = scan_id.split("_")
    subject = parts[0].replace("sub-", "")
    session = None
    task = None
    run = 1
    for p in parts[1:]:
        if p.startswith("ses-"):
            session = p.replace("ses-", "")
        elif p.startswith("task-"):
            task = p.replace("task-", "")
        elif p.startswith("run-"):
            run = int(p.replace("run-", ""))

    tr = float(handler.preprocessing.get("tr"))
    img = nib.load(str(bold_path))
    n_volumes = int(img.shape[3]) if len(img.shape) == 4 else 0

    return {
        "dataset": handler.dataset_name,
        "scan_id": scan_id,
        "subject": subject,
        "session": session,
        "task": task,
        "run": run,
        "tr": tr,
        "n_volumes": n_volumes,
        "bold_trimmed": str(bold_path).replace("\\", "/"),
        "T1w": str(t1_path).replace("\\", "/"),
        "motion_txt": str(motion_path).replace("\\", "/"),
    }


def _rebuild_manifest_from_dirs(dataset, handler, scan_ids, dataset_input_root):
    """Build a manifest dict from already-prepared scan directories."""
    manifest = {"dataset": dataset, "scans": []}
    for scan_id in sorted(scan_ids):
        try:
            entry = _manifest_entry_from_existing(scan_id, handler, dataset_input_root)
            manifest["scans"].append(entry)
        except Exception as exc:  # noqa: BLE001
            print(f"  WARN: could not rebuild manifest entry for {scan_id}: {exc}")
    return manifest


def _prepare_scan_worker(scan, dataset, auto_truncate, dataset_input_root_str, project_dir_str):
    """ProcessPoolExecutor-safe wrapper. Reloads handler/qc inside the worker
    (cheap) rather than pickling complex objects across process boundaries."""
    handler = load_dataset_handler(dataset, project_dir=project_dir_str)
    qc = QCExclusions(dataset, project_dir=Path(project_dir_str))
    return prepare_scan(scan, handler, qc, auto_truncate, Path(dataset_input_root_str))


def prepare_dataset(
    dataset: str,
    subject_filter: list = None,
    session_filter: list = None,
    task_filter: list = None,
    manifest_path: Path = None,
    n_workers: int = 1,
    skip_done: bool = False,
) -> tuple:
    """Process one dataset end-to-end. Returns (n_scans_ready, n_skipped, n_failed)."""
    handler = load_dataset_handler(dataset, project_dir=str(PROJECT_DIR))
    qc = QCExclusions(dataset, project_dir=PROJECT_DIR)
    auto_truncate = handler.config.get("qc", {}).get("auto_truncate", True)

    dataset_input_root = CONN_INPUT_ROOT / dataset
    dataset_input_root.mkdir(parents=True, exist_ok=True)

    if manifest_path is None:
        manifest_path = dataset_input_root / "manifest.json"

    all_scans = handler.discover_all_scans()
    scans = all_scans
    if subject_filter:
        scans = [s for s in scans if s["subject"] in subject_filter]
    if session_filter:
        scans = [s for s in scans if s["session"] in session_filter]
    if task_filter:
        scans = [s for s in scans if s["task"] in task_filter]

    print("=" * 60)
    print(f"CONN input preparation — dataset: {dataset}")
    print(f"Discovered scans: {len(all_scans)}; after filters: {len(scans)}")
    print(f"Output root: {dataset_input_root}")
    print(f"QC auto_truncate: {auto_truncate}")
    print("=" * 60)

    if not scans:
        print(f"No scans selected for dataset {dataset}. Check filters.")
        return 0, 0, 0

    # ---- Skip already-prepared scans ----
    skipped_done = []
    if skip_done:
        remaining = []
        for scan in scans:
            scan_id = make_scan_id(scan["subject"], scan["session"], scan["task"], int(scan["run"]))
            done_marker = dataset_input_root / scan_id / "bold_trimmed.nii"
            if done_marker.exists():
                skipped_done.append(scan_id)
            else:
                remaining.append(scan)
        if skipped_done:
            print(f"  skip-done: {len(skipped_done)} scans already prepared, {len(remaining)} remaining")
        scans = remaining

    if not scans and skipped_done:
        print(f"All {len(skipped_done)} scans already prepared — nothing to do.")
        # Rebuild manifest from existing dirs so it stays complete
        manifest = _rebuild_manifest_from_dirs(dataset, handler, skipped_done, dataset_input_root)
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"Manifest refreshed: {manifest_path} ({len(manifest['scans'])} scans)")
        return len(manifest["scans"]), 0, 0

    manifest = {"dataset": dataset, "scans": []}
    failures = []
    skipped = []

    effective_workers = max(1, min(n_workers, len(scans)))
    if effective_workers > 1:
        print(f"  Parallel mode: {effective_workers} workers (ProcessPoolExecutor)")
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            future_to_scan = {
                executor.submit(
                    _prepare_scan_worker,
                    scan, dataset, auto_truncate,
                    str(dataset_input_root), str(PROJECT_DIR),
                ): make_scan_id(scan["subject"], scan["session"], scan["task"], int(scan["run"]))
                for scan in scans
            }
            for future in as_completed(future_to_scan):
                scan_label = future_to_scan[future]
                try:
                    entry = future.result()
                    if entry is None:
                        skipped.append(scan_label)
                    else:
                        manifest["scans"].append(entry)
                except Exception as exc:  # noqa: BLE001
                    print(f"  ERROR {scan_label}: {type(exc).__name__}: {exc}")
                    failures.append((scan_label, str(exc)))
    else:
        for scan in scans:
            scan_label = make_scan_id(scan["subject"], scan["session"], scan["task"], int(scan["run"]))
            try:
                entry = prepare_scan(scan, handler, qc, auto_truncate, dataset_input_root)
                if entry is None:
                    skipped.append(scan_label)
                else:
                    manifest["scans"].append(entry)
            except Exception as exc:  # noqa: BLE001
                print(f"  ERROR {scan_label}: {type(exc).__name__}: {exc}")
                failures.append((scan_label, str(exc)))

    # Merge previously-done scans into the manifest so it covers all valid scans
    if skipped_done:
        done_manifest = _rebuild_manifest_from_dirs(dataset, handler, skipped_done, dataset_input_root)
        manifest["scans"].extend(done_manifest["scans"])

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print()
    print("-" * 60)
    n_new = len(manifest["scans"]) - len(skipped_done)
    print(f"Dataset {dataset} done: {n_new} new, {len(skipped_done)} already done, "
          f"{len(skipped)} QC-skipped, {len(failures)} failed "
          f"(manifest: {len(manifest['scans'])} total)")
    print(f"Manifest: {manifest_path}")
    print("-" * 60)
    if skipped:
        for s in skipped:
            print(f"  [skipped] {s}")
    if failures:
        for lbl, reason in failures:
            print(f"  [failed]  {lbl}: {reason}")
    return len(manifest["scans"]), len(skipped), len(failures)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare fMRIPrep inputs for CONN denoising (production, dataset-agnostic)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dataset", help="Single dataset name (hcp, london_lsd, basel_lam, basel_lpm, berlin_ketamine)")
    group.add_argument("--all", action="store_true", help="Process all active datasets in turn (per dataset_handler.list_active_datasets)")
    parser.add_argument("--subjects", nargs="+", default=None, help="Optional subject ID filter (applied to every dataset when used with --all)")
    parser.add_argument("--sessions", nargs="+", default=None, help="Optional session ID filter")
    parser.add_argument("--tasks", nargs="+", default=None, help="Optional task name filter")
    parser.add_argument(
        "--manifest_path",
        default=None,
        help="Override manifest path. Only valid when --dataset is given (incompatible with --all).",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=12,
        help=(
            "Number of parallel ProcessPoolExecutor workers (default: 12). "
            "Each worker peaks ~1 GB RAM (BOLD decompression). Use 1 for serial."
        ),
    )
    parser.add_argument(
        "--skip-done",
        action="store_true",
        help="Skip scans whose bold_trimmed.nii already exists in the output directory.",
    )
    args = parser.parse_args()

    if args.all and args.manifest_path:
        parser.error("--manifest_path is incompatible with --all (each dataset needs its own manifest)")

    if args.n_workers > 1:
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
        os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    if args.all:
        datasets = list_active_datasets(project_dir=str(PROJECT_DIR))
        print(f"--all mode: will process {len(datasets)} datasets in sequence: {datasets}")
    else:
        datasets = [args.dataset]

    totals = {"ready": 0, "skipped": 0, "failed": 0}
    per_dataset = {}
    for ds in datasets:
        print()
        print("#" * 60)
        print(f"# DATASET {ds}")
        print("#" * 60)
        try:
            ready, skipped, failed = prepare_dataset(
                ds,
                subject_filter=args.subjects,
                session_filter=args.sessions,
                task_filter=args.tasks,
                manifest_path=Path(args.manifest_path) if args.manifest_path else None,
                n_workers=args.n_workers,
                skip_done=args.skip_done,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"  FATAL: dataset {ds} failed: {type(exc).__name__}: {exc}")
            ready, skipped, failed = 0, 0, 0
        per_dataset[ds] = (ready, skipped, failed)
        totals["ready"] += ready
        totals["skipped"] += skipped
        totals["failed"] += failed

    print()
    print("=" * 60)
    print(f"CONN input prep complete (datasets: {len(datasets)})")
    for ds, (ready, skipped, failed) in per_dataset.items():
        print(f"  {ds:<20} ready={ready}  skipped={skipped}  failed={failed}")
    print(f"  {'TOTAL':<20} ready={totals['ready']}  skipped={totals['skipped']}  failed={totals['failed']}")
    print("=" * 60)
    print()
    print("Next: open MATLAB and run scripts/preprocessing/01b_conn_denoise.m")
    if args.all:
        print('  (set DATASET = "" in the script to auto-loop over all prepared manifests)')


if __name__ == "__main__":
    main()
