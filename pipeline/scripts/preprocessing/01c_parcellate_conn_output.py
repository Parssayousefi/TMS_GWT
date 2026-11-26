#!/usr/bin/env python3
"""
Stage 1c of the CONN-based extraction pipeline.

Dataset-agnostic generalisation of scripts/test_against_luppi/parcellate_conn_output.py.
For each scan in the manifest, apply the 454-ROI atlas to CONN's cleaned 4D
BOLD with NO further processing (no detrending, no bandpass, no confound
regression, no standardisation) — CONN already did all the denoising in stage
1b, so this is pure spatial averaging per atlas label.

Outputs are written in production-identical format:
  results/timeseries/dataset-{name}_sub-{subject}[_ses-{session}]_task-{task}_run-{run}_timeseries.npy
  + matching _metadata.yaml sidecar

This makes the rest of the pipeline (step02 rsHRF, step04 PhiID, 05/06 group
analyses) work without modification.

After successful save, the intermediate cleaned 4D file is deleted to free
disk. Use --keep-cleaned to retain it.

Usage:
    python scripts/preprocessing/01c_parcellate_conn_output.py --dataset hcp
    python scripts/preprocessing/01c_parcellate_conn_output.py --dataset hcp --subjects 100307
    python scripts/preprocessing/01c_parcellate_conn_output.py --dataset basel_lam --keep-cleaned
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse
import json
import os
import sys

import numpy as np
import nibabel as nib
import yaml
from nilearn import maskers

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_DIR / "scripts" / "utils"))
sys.path.insert(0, str(PROJECT_DIR / "scripts" / "preprocessing"))

from dataset_handler import load_dataset_handler, list_active_datasets  # noqa: E402


CONN_INPUT_ROOT = PROJECT_DIR / "data" / "conn_input"
TIMESERIES_DIR = PROJECT_DIR / "results" / "timeseries"

# Canonical 454-ROI atlas built natively in MNI152NLin2009cAsym by
# scripts/archive/build_combined_atlas_nlin2009c.py — matches fMRIPrep's
# output space and therefore CONN's denoised BOLDs (also NLin2009cAsym).
# Loading directly avoids an importlib dance through the archived production
# extraction script and ensures we always parcellate against the right atlas
# for the CONN pipeline, not the legacy NLin6Asym variant.
ATLAS_PATH = PROJECT_DIR / "atlases" / "schaefer_400_tian_54_combined.nii.gz"
N_EXPECTED_ROIS = 454
EXPECTED_ATLAS_SHAPE = (97, 115, 97)  # NLin2009cAsym 2mm grid


def load_canonical_atlas() -> str:
    """Load and sanity-check the canonical NLin2009cAsym 454-ROI atlas.

    Verifies the atlas is on the expected NLin2009cAsym 2mm grid (shape
    97x115x97) — NOT the legacy NLin6Asym grid (shape 91x109x91). The CONN
    cleaned 4D BOLDs are in NLin2009cAsym (inherited from fMRIPrep), so using
    a NLin6Asym atlas would introduce a template mismatch and corrupt the
    parcellation.
    """
    if not ATLAS_PATH.exists():
        raise FileNotFoundError(
            f"Canonical atlas not found at {ATLAS_PATH}. "
            f"Run scripts/archive/build_combined_atlas_nlin2009c.py to build it."
        )
    img = nib.load(str(ATLAS_PATH))
    if img.shape != EXPECTED_ATLAS_SHAPE:
        raise RuntimeError(
            f"Atlas at {ATLAS_PATH} has shape {img.shape}, expected "
            f"{EXPECTED_ATLAS_SHAPE} (NLin2009cAsym 2mm). The legacy NLin6Asym "
            f"atlas (91x109x91) is incompatible with CONN's NLin2009cAsym BOLDs. "
            f"Rebuild via scripts/archive/build_combined_atlas_nlin2009c.py --force."
        )
    data = img.get_fdata()
    n_labels = len(np.unique(data[data > 0]))
    if n_labels != N_EXPECTED_ROIS:
        raise RuntimeError(
            f"Atlas at {ATLAS_PATH} has {n_labels} unique labels, "
            f"expected {N_EXPECTED_ROIS}."
        )
    return str(ATLAS_PATH)


def parcellate_one_scan(entry: dict, handler, atlas_path: str, intermediate_root: Path,
                         output_dir: Path, keep_cleaned: bool, keep_inputs: bool) -> dict:
    """Apply atlas to one CONN-cleaned 4D, save production-format timeseries."""
    scan_id = entry["scan_id"]
    subject = entry["subject"]
    session = entry["session"]
    task = entry["task"]
    run = int(entry["run"])

    print(f"\n{'=' * 60}")
    print(scan_id)
    print(f"{'=' * 60}")

    cleaned_path = intermediate_root / f"{scan_id}_conn_cleaned.nii"
    if not cleaned_path.exists():
        raise FileNotFoundError(
            f"Cleaned 4D BOLD not found at {cleaned_path}. "
            f"Has 01b_conn_denoise.m been run for this dataset?"
        )

    cleaned_img = nib.load(str(cleaned_path))
    print(f"  cleaned 4D: {cleaned_img.shape}  ({cleaned_path.stat().st_size / 1e9:.2f} GB)")

    # Pure spatial averaging — CONN has already done detrend + confound regression
    # + bandpass. Cleaning args off.
    masker = maskers.NiftiLabelsMasker(
        labels_img=atlas_path,
        standardize=False,
        detrend=False,
        low_pass=None,
        high_pass=None,
        t_r=None,
        smoothing_fwhm=None,
        verbose=1,
    )
    roi_ts = masker.fit_transform(cleaned_img)
    print(f"  ROI timeseries: {roi_ts.shape}")

    if roi_ts.shape[1] != 454:
        raise RuntimeError(
            f"{scan_id}: expected 454 ROIs, got {roi_ts.shape[1]}"
        )

    # ---- Production-format output filename via dataset_handler ----
    output_filename = handler.get_output_filename(
        subject=subject,
        task=task,
        run=run,
        session=session,
        suffix="timeseries",
        extension="npy",
    )
    output_path = output_dir / output_filename
    np.save(output_path, roi_ts)
    print(f"  saved: {output_path}")

    # ---- Metadata sidecar (matches production schema) ----
    art_scrubbing = bool(handler.preprocessing.get("art_scrubbing", False))
    denoising_confounds = "5 WM PCs + 5 CSF PCs + 6 motion + 6 first derivatives"
    alignment_reference = (
        "Luppi 2022 Nat Neurosci (paper 3) HCP methodology, applied uniformly"
    )
    if art_scrubbing:
        denoising_confounds += (
            " + per-outlier-volume ART scrubbing regressors (regressed, not censored)"
        )
        alignment_reference = (
            "Luppi 2022 Nat Neurosci (paper 3) HCP methodology, plus ART spike "
            "scrubbing regressors (deliberate stricter-than-Luppi deviation; "
            "see step01b_conn_denoise_cambridge_art.m)"
        )
    metadata = {
        "subject": subject,
        "session": session,
        "task": task,
        "run": run,
        "dataset": handler.dataset_name,
        "tr": float(entry["tr"]),
        "n_timepoints": int(roi_ts.shape[0]),
        "n_rois": int(roi_ts.shape[1]),
        "atlas": "schaefer_400_tian_54_combined.nii.gz",
        "atlas_space": "MNI152NLin2009cAsym",
        "atlas_resolution_mm": 2,
        "atlas_build_script": "scripts/archive/build_combined_atlas_nlin2009c.py",
        "denoising_method": "CONN_V25.b_aCompCor",
        "denoising_confounds": denoising_confounds,
        "bandpass": [0.008, 0.09],
        "detrending": "linear",
        "art_spike_regressors": art_scrubbing,
        "global_signal_regression": False,
        "alignment_reference": alignment_reference,
    }
    metadata_filename = handler.get_output_filename(
        subject=subject,
        task=task,
        run=run,
        session=session,
        suffix="timeseries_metadata",
        extension="yaml",
    )
    metadata_path = output_dir / metadata_filename
    with open(metadata_path, "w") as f:
        yaml.safe_dump(metadata, f, default_flow_style=False, sort_keys=False)
    print(f"  metadata: {metadata_path.name}")

    # ---- Cleanup intermediate cleaned 4D (keep ~3KB .npy, drop ~4GB .nii) ----
    if not keep_cleaned:
        try:
            cleaned_path.unlink()
            print(f"  removed intermediate: {cleaned_path.name}")
        except OSError as exc:
            print(f"  WARN: could not remove {cleaned_path}: {exc}")

    # ---- Cleanup per-scan input dir (~5 GB scratch) ----
    if not keep_inputs:
        scan_input_dir = CONN_INPUT_ROOT / handler.dataset_name / scan_id
        if scan_input_dir.exists():
            import shutil
            try:
                shutil.rmtree(scan_input_dir)
                print(f"  removed input scratch: {scan_input_dir}")
            except OSError as exc:
                print(f"  WARN: could not remove {scan_input_dir}: {exc}")

    return {
        "scan_id": scan_id,
        "n_timepoints": int(roi_ts.shape[0]),
        "n_rois": int(roi_ts.shape[1]),
        "output_file": str(output_path),
    }


def _parcellate_scan_worker(
    entry: dict,
    dataset: str,
    atlas_path: str,
    intermediate_root_str: str,
    output_dir_str: str,
    keep_cleaned: bool,
    keep_inputs: bool,
) -> dict:
    """ProcessPoolExecutor-safe wrapper. Reloads handler inside the worker
    (cheap) rather than pickling it across process boundaries."""
    handler = load_dataset_handler(dataset, project_dir=str(PROJECT_DIR))
    return parcellate_one_scan(
        entry, handler, atlas_path, Path(intermediate_root_str),
        Path(output_dir_str), keep_cleaned, keep_inputs,
    )


def parcellate_dataset(
    dataset: str,
    atlas_path: str,
    subject_filter: list = None,
    session_filter: list = None,
    task_filter: list = None,
    keep_cleaned: bool = False,
    keep_inputs: bool = False,
    n_workers: int = 1,
) -> tuple:
    """Parcellate one dataset end-to-end. Returns (n_done, n_failed)."""
    manifest_path = CONN_INPUT_ROOT / dataset / "manifest.json"
    if not manifest_path.exists():
        print(f"  SKIP dataset {dataset}: manifest not found at {manifest_path}")
        print(f"  (run scripts/preprocessing/01a_prepare_conn_inputs.py --dataset {dataset} first)")
        return 0, 0
    with open(manifest_path) as f:
        manifest = json.load(f)

    intermediate_root = PROJECT_DIR / "results" / "timeseries_intermediate" / dataset
    if not intermediate_root.exists():
        print(f"  SKIP dataset {dataset}: no intermediate dir at {intermediate_root}")
        print("  (run scripts/preprocessing/01b_conn_denoise.m first)")
        return 0, 0

    handler = load_dataset_handler(dataset, project_dir=str(PROJECT_DIR))

    # The ART 1b variant leaves art_outlier_summary.csv in the intermediate dir.
    # If that evidence is present but the config does not declare art_scrubbing,
    # every sidecar written below would falsely claim no ART — refuse to write
    # wrong provenance rather than warn and scroll past.
    art_marker = intermediate_root / "art_outlier_summary.csv"
    if art_marker.exists() and not handler.preprocessing.get("art_scrubbing", False):
        raise RuntimeError(
            f"{art_marker} exists (stage 1b ART variant was run for {dataset}) but "
            f"'preprocessing.art_scrubbing' is not set in the dataset config. Add "
            f"'art_scrubbing: true' to configs/datasets/{dataset}.yaml (or remove a "
            f"stale ART summary) so metadata sidecars record the true denoising."
        )

    scans_to_process = manifest["scans"]
    if subject_filter:
        scans_to_process = [s for s in scans_to_process if s["subject"] in subject_filter]
    if session_filter:
        scans_to_process = [s for s in scans_to_process if s["session"] in session_filter]
    if task_filter:
        scans_to_process = [s for s in scans_to_process if s["task"] in task_filter]

    print("=" * 60)
    print(f"Parcellating dataset: {dataset}")
    print(f"Scans in manifest: {len(manifest['scans'])}; after filters: {len(scans_to_process)}")
    print("=" * 60)

    if not scans_to_process:
        print(f"  No scans selected for dataset {dataset}.")
        return 0, 0

    results = []
    failures = []
    effective_workers = max(1, min(n_workers, len(scans_to_process)))
    if effective_workers > 1:
        print(f"  Parallel mode: {effective_workers} workers (ProcessPoolExecutor)")
        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            future_to_scan = {
                executor.submit(
                    _parcellate_scan_worker,
                    entry, dataset, atlas_path,
                    str(intermediate_root), str(TIMESERIES_DIR),
                    keep_cleaned, keep_inputs,
                ): entry["scan_id"]
                for entry in scans_to_process
            }
            for future in as_completed(future_to_scan):
                scan_id = future_to_scan[future]
                try:
                    results.append(future.result())
                except Exception as exc:  # noqa: BLE001
                    print(f"  ERROR {scan_id}: {type(exc).__name__}: {exc}")
                    failures.append((scan_id, str(exc)))
    else:
        for entry in scans_to_process:
            try:
                result = parcellate_one_scan(
                    entry, handler, atlas_path, intermediate_root,
                    TIMESERIES_DIR, keep_cleaned, keep_inputs,
                )
                results.append(result)
            except Exception as exc:  # noqa: BLE001
                print(f"  ERROR {entry['scan_id']}: {type(exc).__name__}: {exc}")
                failures.append((entry["scan_id"], str(exc)))

    print()
    print("-" * 60)
    print(f"Dataset {dataset}: {len(results)}/{len(scans_to_process)} parcellated, {len(failures)} failed")
    print("-" * 60)
    if failures:
        for scan_id, reason in failures:
            print(f"  [failed] {scan_id}: {reason}")
    return len(results), len(failures)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parcellate CONN-denoised 4D BOLD into production-format 454-ROI timeseries"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dataset", help="Single dataset name")
    group.add_argument("--all", action="store_true", help="Parcellate all active datasets in turn")
    parser.add_argument("--subjects", nargs="+", default=None, help="Optional subject ID filter")
    parser.add_argument("--sessions", nargs="+", default=None, help="Optional session ID filter")
    parser.add_argument("--tasks", nargs="+", default=None, help="Optional task name filter")
    parser.add_argument(
        "--keep-cleaned",
        action="store_true",
        help="Keep intermediate CONN-cleaned 4D files (default: delete after parcellation)",
    )
    parser.add_argument(
        "--keep-inputs",
        action="store_true",
        help="Keep per-scan input scratch dirs (default: delete after parcellation)",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=12,
        help=(
            "Number of parallel ProcessPoolExecutor workers (default: 12). "
            "Each worker peaks ~5-8 GB RAM (4.3 GB 4D nifti + masker working set). "
            "Use 1 for serial."
        ),
    )
    args = parser.parse_args()

    if args.n_workers > 1:
        # Prevent BLAS thread oversubscription — each worker process will inherit
        # these env vars at spawn time (Windows) and limit internal threading to 1.
        # Must be set BEFORE ProcessPoolExecutor spawns children.
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
        os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    if args.all:
        datasets = list_active_datasets(project_dir=str(PROJECT_DIR))
        print(f"--all mode: will parcellate {len(datasets)} datasets in sequence: {datasets}")
    else:
        datasets = [args.dataset]

    TIMESERIES_DIR.mkdir(parents=True, exist_ok=True)
    atlas_path = load_canonical_atlas()
    print(f"Atlas: {atlas_path} (NLin2009cAsym, 454 ROIs)")
    print(f"Output dir: {TIMESERIES_DIR}")
    print(f"Keep cleaned 4D: {args.keep_cleaned}")
    print(f"Keep input scratch: {args.keep_inputs}")
    print(f"Parallel workers: {args.n_workers}")

    totals = {"done": 0, "failed": 0}
    per_dataset = {}
    for ds in datasets:
        print()
        print("#" * 60)
        print(f"# DATASET {ds}")
        print("#" * 60)
        try:
            done, failed = parcellate_dataset(
                ds,
                atlas_path=atlas_path,
                subject_filter=args.subjects,
                session_filter=args.sessions,
                task_filter=args.tasks,
                keep_cleaned=args.keep_cleaned,
                keep_inputs=args.keep_inputs,
                n_workers=args.n_workers,
            )
        except Exception as exc:  # noqa: BLE001
            print(f"  FATAL: dataset {ds} failed: {type(exc).__name__}: {exc}")
            done, failed = 0, 0
        per_dataset[ds] = (done, failed)
        totals["done"] += done
        totals["failed"] += failed

    print()
    print("=" * 60)
    print(f"Stage 1c complete (datasets: {len(datasets)})")
    for ds, (done, failed) in per_dataset.items():
        print(f"  {ds:<20} done={done}  failed={failed}")
    print(f"  {'TOTAL':<20} done={totals['done']}  failed={totals['failed']}")
    print("=" * 60)
    print()
    print("Next step:")
    print("  MATLAB: run scripts/preprocessing/step02_hrf_deconvolution_matlab.m")


if __name__ == "__main__":
    main()
