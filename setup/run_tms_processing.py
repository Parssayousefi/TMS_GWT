"""Resume the frozen TMS processing chain in disk-bounded subject batches.

Calls production scripts without changing their analytical settings. MATLAB run
copies differ only in DATASET/N_WORKERS. No group contrasts or hypothesis tests
are run. Technical validation precedes deletion of regenerable CONN scratch.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import time
from datetime import datetime, timezone

import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "pipeline"
ALLOWED = ("aachen_itbs", "sudmex_tms")
MATLAB = Path(r"C:\Program Files\MATLAB\R2025a\bin\matlab.exe")
CONN_MASK = Path(r"C:\Users\YousefiP\software\conn25b\conn\utils\surf\mask.volume.brainmask.nii")
sys.path.insert(0, str(PIPE / "scripts" / "utils"))
from dataset_handler import load_dataset_handler
from qc_loader import QCExclusions


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    # Windows readers/antivirus can briefly hold the target without delete
    # sharing. Keep the old JSON intact and retry the atomic replacement.
    for attempt in range(20):
        try:
            temporary.replace(path)
            return
        except PermissionError:
            if attempt == 19:
                raise
            time.sleep(min(0.05 * 2**attempt, 0.5))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def below(path, parent):
    resolved, base = path.resolve(), parent.resolve()
    if resolved == base or base not in resolved.parents:
        raise RuntimeError(f"Path outside permitted scratch tree: {resolved}")
    return resolved


def scratch_target(path, parent):
    # Reject junctions/symlinks before resolving them, including parent links.
    lexical, base = Path(os.path.abspath(path)), Path(os.path.abspath(parent))
    lexical.relative_to(base)
    for component in (lexical, *lexical.parents):
        if component.exists() or component.is_symlink():
            attributes = getattr(component.lstat(), "st_file_attributes", 0)
            if component.is_symlink() or attributes & stat.FILE_ATTRIBUTE_REPARSE_POINT:
                raise RuntimeError(f"Refusing linked scratch path: {component}")
        if component == base:
            break
    return below(lexical, base)


def scan_id(scan):
    return (f"sub-{scan['subject']}_ses-{scan['session']}_"
            f"task-{scan['task']}_run-{int(scan['run'])}")


def series_path(handler, scan):
    return PIPE / "results" / "timeseries" / handler.get_output_filename(
        subject=scan["subject"], session=scan["session"], task=scan["task"],
        run=int(scan["run"]), suffix="timeseries", extension="npy")


def validate_series(path, n_volumes=None):
    values = np.load(path, mmap_mode="r")
    if values.ndim != 2 or values.shape[1] != 454:
        raise RuntimeError(f"Unexpected regional timeseries shape: {path}: {values.shape}")
    if n_volumes is not None and values.shape[0] != n_volumes:
        raise RuntimeError(f"Timepoint mismatch: {path}")
    if values.shape[0] < 240 or not np.isfinite(values).all():
        raise RuntimeError(f"Invalid or unexpectedly short regional timeseries: {path}")
    zero_count = int((np.std(values, axis=0) < 1e-10).sum())
    if zero_count == 454:
        raise RuntimeError(f"All regional signals are constant: {path}")
    if path.name.endswith("_timeseries.npy"):
        metadata = path.with_name(path.name.replace("_timeseries.npy", "_timeseries_metadata.yaml"))
        if not metadata.is_file():
            raise RuntimeError(f"Missing regional timeseries metadata: {metadata}")
    elif path.name.endswith("_timeseries_deconv.npy"):
        prefix = path.name.removesuffix("_timeseries_deconv.npy")
        for suffix in ("_hrf.mat", "_timeseries_deconv_metadata.yaml"):
            if not path.with_name(prefix + suffix).is_file():
                raise RuntimeError(f"Incomplete deconvolution sidecars: {path}")
    return {"file": str(path), "shape": list(values.shape),
            "zero_variance_rois": zero_count, "sha256": sha(path)}


def validate_cleaned(path, entry):
    # Production Setup.voxelresolution=1 resamples to CONN's mask grid;
    # it does not preserve the fMRIPrep voxel lattice or perform a new warp.
    cleaned, reference = nib.load(path), nib.load(CONN_MASK)
    if cleaned.shape != (*reference.shape[:3], entry["n_volumes"]):
        raise RuntimeError(f"Unexpected cleaned image dimensions: {path}: {cleaned.shape}")
    if not np.allclose(cleaned.affine, reference.affine):
        raise RuntimeError(f"Cleaned image differs from CONN analysis grid: {path}")
    return {"file": str(path), "shape": list(cleaned.shape),
            "affine": cleaned.affine.tolist(), "reference_mask": str(CONN_MASK)}


def run_command(command, log, env, state_path, stage):
    write_json(state_path, {"updated_utc": stamp(), "status": "running",
                           "stage": stage, "log": str(log), "command": command})
    print(f"[{stamp()}] {stage}; log={log}", flush=True)
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=str(PIPE), env=env,
                                stdout=stream, stderr=subprocess.STDOUT,
                                creationflags=subprocess.CREATE_NO_WINDOW)
    if result.returncode:
        raise RuntimeError(f"{stage} exited {result.returncode}; see {log}")


def matlab_copy(source, destination, workers, dataset=None):
    original = source.read_text(encoding="utf-8")
    updated, count = re.subn(r"(?m)^N_WORKERS = \d+;", f"N_WORKERS = {workers};", original)
    if count != 1:
        raise RuntimeError(f"Expected one N_WORKERS setting: {source}")
    changes = {"N_WORKERS": workers}
    if dataset is not None:
        if dataset not in ALLOWED:
            raise RuntimeError("Only TMS datasets are permitted")
        updated, count = re.subn(r"(?m)^DATASET\s*= '[^']*';", f"DATASET   = '{dataset}';", updated)
        if count != 1:
            raise RuntimeError(f"Expected one DATASET setting: {source}")
        changes["DATASET"] = dataset
    destination.write_text(updated, encoding="utf-8")
    dependencies = []
    if source.name == "step02_hrf_deconvolution_matlab.m":
        # The production script finds this sibling helper via its directory.
        # Preserve that dependency when the runtime copy lives in a log folder.
        helper = source.with_name("parsave_hrf.m")
        copied_helper = destination.with_name(helper.name)
        shutil.copy2(helper, copied_helper)
        dependencies.append({"source": str(helper), "copy": str(copied_helper),
                             "sha256": sha(copied_helper)})
    write_json(destination.with_suffix(".provenance.json"), {
        "source": str(source), "source_sha256": sha(source),
        "run_copy": str(destination), "run_sha256": sha(destination),
        "runtime_settings": changes, "analytical_settings_changed": False,
        "dependencies": dependencies})


def matlab_command(script):
    return [str(MATLAB), "-wait", "-singleCompThread", "-batch",
            "set(0,'DefaultFigureVisible','off'); run('" + script.as_posix() + "');"]


def clean_batch(dataset, entries, record):
    if dataset not in ALLOWED:
        raise RuntimeError("Only TMS scratch may be removed")
    input_root = PIPE / "data" / "conn_input" / dataset
    intermediate = PIPE / "results" / "timeseries_intermediate" / dataset
    targets = []
    for entry in entries:
        expected_dir = scratch_target(input_root / entry["scan_id"], input_root)
        if below(Path(entry["bold_trimmed"]), input_root).parent != expected_dir:
            raise RuntimeError("Manifest scratch path does not match its scan ID")
        targets.append(expected_dir)
        targets.append(scratch_target(intermediate / (entry["scan_id"] + "_conn_cleaned.nii"), intermediate))
    # These files are created by CONN and can exceed the original input size.
    targets.append(scratch_target(intermediate / "conn_project" / f"conn_{dataset}", intermediate))
    write_json(record, {"verified_utc": stamp(), "dataset": dataset,
                        "targets": [str(p) for p in targets], "status": "validated"})
    for target in targets:
        if target.is_symlink():
            raise RuntimeError(f"Refusing linked scratch target: {target}")
        if target.is_dir():
            shutil.rmtree(target)
        elif target.exists():
            target.unlink()
    write_json(record, {"completed_utc": stamp(), "dataset": dataset,
                        "targets": [str(p) for p in targets], "status": "removed"})


def process_dataset(dataset, args, env, run_dir, state_path):
    handler = load_dataset_handler(dataset, project_dir=str(PIPE))
    qc = QCExclusions(dataset, project_dir=PIPE)
    scans = [s for s in handler.discover_all_scans()
             if qc.get_status(s["subject"], s["session"], s["task"]) != "EXCLUDE"]
    expected_count = {"aachen_itbs": 60, "sudmex_tms": 133}[dataset]
    if len(scans) != expected_count:
        raise RuntimeError(f"QC selection changed for {dataset}: {len(scans)} != {expected_count}")
    write_json(run_dir / f"{dataset}_expected_scans.json", [
        {"scan_id": scan_id(s), "timeseries": str(series_path(handler, s))} for s in scans])
    subjects = []
    resumed = []
    for subject in sorted(set(s["subject"] for s in scans)):
        subject_scans = [s for s in scans if s["subject"] == subject]
        if all(series_path(handler, s).exists() for s in subject_scans):
            resumed.extend(validate_series(series_path(handler, s)) for s in subject_scans)
        else:
            subjects.append(subject)
    write_json(run_dir / f"{dataset}_resumed_scans.json", resumed)
    print(f"{dataset}: {len(resumed)} scans already verified; {len(subjects)} subjects remaining", flush=True)
    for start in range(0, len(subjects), args.batch_subjects):
        batch_subjects = subjects[start:start + args.batch_subjects]
        selected = [s for s in scans if s["subject"] in batch_subjects]
        if all(series_path(handler, s).exists() for s in selected):
            for scan in selected:
                validate_series(series_path(handler, scan))
            print(f"{dataset}: verified existing batch {batch_subjects}", flush=True)
            continue
        tag = dataset + "_" + "-".join(batch_subjects)
        batch_dir = run_dir / tag
        batch_dir.mkdir(exist_ok=True)
        # Conservative allowance for input, CONN data, cleaned copy and margin.
        input_bytes = sum(int(np.prod(nib.load(s["func_file"]).shape)) * 4 for s in selected)
        required = input_bytes * 5 + 30 * 2**30
        if shutil.disk_usage(ROOT).free < required:
            raise RuntimeError(f"Insufficient scratch space for {tag}: need {required / 2**30:.1f} GiB")
        command = [sys.executable, "-u", str(PIPE / "scripts/preprocessing/01a_prepare_conn_inputs.py"),
                   "--dataset", dataset, "--subjects", *batch_subjects,
                   "--n-workers", str(min(args.workers, 4)), "--skip-done"]
        run_command(command, batch_dir / "prepare.log", env, state_path, tag + ": prepare")
        manifest_file = PIPE / "data" / "conn_input" / dataset / "manifest.json"
        manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
        entries = manifest["scans"]
        if {e["scan_id"] for e in entries} != {scan_id(s) for s in selected}:
            raise RuntimeError(f"Incomplete or unexpected manifest: {tag}")
        write_json(batch_dir / "manifest.json", manifest)
        for entry in entries:
            for key in ("bold_trimmed", "T1w", "motion_txt"):
                below(Path(entry[key]), PIPE / "data" / "conn_input" / dataset)
        run_script = batch_dir / "conn_run.m"
        matlab_copy(PIPE / "scripts/preprocessing/step01b_conn_denoise.m",
                    run_script, min(args.workers, len(entries)), dataset)
        run_command(matlab_command(run_script), batch_dir / "conn.log", env,
                    state_path, tag + ": CONN")
        intermediate = PIPE / "results" / "timeseries_intermediate" / dataset
        shutil.copy2(intermediate / "conn_project" / f"conn_{dataset}.mat",
                     batch_dir / f"conn_{dataset}.mat")
        geometry = []
        for entry in entries:
            cleaned = intermediate / (entry["scan_id"] + "_conn_cleaned.nii")
            geometry.append(validate_cleaned(cleaned, entry))
        write_json(batch_dir / "geometry_validation.json", geometry)
        run_command([sys.executable, "-u", str(PIPE / "scripts/preprocessing/01c_parcellate_conn_output.py"),
                     "--dataset", dataset, "--n-workers", str(min(args.workers, 4)),
                     "--keep-cleaned", "--keep-inputs"], batch_dir / "parcellate.log",
                    env, state_path, tag + ": parcellation")
        checks = [validate_series(series_path(handler, e), e["n_volumes"]) for e in entries]
        write_json(batch_dir / "validation.json", {"completed_utc": stamp(), "scans": checks})
        clean_batch(dataset, entries, batch_dir / "scratch_cleanup.json")
        print(f"[{stamp()}] {tag}: {len(checks)} scans verified; scratch released", flush=True)
    return [series_path(handler, s) for s in scans]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", choices=ALLOWED, default=list(ALLOWED))
    parser.add_argument("--batch-subjects", type=int, default=4)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--matlab-workers", type=int, default=12)
    parser.add_argument("--through", choices=("timeseries", "phiid"), default="phiid")
    parser.add_argument("--run-dir", type=Path, default=ROOT / "logs/tms_pipeline_20261125/bulk")
    args = parser.parse_args()
    if min(args.batch_subjects, args.workers, args.matlab_workers) < 1:
        parser.error("Batch and worker counts must be positive")
    run_dir = below(args.run_dir, ROOT)
    run_dir.mkdir(parents=True, exist_ok=True)
    state_path = run_dir / "status.json"
    lock = run_dir / "runner.lock"
    with lock.open("x") as stream:
        stream.write(str(os.getpid()))
    env = os.environ.copy()
    env["MATLAB_PREFDIR"] = str(ROOT / "tmp/matlab_prefs")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[key] = "1"
    env["PYTHONUNBUFFERED"] = "1"
    try:
        expected = []
        for dataset in args.datasets:
            expected.extend(process_dataset(dataset, args, env, run_dir, state_path))
        if args.through == "phiid":
            # MATLAB production scripts scan this directory. Reject foreign input.
            for path in (PIPE / "results/timeseries").glob("*_timeseries*.npy"):
                if not any(path.name.startswith(f"dataset-{d}_") for d in ALLOWED):
                    raise RuntimeError(f"Frozen/non-TMS timeseries must not be processed: {path}")
            for stage, relative in (("hrf", "preprocessing/step02_hrf_deconvolution_matlab.m"),
                                    ("phiid", "analysis/step04_compute_phiid_matlab.m")):
                run_script = run_dir / (stage + "_run.m")
                matlab_copy(PIPE / "scripts" / relative, run_script, args.matlab_workers)
                run_command(matlab_command(run_script), run_dir / (stage + ".log"),
                            env, state_path, stage)
                if stage == "hrf":
                    for path in expected:
                        validate_series(path.with_name(path.name.replace("_timeseries.npy", "_timeseries_deconv.npy")),
                                        np.load(path, mmap_mode="r").shape[0])
            for path in expected:
                prefix = path.name.removesuffix("_timeseries.npy")
                valid = np.std(np.load(path.with_name(prefix + "_timeseries_deconv.npy")), axis=0) > 0
                valid_pairs = np.outer(valid, valid) & ~np.eye(454, dtype=bool)
                for suffix in ("synergy", "redundancy", "synergy_total", "phi", "phi_r"):
                    matrix = np.load(PIPE / "results/phiid_matrices" / (prefix + "_" + suffix + ".npy"))
                    if matrix.shape != (454, 454) or not np.isfinite(matrix[valid_pairs]).all():
                        raise RuntimeError(f"Invalid {suffix} matrix for {prefix}")
                    if not np.allclose(matrix, matrix.T, equal_nan=True):
                        raise RuntimeError(f"Asymmetric {suffix} matrix for {prefix}")
                if not (PIPE / "results/phiid_matrices" / (prefix + "_summary.yaml")).exists():
                    raise RuntimeError(f"Missing PhiID metadata: {prefix}")
        write_json(state_path, {"updated_utc": stamp(), "status": "complete",
                               "through": args.through, "verified_scans": len(expected),
                               "datasets": args.datasets, "group_comparisons_run": False})
        print(f"[{stamp()}] COMPLETE: {len(expected)} scans through {args.through}", flush=True)
    except Exception as exc:
        write_json(state_path, {"updated_utc": stamp(), "status": "failed", "error": str(exc)})
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
