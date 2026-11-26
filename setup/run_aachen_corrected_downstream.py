"""Run the unchanged Aachen downstream chain in bounded disk-space batches.

Each five-person batch waits for successful fMRIPrep and independent visual QA.
Only runtime worker counts change in logged MATLAB execution copies. Numerical
processing remains identical to the staged production scripts. HCP and all other
cohorts are excluded by the staged registry and explicit input checks.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import yaml

from audit_tms_hrf import audit as audit_hrf

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "reruns/aachen_skullstrip_skip_20261125/pipeline"
RUN = ROOT / "logs/aachen_corrected_20261125/downstream"
QUEUE = ROOT / "logs/aachen_corrected_20261125/queue_status.json"
QA = ROOT / "output/aachen_alignment_review/corrected_cohort"
PYTHON = Path(sys.executable)
MATLAB = Path("C:/Program Files/MATLAB/R2025a/bin/matlab.exe")
BATCHES = [["01", "03", "04", "05", "06"], ["07", "08", "09", "10", "11"], ["12", "13", "14", "15", "16"]]
SESSIONS = ["pre90", "post90", "pre120", "post120"]
MEASURES = ["synergy", "redundancy", "synergy_total", "phi", "phi_r"]


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def read_shared_json(path):
    for attempt in range(20):
        try:
            return json.loads(path.read_text(encoding="utf-8-sig"))
        except (json.JSONDecodeError, PermissionError):
            if attempt == 19:
                raise
            time.sleep(.25)


def confined(path, parent=STAGE):
    resolved = Path(path).resolve()
    assert resolved.is_relative_to(parent.resolve()) and resolved != parent.resolve(), str(path)
    return resolved


def prefixes(subjects):
    return [f"dataset-aachen_itbs_sub-{s}_ses-{v}_task-rest_run-1" for s in subjects for v in SESSIONS]


def validate_matrix_support(matrices, valid_rois):
    assert valid_rois.sum() >= 2, "Insufficient nonconstant deconvolved regions"
    expected = np.outer(valid_rois, valid_rois) & ~np.eye(len(valid_rois), dtype=bool)
    for name, matrix in matrices.items():
        assert matrix.shape == expected.shape and not np.isinf(matrix).any(), name
        assert np.array_equal(np.isfinite(matrix), expected), f"Incomplete or inconsistent finite edge support: {name}"
        assert np.allclose(matrix, matrix.T, equal_nan=True), f"Asymmetry: {name}"
    return expected


def check_support_guard():
    valid = np.array([True, True, False])
    good = np.full((3, 3), np.nan)
    good[0, 1] = good[1, 0] = .01
    validate_matrix_support({"test": good}, valid)
    for bad in [np.full((3, 3), np.nan), np.where(np.eye(3), np.nan, 0.)]:
        try:
            validate_matrix_support({"test": bad}, valid)
        except AssertionError:
            continue
        raise AssertionError("Finite-edge guard failed to reject corrupted support")


def validate_input_scope(subjects, require_deconvolved=False):
    folder = STAGE / "results/timeseries"
    expected = {p + "_timeseries.npy" for p in prefixes(subjects)}
    assert {p.name for p in folder.glob("dataset-*_timeseries.npy")} == expected, "Unexpected raw inputs for global MATLAB stage"
    expected_deconv = {p + "_timeseries_deconv.npy" for p in prefixes(subjects)}
    present_deconv = {p.name for p in folder.glob("dataset-*_timeseries_deconv.npy")}
    assert present_deconv <= expected_deconv, "Non-Aachen or extra deconvolved input refused"
    if require_deconvolved:
        assert present_deconv == expected_deconv, "Incomplete deconvolution cohort"


def validate_timeseries(subjects):
    rows = []
    for prefix in prefixes(subjects):
        path = STAGE / "results/timeseries" / (prefix + "_timeseries.npy")
        data = np.load(path)
        meta = yaml.safe_load(path.with_name(prefix + "_timeseries_metadata.yaml").read_text())
        assert data.ndim == 2 and data.shape[1] == 454 and data.shape[0] >= 240
        assert np.isfinite(data).all() and meta["n_timepoints"] == data.shape[0]
        assert meta["denoising_method"] == "CONN_V25.b_aCompCor" and meta["tr"] == 2.0
        assert meta["bandpass"] == [0.008, 0.09] and meta["global_signal_regression"] is False
        assert meta["art_spike_regressors"] is False and meta["dataset"] == "aachen_itbs"
        rows.append({"file": str(path), "shape": list(data.shape), "sha256": sha256(path),
                     "zero_variance_rois": int((data.std(axis=0) < 1e-10).sum())})
    return rows


def verify_staged_scripts():
    manifest = json.loads((STAGE.parent / "stage_manifest.json").read_text())
    for item in manifest["changes"]:
        assert sha256(STAGE / item["file"]) == item["staged_sha256"]
    for name in manifest["unchanged_processing_scripts_verified"]:
        assert sha256(STAGE / name) == sha256(ROOT / "pipeline" / name)
    registry = yaml.safe_load((STAGE / "configs/dataset_configs.yaml").read_text())
    assert set(registry["datasets"]) == {"aachen_itbs"}
    config = yaml.safe_load((STAGE / "configs/datasets/aachen_itbs.yaml").read_text())
    assert config["paths"]["local"]["fmriprep_dir"].endswith("fmriprep_skullstrip_skip_20261125")
    return {"scripts_verified": len(manifest["unchanged_processing_scripts_verified"]) + len(manifest["changes"]),
            "registry": "Aachen only", "numerical_processing_changed": False}


def wait_for_subjects(subjects, state):
    state.update(phase="waiting_for_anatomical_release", waiting_subjects=subjects, updated_utc=now())
    write_json(RUN / "status.json", state)
    while True:
        queue = read_shared_json(QUEUE)
        failures = [s for s in subjects if s in queue["finished"] and not queue["finished"][s]["success"]]
        if failures:
            raise RuntimeError("fMRIPrep failures: " + ", ".join(failures))
        path = QA / "subject_status.csv"
        if path.exists():
            try:
                qa = pd.read_csv(path, dtype={"subject": str}).set_index("subject")
            except (pd.errors.EmptyDataError, pd.errors.ParserError, PermissionError):
                time.sleep(10)
                continue
            if all(s in qa.index for s in subjects):
                flags = [s for s in subjects if qa.loc[s, "status"] == "flag"]
                if flags:
                    raise RuntimeError("Anatomical review flagged: " + ", ".join(flags))
                ready = all(qa.loc[s, "status"] == "pass" and
                            queue["finished"].get(s, {}).get("success") is True for s in subjects)
                if ready:
                    for subject in subjects:
                        review = read_shared_json(QA / f"sub-{subject}/visual_review.json")
                        for item in review["evidence"]:
                            assert sha256(Path(item["path"])) == item["sha256"], "Reviewed image changed"
                    return
        time.sleep(30)


def run_command(label, args, state):
    state.update(phase=label, updated_utc=now())
    write_json(RUN / "status.json", state)
    log = RUN / (label + ".log")
    started = now()
    print(f"{started} Starting {label}", flush=True)
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run([str(x) for x in args], cwd=STAGE, stdout=stream, stderr=subprocess.STDOUT,
                                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    receipt = {"label": label, "command": [str(x) for x in args], "started_utc": started,
               "finished_utc": now(), "exit_code": result.returncode, "log": str(log)}
    write_json(RUN / (label + "_receipt.json"), receipt)
    if result.returncode:
        raise RuntimeError(f"{label} failed with exit {result.returncode}; see {log}")


def matlab_copy(relative, label, workers):
    source = STAGE / relative
    content = source.read_bytes()
    old = re.search(rb"N_WORKERS = (\d+);", content)
    assert old
    new = b"N_WORKERS = " + str(workers).encode() + b";"
    modified = content[:old.start()] + new + content[old.end():]
    assert modified[:old.start()] + old.group() + modified[old.start() + len(new):] == content
    target = RUN / (label + ".m")
    target.write_bytes(modified)
    write_json(RUN / (label + "_source.json"), {"source": str(source), "source_sha256": sha256(source),
               "execution_copy": str(target), "execution_sha256": sha256(target),
               "runtime_workers": workers, "numerical_processing_changed": False})
    # parsave_hrf is an unchanged neighbouring helper in the production tree.
    return [MATLAB, "-nosplash", "-nodesktop", "-batch",
            f"addpath('{(STAGE / 'scripts/preprocessing').as_posix()}'); run('{target.as_posix()}');"]


def load_review_resume():
    """Resume only a QA hold before the next batch, never a partial stage."""
    state = read_shared_json(RUN / "status.json")
    assert state["phase"] == "failed" and state["error"].startswith("RuntimeError: Anatomical review flagged: ")
    assert Path(state["pipeline"]).resolve() == STAGE.resolve() and state["batches"] == BATCHES
    completed = state["completed_batches"]
    assert completed == list(range(1, len(completed) + 1)) and len(completed) < len(BATCHES)
    assert state["waiting_subjects"] == BATCHES[len(completed)]
    flagged = state["error"].split("Anatomical review flagged: ", 1)[1].split(", ")
    assert set(flagged) <= set(state["waiting_subjects"])
    qa = pd.read_csv(QA / "subject_status.csv", dtype={"subject": str}).set_index("subject")
    assert all(qa.loc[s, "status"] == "pass" for s in flagged), "Flagged participant has not been released"
    entries = []
    complete_subjects = []
    for number in completed:
        subjects = BATCHES[number - 1]
        complete_subjects.extend(subjects)
        for step in ["motion_qc", "prepare", "conn", "parcellate"]:
            receipt = read_shared_json(RUN / f"batch{number}_{step}_receipt.json")
            assert receipt["exit_code"] == 0
        saved = read_shared_json(RUN / f"batch{number}_timeseries_verification.json")
        assert validate_timeseries(subjects) == saved, "Completed batch inputs changed"
        manifest = read_shared_json(RUN / f"batch{number}_manifest.json")
        assert {x["scan_id"] for x in manifest["scans"]} == {p.removeprefix("dataset-aachen_itbs_") for p in prefixes(subjects)}
        assert len(manifest["scans"]) == 20
        entries.extend(manifest["scans"])
    validate_input_scope(complete_subjects)
    next_batch = len(completed) + 1
    assert not list(RUN.glob(f"batch{next_batch}_*_receipt.json")), "Next batch already started; manual recovery required"
    return state, entries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--resume-after-review", action="store_true")
    args = parser.parse_args()
    verification = verify_staged_scripts()
    check_support_guard()
    verification["finite_edge_guard_positive_and_negative_cases"] = "passed"
    resumed = load_review_resume() if args.resume_after_review else None
    if args.validate_only:
        if resumed:
            verification["resume_completed_batches_verified"] = resumed[0]["completed_batches"]
        print(json.dumps(verification))
        return
    RUN.mkdir(parents=True, exist_ok=True)
    if resumed:
        state, all_manifest_entries = resumed
        snapshot = RUN / ("status_review_hold_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ") + ".json")
        shutil.copy2(RUN / "status.json", snapshot)
        state.setdefault("review_resumes", []).append({"utc": now(), "prior_state": str(snapshot),
                    "prior_error": state.pop("error"), "completed_batches_revalidated": state["completed_batches"].copy()})
        state["verification"] = verification
    else:
        assert not (RUN / "status.json").exists(), "Inspect existing run before restarting; no automatic overwrite"
        state = {"started_utc": now(), "pipeline": str(STAGE), "batches": BATCHES, "completed_batches": [],
                 "verification": verification,
                 "batch_rationale": "Bound temporary uncompressed BOLD/CONN disk use; per-scan processing and settings unchanged."}
        all_manifest_entries = []
    try:
        for number, subjects in enumerate(BATCHES, 1):
            if number in state["completed_batches"]:
                continue
            wait_for_subjects(subjects, state)
            assert shutil.disk_usage(ROOT).free > 95 * 1024**3, "Insufficient scratch space for a 20-scan CONN batch"
            batch = f"batch{number}"
            run_command(batch + "_motion_qc", [PYTHON, "-X", "utf8", STAGE / "scripts/qc/fmriprep_motion_qc.py",
                        "--dataset", "aachen_itbs"], state)
            qc = pd.read_csv(STAGE / "results/qc/aachen_itbs/exclusions.csv", dtype={"subject": str})
            selected = qc.loc[qc.subject.isin(subjects)]
            assert len(selected) == 20 and selected.status.isin(["INCLUDE", "TRUNCATE"]).all(), "Changed QC cohort needs review"
            run_command(batch + "_prepare", [PYTHON, "-X", "utf8", STAGE / "scripts/preprocessing/01a_prepare_conn_inputs.py",
                        "--dataset", "aachen_itbs", "--subjects", *subjects, "--n-workers", "4"], state)
            manifest_path = STAGE / "data/conn_input/aachen_itbs/manifest.json"
            manifest = json.loads(manifest_path.read_text())
            assert len(manifest["scans"]) == 20
            expected_scan_ids = {p.removeprefix("dataset-aachen_itbs_") for p in prefixes(subjects)}
            assert {x["scan_id"] for x in manifest["scans"]} == expected_scan_ids
            # Preflight all recursive cleanup targets in unchanged CONN/parcellation.
            for scan in manifest["scans"]:
                confined(STAGE / "data/conn_input/aachen_itbs" / scan["scan_id"])
                confined(scan["bold_trimmed"])
                confined(scan["T1w"])
            confined(STAGE / "results/timeseries_intermediate/aachen_itbs/conn_project/conn_aachen_itbs")
            write_json(RUN / (batch + "_manifest.json"), manifest)
            queue = read_shared_json(QUEUE)
            workers = 4 if queue["status"] != "complete" else 12
            run_command(batch + "_conn", matlab_copy("scripts/preprocessing/step01b_conn_denoise.m",
                                                    batch + "_conn_runtime", workers), state)
            for scan in manifest["scans"]:
                cleaned = STAGE / "results/timeseries_intermediate/aachen_itbs" / (scan["scan_id"] + "_conn_cleaned.nii")
                assert cleaned.is_file() and cleaned.stat().st_size > 100000000, str(cleaned)
            run_command(batch + "_parcellate", [PYTHON, "-X", "utf8", STAGE / "scripts/preprocessing/01c_parcellate_conn_output.py",
                        "--dataset", "aachen_itbs", "--subjects", *subjects, "--n-workers", "4"], state)
            rows = validate_timeseries(subjects)
            write_json(RUN / (batch + "_timeseries_verification.json"), rows)
            all_manifest_entries.extend(manifest["scans"])
            state["completed_batches"].append(number)
            print(f"Verified {batch}: {len(rows)} parcellated scans", flush=True)
        # Wait for combined mapping/QA sign-off before final dataset-level stages.
        release_path = QA / "cohort_release_receipt.json"
        state.update(phase="waiting_for_cohort_release", updated_utc=now())
        write_json(RUN / "status.json", state)
        while not release_path.exists():
            time.sleep(30)
        release = read_shared_json(release_path)
        assert release["passed"] is True and release["participants_passed"] == 15
        assert release["functional_runs_audited"] == 60
        assert set(release["subjects"]) == set(sum(BATCHES, []))
        assert set(release["required_sessions"]) == set(SESSIONS)
        assert Path(release["derivative_root"]).resolve() == (ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125").resolve()
        for item in release["evidence"]:
            assert sha256(Path(item["path"])) == item["sha256"], "Cohort release evidence changed"
        manifest["scans"] = all_manifest_entries
        manifest["note_corrected_batch_execution"] = "Three archived 20-scan batches; complete inventory restored after parcellation. Temporary inputs were removed by unchanged production cleanup."
        write_json(STAGE / "data/conn_input/aachen_itbs/manifest.json", manifest)
        write_json(RUN / "complete_manifest.json", manifest)
        subjects = sum(BATCHES, [])
        before = validate_timeseries(subjects)
        write_json(RUN / "inputs.json", before)
        validate_input_scope(subjects)
        run_command("hrf", matlab_copy("scripts/preprocessing/step02_hrf_deconvolution_matlab.m", "hrf_runtime", 24), state)
        audits = [audit_hrf(Path(row["file"])) for row in before]
        assert len({r["deconv_sha256"] for r in audits}) == 60
        assert [r["source_sha256"] for r in audits] == [r["sha256"] for r in before]
        pd.DataFrame(audits).to_csv(RUN / "hrf_output_verification.csv", index=False)
        validate_input_scope(subjects, require_deconvolved=True)
        run_command("phiid", matlab_copy("scripts/analysis/step04_compute_phiid_matlab.m", "phiid_runtime", 24), state)
        matrix_checks = []
        for prefix in prefixes(subjects):
            folder = STAGE / "results/phiid_matrices"
            meta = yaml.safe_load((folder / (prefix + "_summary.yaml")).read_text())
            assert meta["method"] == "gaussian_matlab" and meta["tau"] == 1 and meta["n_rois"] == 454
            deconv = np.load(STAGE / "results/timeseries" / (prefix + "_timeseries_deconv.npy"))
            assert np.isfinite(deconv).all() and meta["n_timepoints"] == deconv.shape[0]
            valid_rois = deconv.std(axis=0, ddof=1) > 0  # identical validity rule to production PhiID
            assert meta["n_zero_variance_rois"] == int((~valid_rois).sum())
            identity = re.fullmatch(r"dataset-aachen_itbs_sub-(\d+)_ses-(\w+)_task-rest_run-1", prefix)
            assert meta["dataset"] == "aachen_itbs" and str(meta["subject"]).zfill(2) == identity[1]
            assert meta["session"] == identity[2] and meta["task"] == "rest" and int(meta["run"]) == 1
            matrices = {m: np.load(folder / (prefix + "_" + m + ".npy")) for m in MEASURES}
            validate_matrix_support(matrices, valid_rois)
            np.testing.assert_allclose(matrices["phi_r"], matrices["phi"] + matrices["redundancy"], atol=1e-10, equal_nan=True)
            upper = np.triu_indices(454, 1)
            for name, matrix in matrices.items():
                np.testing.assert_allclose(meta["mean_" + name], np.nanmean(matrix[upper]), atol=1e-12, rtol=1e-10)
            finite_phi_r = matrices["phi_r"][upper]
            finite_phi_r = finite_phi_r[np.isfinite(finite_phi_r)]
            assert bool(meta["positive_phi_r_check"]) == bool((finite_phi_r >= 0).all())
            matrix_checks.append({"scan": prefix, "valid_pairs": int(np.isfinite(matrices["phi"][np.triu_indices(454, 1)]).sum()),
                                  "files": {m: sha256(folder / (prefix + "_" + m + ".npy")) for m in MEASURES}})
        write_json(RUN / "phiid_output_verification.json", matrix_checks)
        state.update(phase="complete", updated_utc=now(), parcellated_scans=60, hrf_scans=60, phiid_scans=60,
                     anatomy_release=str(release_path), anatomy_release_sha256=sha256(release_path))
        write_json(RUN / "status.json", state)
        print("Corrected Aachen pipeline verified: 60 scans through PhiID.", flush=True)
    except Exception as exc:
        state.update(phase="failed", updated_utc=now(), error=f"{type(exc).__name__}: {exc}")
        write_json(RUN / "status.json", state)
        raise


if __name__ == "__main__":
    main()
