"""Bounded sequential fMRIPrep queue; no automatic anatomical/PID approval.

The two pilots must be reviewed, copied storage verified, and the analysis plan
recorded before dispatch. Existing containers and receipts are never replaced.
Scratch is preserved. Each completed subject receives QC tables and review
figures; inference remains gated on human/agent image inspection and signal QC.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROV = ROOT / "PROVENANCE/ds001927_preprocessing_20261126"
IMAGE = "nipreps/fmriprep:23.2.1"
IMAGE_ID = "sha256:9189837e68b0a0c912fad8f0948d74e6280a2d0617b234c8dee738127a035233"
SUBJECTS = [f"s{i:03d}" for i in range(1, 24) if i not in (1, 22)]


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def docker(*args):
    return subprocess.check_output(["docker", *args], text=True, encoding="utf-8").strip()


def container_state(name):
    return json.loads(docker("inspect", name))[0]


def build_command(subject, derivatives, work, settings):
    if subject not in SUBJECTS:
        raise ValueError("This queue processes only the 21 reviewed-design nonpilot participants")
    return ["docker", "run", "--detach", "--name", f"ds001927_cohort_{subject}_20261126",
            "-v", f"{(ROOT / 'data/ds001927').as_posix()}:/data:ro",
            "-v", f"{derivatives.as_posix()}:/out", "-v", f"{work.as_posix()}:/work",
            "-v", f"{settings.as_posix()}:/settings:ro",
            "-v", f"{(ROOT / 'license.txt').as_posix()}:/opt/freesurfer/license.txt:ro",
            "-v", f"{(ROOT / 'tmp/templateflow').as_posix()}:/templateflow",
            "-e", "TEMPLATEFLOW_HOME=/templateflow", IMAGE, "/data", "/out", "participant",
            "--participant-label", subject, "--bids-filter-file", "/settings/cohort_bids_filter.json",
            "-w", "/work", "--fs-license-file", "/opt/freesurfer/license.txt",
            "--output-spaces", "MNI152NLin2009cAsym:res-2", "--nprocs", "12",
            "--omp-nthreads", "4", "--mem", "49152", "--fs-no-reconall",
            "--skull-strip-t1w", "force", "--longitudinal", "--ignore", "fieldmaps", "slicetiming",
            "--random-seed", "20261126", "--stop-on-first-crash", "--notrack"]


def validate_receipts():
    release = json.loads((PROV / "pilot_coverage_release.json").read_text())
    storage = json.loads((PROV / "storage.json").read_text())
    plan = json.loads((PROV / "analysis_plan_v1.json").read_text())
    if release.get("cohort_fmriprep_release") is not True or release.get("analysis_release") is not False:
        raise RuntimeError("Explicit preprocessing-only pilot release is required")
    for relative, digest in release["reviewed_file_sha256"].items():
        if sha(ROOT / relative) != digest:
            raise RuntimeError(f"Pilot review evidence changed: {relative}")
    if storage.get("status") != "copied_and_all_hashes_verified":
        raise RuntimeError("Verified storage copy required")
    if sha(ROOT / plan["plan_path"]) != plan["plan_sha256"] or plan["pid_results_exist_at_decision"] is not False:
        raise RuntimeError("Recorded pre-outcome plan differs")
    derivatives, storage_root = Path(storage["derivatives"]), Path(storage["storage_root"])
    derivatives.resolve().relative_to(storage_root.resolve())
    if not (derivatives / "sourcedata/freesurfer/fsaverage").is_dir():
        raise RuntimeError("Verified fsaverage output setup missing")
    if docker("image", "inspect", IMAGE, "--format", "{{.Id}}") != IMAGE_ID:
        raise RuntimeError("Pinned image identity changed")
    return storage_root, derivatives


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--subjects", nargs="+", choices=SUBJECTS, default=SUBJECTS)
    args = ap.parse_args()
    if len(set(args.subjects)) != len(args.subjects):
        raise RuntimeError("Duplicate subjects refused")
    storage_root, derivatives = validate_receipts()
    settings = ROOT / "reruns/ds001927_20261126/settings"
    settings.mkdir(parents=True, exist_ok=True)
    filter_path = settings / "cohort_bids_filter.json"
    filt = {"t1w": {"datatype": "anat", "suffix": "T1w", "session": ["ctr", "fro", "occ"]},
            "bold": {"datatype": "func", "suffix": "bold", "session": ["ctr", "fro", "occ"],
                     "task": ["restpre", "restpost"]}}
    if filter_path.exists() and json.loads(filter_path.read_text()) != filt:
        raise RuntimeError("Different cohort filter preserved")
    if not filter_path.exists():
        write_json(filter_path, filt)
    state_path = PROV / "cohort_queue.json"
    if state_path.exists():
        raise RuntimeError("Queue already recorded; inspect rather than launch a duplicate")
    state = dict(created_utc=now(), status="prepared", selected_subjects=args.subjects,
                 completed_subjects=[], active_subject=None, image_id=IMAGE_ID,
                 storage_root=str(storage_root), derivatives=str(derivatives),
                 runner_sha256=sha(Path(__file__)), filter_sha256=sha(filter_path),
                 plan_receipt_sha256=sha(PROV / "analysis_plan_v1.json"),
                 pilot_release_sha256=sha(PROV / "pilot_coverage_release.json"),
                 no_downstream_analysis_started=True, analysis_release=False,
                 scratch_cleanup=False, maximum_simultaneous_subjects=1)
    if args.prepare_only:
        print(json.dumps(state, indent=2))
        return
    # Exclusive creation avoids duplicate orchestrators.
    with state_path.open("x", encoding="utf-8") as stream:
        json.dump(state, stream, indent=2)
    logs = ROOT / "logs/ds001927_cohort_20261126"
    logs.mkdir(parents=True, exist_ok=True)
    try:
        for subject in args.subjects:
            validate_receipts()
            if shutil.disk_usage(storage_root).free < 100 * 1024**3:
                raise RuntimeError("Storage reserve below 100 GiB; queue stopped without deleting files")
            active = [x for x in docker("ps", "--format", "{{.Names}}").splitlines()
                      if x.startswith(("ds001927_pilot_", "ds001927_cohort_"))]
            if active:
                raise RuntimeError(f"Another dataset preprocessing container is active: {active}")
            name = f"ds001927_cohort_{subject}_20261126"
            if name in docker("ps", "-a", "--format", "{{.Names}}").splitlines():
                raise RuntimeError(f"Existing container preserved: {name}")
            if (derivatives / f"sub-{subject}").exists() or (derivatives / f"sub-{subject}.html").exists():
                raise RuntimeError(f"Existing participant output requires review: {subject}")
            work = storage_root / "work" / f"sub-{subject}"
            work.resolve().relative_to(storage_root.resolve())
            work.mkdir(parents=True, exist_ok=False)
            command = build_command(subject, derivatives, work, settings)
            receipt_path = PROV / f"cohort_{subject}.json"
            if receipt_path.exists():
                raise RuntimeError(f"Existing subject receipt preserved: {subject}")
            receipt = dict(subject=subject, container=name, created_utc=now(), command=command,
                           image_id=IMAGE_ID, output=str(derivatives), work=str(work),
                           expected_bold_runs=6, analysis_release=False, status="prepared")
            write_json(receipt_path, receipt)
            container_id = subprocess.check_output(command, text=True).strip()
            receipt.update(container_id=container_id, status="running")
            write_json(receipt_path, receipt)
            state.update(status="running", active_subject=subject, updated_utc=now())
            write_json(state_path, state)
            print(f"Started {name}", flush=True)
            while True:
                current = container_state(name)["State"]
                if not current["Running"]:
                    break
                time.sleep(20)
            with (logs / f"{subject}_fmriprep.log").open("w", encoding="utf-8") as stream:
                subprocess.run(["docker", "logs", name], stdout=stream, stderr=subprocess.STDOUT, check=True)
            receipt.update(finished_utc=current["FinishedAt"], exit_code=current["ExitCode"])
            if current["ExitCode"] != 0:
                receipt["status"] = "fmriprep_failed"
                write_json(receipt_path, receipt)
                raise RuntimeError(f"fMRIPrep failed for {subject}; inspect saved logs")
            qc_out = ROOT / "output/ds001927_cohort_qc_20261126" / f"sub-{subject}"
            qc = subprocess.run([sys.executable, "-X", "utf8", str(ROOT / "setup/qc_ds001927_derivatives.py"),
                                 "--subject", subject, "--derivatives", str(derivatives), "--out", str(qc_out)],
                                capture_output=True, text=True)
            (logs / f"{subject}_qc.log").write_text(qc.stdout + qc.stderr, encoding="utf-8")
            if not (qc_out / "derivative_qc.json").exists():
                raise RuntimeError(f"QC did not produce a receipt: {subject}")
            qc_record = json.loads((qc_out / "derivative_qc.json").read_text())
            if qc.returncode not in (0, 1) or not qc_record["six_run_files_complete"]:
                raise RuntimeError(f"Incomplete subject derivatives: {subject}")
            review = subprocess.run([sys.executable, "-X", "utf8", str(ROOT / "setup/review_ds001927_pilot.py"),
                                     "--subject", subject, "--derivatives", str(derivatives)], capture_output=True, text=True)
            (logs / f"{subject}_review.log").write_text(review.stdout + review.stderr, encoding="utf-8")
            if review.returncode != 0:
                raise RuntimeError(f"Review material generation failed: {subject}")
            receipt.update(status="fmriprep_complete_visual_review_pending", motion_pass_run_count=qc_record["motion_pass_run_count"],
                           qc_receipt_sha256=sha(qc_out / "derivative_qc.json"), anatomical_release=False)
            write_json(receipt_path, receipt)
            state["completed_subjects"].append(subject)
            state.update(active_subject=None, updated_utc=now())
            write_json(state_path, state)
            print(f"Completed preprocessing and generated QC for {subject}; visual/analysis release pending", flush=True)
        state.update(status="fmriprep_queue_complete_visual_review_pending", updated_utc=now())
        write_json(state_path, state)
    except BaseException as error:
        state.update(status="stopped_for_review", error=f"{type(error).__name__}: {error}", updated_utc=now())
        write_json(state_path, state)
        raise


if __name__ == "__main__":
    main()
