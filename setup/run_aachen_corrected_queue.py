"""Finite Aachen-only fMRIPrep queue after independent pilot release.

Never deletes or replaces containers, work directories, or old derivatives.
Completion here means process/output completion; cohort anatomical QA is separate.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
SUBJECTS = ["01", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]
SESSIONS = ["pre90", "post90", "pre120", "post120"]
DERIVATIVES = ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125"
LOGS = ROOT / "logs/aachen_corrected_20261125"


def command(args, check=True):
    return subprocess.run(args, capture_output=True, text=True, encoding="utf-8", errors="replace", check=check)


def inventory(subject):
    base = DERIVATIVES / f"sub-{subject}"
    required = []
    for session in SESSIONS:
        func = base / f"ses-{session}/func"
        stem = f"sub-{subject}_ses-{session}_task-rest"
        required.extend([func / f"{stem}_desc-confounds_timeseries.tsv"])
        for suffix in ["desc-preproc_bold.nii.gz", "desc-brain_mask.nii.gz", "boldref.nii.gz"]:
            required.append(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_{suffix}")
    anat = base / "ses-anat/anat"
    stem = f"sub-{subject}_ses-anat"
    required += [anat / f"{stem}_desc-brain_mask.nii.gz",
                 anat / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz",
                 anat / f"{stem}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5",
                 DERIVATIVES / f"sub-{subject}.html"]
    missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]
    return {"expected_files": len(required), "missing": missing, "complete": not missing}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=Path, required=True)
    parser.add_argument("--parallel", type=int, choices=[1, 2, 3, 4], default=4)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--restart-note")
    args = parser.parse_args()
    release = json.loads(args.release.read_text(encoding="utf-8-sig"))
    assert release.get("passed") is True, "Pilot anatomical release has not passed"
    # Receipt schema may include more evidence, but this queue needs explicit root/subject.
    assert str(release["subject"]).zfill(2) == "06"
    assert Path(release["derivative_root"]).resolve() == DERIVATIVES.resolve()
    assert set(release["observed_sessions"]) == set(SESSIONS)
    assert all(release["checks"].values()), "Pilot checks contain a failure"
    pilot = json.loads(command(["docker", "inspect", "aachen_skip_06_20261125", "--format", "{{json .State}}"] ).stdout)
    assert pilot["Status"] == "exited" and pilot["ExitCode"] == 0
    assert inventory("06")["complete"], inventory("06")
    if args.validate_only:
        print("Pilot release, successful exit, output inventory, and project paths verified.")
        return
    LOGS.mkdir(parents=True, exist_ok=True)
    status_path = LOGS / "queue_status.json"
    prior = json.loads(status_path.read_text()) if status_path.exists() else {}
    restarts = list(prior.get("supervisor_restarts", []))
    if prior:
        restarts.append({"utc": datetime.now(timezone.utc).isoformat(),
                         "note": args.restart_note or "Supervisor resumed from actual container states"})
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "subjects": SUBJECTS,
              "release": str(args.release), "release_sha256": hashlib.sha256(args.release.read_bytes()).hexdigest(),
              "derivatives": str(DERIVATIVES), "parallel": args.parallel,
              "resources_per_new_container": {"nprocs": 5, "omp_threads": 4, "memory_mb": 32768},
              "note": "Process/output completion only. Full cohort anatomical review required before interpreting downstream results."}
    record["started_utc"] = prior.get("started_utc", record["started_utc"])
    record["supervisor_restarts"] = restarts
    finished = {}
    attention = {}
    previous = None
    while True:
        names = set(command(["docker", "ps", "-a", "--format", "{{.Names}}"] ).stdout.splitlines())
        pending = []
        running = []
        for subject in SUBJECTS:
            name = f"aachen_skip_{subject}_20261125"
            if subject in finished:
                continue
            if name not in names:
                pending.append(subject)
                continue
            state = json.loads(command(["docker", "inspect", name, "--format", "{{json .State}}"] ).stdout)
            if state["Status"] in ["running", "created", "restarting"]:
                running.append(subject)
                if state["Status"] == "running":
                    recent = command(["docker", "logs", "--since", state["StartedAt"], "--tail", "250", name], check=False)
                    lines = (recent.stdout + recent.stderr).splitlines()
                    errors = [line for line in lines if "exception calling callback" in line or "Traceback (most recent call last)" in line or "PermissionError:" in line]
                    if errors:
                        attention[subject] = {"started_at": state["StartedAt"], "messages": errors[-5:]}
            else:
                check = inventory(subject)
                success = state["Status"] == "exited" and state["ExitCode"] == 0 and check["complete"]
                finished[subject] = {"success": success, "state": state, "inventory": check}
                logs = command(["docker", "logs", name], check=False)
                (LOGS / f"sub-{subject}_fmriprep.log").write_text(logs.stdout + logs.stderr, encoding="utf-8")
        while pending and len(running) < args.parallel:
            subject = pending.pop(0)
            start = command(["powershell.exe", "-NoProfile", "-File",
                             str(ROOT / "setup/run_aachen_skullstrip_pilot.ps1"),
                             "-Subject", subject, "-NProcs", "5", "-OmpThreads", "4", "-MemoryMB", "32768"], check=False)
            (LOGS / f"sub-{subject}_launch.log").write_text(start.stdout + start.stderr, encoding="utf-8")
            if start.returncode:
                raise RuntimeError(f"Could not launch sub-{subject}; see its launch log. Existing jobs remain running.")
            running.append(subject)
        record.update({"updated_utc": datetime.now(timezone.utc).isoformat(), "running": running,
                       "pending": pending, "finished": finished,
                       "attention": attention,
                       "status": "running" if running or pending else "complete"})
        temporary = status_path.with_name("queue_status.json.tmp")
        temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
        temporary.replace(status_path)
        compact = {"running": running, "pending": pending,
                   "successful": [s for s, r in finished.items() if r["success"]],
                   "failed": [s for s, r in finished.items() if not r["success"]],
                   "attention": attention}
        if compact != previous:
            print(json.dumps(compact), flush=True)
            previous = compact
        if not running and not pending:
            if any(not r["success"] for r in finished.values()):
                raise RuntimeError("One or more subjects failed; no downstream analysis was launched")
            break
        time.sleep(30)


if __name__ == "__main__":
    main()
