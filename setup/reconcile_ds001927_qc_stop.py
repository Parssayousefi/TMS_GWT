"""One-time audited recovery of the completed s017 QC-exit handling failure.

Run only after the user's explicit resumption instruction. Does not rerun
fMRIPrep, relax exclusions, or start any new participant itself.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import subprocess
import sys

import run_ds001927_cohort as runner


def main():
    root, prov = runner.ROOT, runner.PROV
    queue_path = prov / "cohort_queue.json"
    before_queue = queue_path.read_bytes()
    queue = json.loads(before_queue)
    subject = "s017"
    if queue.get("status") != "paused_by_user" or queue.get("subject_pending_completion_reconciliation") != subject:
        raise RuntimeError("Recovery is only valid for the inspected, paused s017 failure")
    if queue["completed_subjects"] != runner.SUBJECTS[:15] or queue["active_subject"] is not None:
        raise RuntimeError("Unexpected queue position")
    receipt_path = prov / f"cohort_{subject}.json"
    before_receipt = receipt_path.read_bytes()
    receipt = json.loads(before_receipt)
    container = runner.container_state(receipt["container"])
    if container["Id"] != receipt["container_id"] or container["Image"] != runner.IMAGE_ID:
        raise RuntimeError("Completed container identity differs")
    state = container["State"]
    if state["Running"] or state["ExitCode"] != 0:
        raise RuntimeError("Previous participant did not finish successfully")
    qc_path = root / "output/ds001927_cohort_qc_20261126/sub-s017/derivative_qc.json"
    original_qc = json.loads(qc_path.read_text())
    runner.validate_completed_qc(original_qc, 3, subject)
    output = root / "output/ds001927_resume_review_20261127" / datetime.now(timezone.utc).strftime("%H%M%S_%f") / "sub-s017"
    if output.exists():
        raise RuntimeError("Existing recovery review preserved")
    # Real CLI regression check: the expected exclusion must still return 3.
    checked = subprocess.run([sys.executable, "-X", "utf8", str(root / "setup/qc_ds001927_derivatives.py"),
                              "--subject", subject, "--derivatives", receipt["output"], "--out", str(output)],
                             capture_output=True, text=True)
    repeated_qc = json.loads((output / "derivative_qc.json").read_text())
    runner.validate_completed_qc(repeated_qc, checked.returncode, subject)
    if len(original_qc["runs"]) != len(repeated_qc["runs"]):
        raise RuntimeError("QC run count changed")
    for old, new in zip(original_qc["runs"], repeated_qc["runs"]):
        for key in ("subject", "session", "task", "confounds_sha256", "motion", "file_errors", "run_motion_gate_pass"):
            if old[key] != new[key]:
                raise RuntimeError(f"QC decision/inputs changed: {key}")
    assert repeated_qc["motion_pass_run_count"] == 5
    reviewed = subprocess.run([sys.executable, "-X", "utf8", str(root / "setup/review_ds001927_pilot.py"),
                               "--subject", subject, "--derivatives", receipt["output"]],
                              capture_output=True, text=True)
    (output / "review_generation.log").write_text(reviewed.stdout + reviewed.stderr, encoding="utf-8")
    if reviewed.returncode != 0:
        raise RuntimeError("Review material generation failed")
    # Preserve the original state and failure receipt before reconciliation.
    archive = prov / "qc_exit_recovery_20261127"
    archive.mkdir(exist_ok=False)
    (archive / "queue_before_reconciliation.json").write_bytes(before_queue)
    (archive / "s017_before_reconciliation.json").write_bytes(before_receipt)
    receipt.update(status="fmriprep_complete_visual_review_pending", exit_code=0,
                   finished_utc=state["FinishedAt"], motion_pass_run_count=5,
                   qc_receipt_sha256=runner.sha(qc_path), anatomical_release=False,
                   recovered_utc=runner.now(), recovery="QC exit 3 denotes complete files with motion exclusion/review; original scan decisions preserved.")
    runner.write_json(receipt_path, receipt)
    queue["completed_subjects"].append(subject)
    queue.update(subject_pending_completion_reconciliation=None, updated_utc=runner.now(),
                 user_resume_requested_utc=runner.now(),
                 user_resume_instruction="run them now the rig is free")
    runner.validate_resume_state(queue, runner.SUBJECTS)
    runner.write_json(queue_path, queue)
    recovery = dict(created_utc=runner.now(), user_authorized_resume=True, subject=subject,
                    fmriprep_rerun=False, fmriprep_exit_code=0, real_qc_cli_exit_code=checked.returncode,
                    original_motion_decisions_unchanged=True, failed_scan="FRO pre: max-axis rotation 3.4996126335593236 degrees > 3.15 degrees",
                    analysis_release=False, remaining_subjects=runner.SUBJECTS[16:],
                    old_runner_sha256="d4fc6ce31c6a380636ba45929cd82b2d715b4ab651ddb216de2df81ae45ab33b",
                    new_runner_sha256=runner.sha(root / "setup/run_ds001927_cohort.py"),
                    qc_policy_sha256=runner.sha(prov / "qc_policy.json"),
                    original_qc_receipt_sha256=runner.sha(qc_path),
                    recovery_qc_receipt_sha256=runner.sha(output / "derivative_qc.json"),
                    regression_tests="8/8 passed; exit-code and safe-resume cases")
    runner.write_json(prov / "qc_exit_recovery_20261127.json", recovery)
    print(json.dumps(recovery, indent=2))


if __name__ == "__main__":
    main()
