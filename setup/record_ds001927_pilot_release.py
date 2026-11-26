"""Record the completed, visually inspected two-pilot coverage investigation.

The visual statements below record the interactive review on 26 November 2026;
they are not an automated classifier or an approval of future participants.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PROV = ROOT / "PROVENANCE/ds001927_preprocessing_20261126"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, record):
    if path.exists():
        raise RuntimeError(f"Existing decision receipt preserved: {path}")
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def main():
    stamp = datetime.now(timezone.utc).isoformat()
    trace_path = ROOT / "output/ds001927_coverage_trace_20261126/coverage_trace_summary.json"
    trace = json.loads(trace_path.read_text())
    assert trace["run_count"] == 12 and trace["max_transform_parity_error"] == 0
    assert trace["below_half_mask_parcel_runs"] == trace["below_half_mask_parcel_runs_with_fov_above_95pct"] == 220
    assert trace["zero_variance_parcel_runs"] == 0
    evidence = [trace_path, ROOT / "output/ds001927_coverage_trace_20261126/raw_payload_audit.json",
                ROOT / "output/ds001927_coverage_trace_20261126/run_coverage_trace.csv",
                ROOT / "output/ds001927_coverage_trace_20261126/parcel_coverage_trace.csv",
                ROOT / "output/ds001927_pilot_qc_20261126/derivative_qc.json",
                ROOT / "output/ds001927_pilot_qc_s022_20261126/derivative_qc.json"]
    for subject in ("s001", "s022"):
        state = json.loads(subprocess.check_output(["docker", "inspect", f"ds001927_pilot_{subject}_20261126"], text=True))[0]["State"]
        assert state["Running"] is False and state["ExitCode"] == 0
        if subject == "s022":
            qc = json.loads((ROOT / "output/ds001927_pilot_qc_s022_20261126/derivative_qc.json").read_text())
            assert qc["motion_pass_run_count"] == 5 and qc["six_run_files_complete"]
            record = dict(subject="sub-s022", reviewed_utc=stamp, fmriprep_exit_code=0,
                          fmriprep_finished_utc=state["FinishedAt"],
                          native_extraction_visual_review="No gross cortical clipping on reviewed sagittal/coronal/axial sections; usable FRO/OCC T1s produce an anatomically plausible common reference.",
                          standard_alignment_visual_review="Five axial levels show plausible T1/template alignment and GM segmentation.",
                          functional_alignment_visual_review="All six references reviewed; gross alignment plausible; inferior frontal/temporal signal loss is already visible in raw BOLD.",
                          six_run_motion_qc="Five primary passes; CTR pre fails max-axis translation (3.1134874 mm > 3 mm). Both CTR runs fail the stricter mean-FD rule.",
                          frontal_control_primary_eligible=False,
                          preprocessing_pilot_pass=True, full_cohort_release=False, analysis_release=False,
                          note="Failure of the control baseline removes this participant from the main FRO-versus-CTR contrast. Other valid sessions remain separately available. No automatic motion salvage.")
            save(PROV / "pilot_s022_completion_review.json", record)
            receipt_path = PROV / "pilot_s022.json"
            receipt = json.loads(receipt_path.read_text(encoding="utf-8-sig"))
            receipt.update(status="fmriprep_complete_reviewed; one motion failure; no PID release", finished_utc=state["FinishedAt"], exit_code=0)
            receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        directory = ROOT / "output/ds001927_pilot_review_20261126" / f"sub-{subject}"
        evidence.extend(directory / name for name in ("native_extraction_and_segmentation.png", "standard_alignment.png", "six_run_bold_alignment_and_coverage.png", "review_summary.json"))
        evidence.append(ROOT / "output/ds001927_coverage_trace_20261126" / f"sub-{subject}_raw_to_preprocessed_coverage.png")
    plan_path = ROOT / "docs/DS001927_SITE_COMPARISON_PLAN_20261126.md"
    outputs = ROOT / "reruns/ds001927_20261126/pipeline/results"
    assert not list((outputs / "phiid_matrices").glob("*"))
    assert not list((outputs / "timeseries").glob("*"))
    save(PROV / "analysis_plan_v1.json", dict(created_utc=stamp,
         plan_path=plan_path.relative_to(ROOT).as_posix(), plan_sha256=sha(plan_path),
         user_selected_question="Does prefrontal stimulation change information processing differently from an active control site?",
         contrast="(FRO_post-FRO_pre)-(CTR_post-CTR_pre)", primary_endpoints=["phi_r", "six_mode_storage"],
         correction="two-sided participant joint-sign-flip maximum-|t| across two primary endpoints",
         pid_results_exist_at_decision=False, denoising_started_at_decision=False,
         prior_cohorts_results_known=True, pilot_anatomy_motion_coverage_known=True,
         external_preregistration=False, group_fd_covariate=False))
    save(PROV / "pilot_coverage_release.json", dict(created_utc=stamp,
         cohorts_reviewed=["sub-s001", "sub-s022"], cohort_fmriprep_release=True,
         pilot_denoising_for_signal_qc_allowed=True, analysis_release=False,
         interpretation="No evidence of the Aachen repeat-extraction failure in these pilots. Weak inferior signal is visible in the unprocessed data and persists after processing; acquisition-related dropout is the leading interpretation, not demonstrated causal proof of a specific artifact mechanism.",
         source_image_support="All 220 parcel/run cases below half mask overlap retain >95% geometric acquired-FOV support; median raw relative signal in those cases is 0.189974. Complete raw and processed 4D payloads are finite; finite noise/variance alone is not regional validity.",
         transform_reproduction="Saved transform chain reproduces all twelve standard-space references exactly (max voxel error 0). This checks transform application, not registration ground truth.",
         release_scope="Same fMRIPrep configuration may process remaining subjects. Each future subject still needs anatomical review; no group or PID outcome release. Common coverage set and denoised time-series QC remain required.",
         reviewed_file_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in evidence}))
    print("Recorded pre-outcome plan and preprocessing-only pilot release.")


if __name__ == "__main__":
    main()
