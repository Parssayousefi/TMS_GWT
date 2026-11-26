"""Verify the frozen corrected-cohort release and all linked hash records.

Reads anatomy, source masks, and already documented sampled BOLD inputs only as
bytes for integrity verification. It does not recompute functional statistics.
Writes a separate verification receipt; never changes the reviewed release.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COHORT = ROOT / "output/aachen_alignment_review/corrected_cohort"
RECEIPT = COHORT / "cohort_release_receipt.json"
SUBJECTS = [f"{i:02d}" for i in range(1, 17) if i != 2]
SESSIONS = {"pre90", "post90", "pre120", "post120"}


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    release = json.loads(RECEIPT.read_text())
    if (release.get("passed") is not True or release.get("subjects") != SUBJECTS
            or set(release.get("required_sessions", [])) != SESSIONS
            or release.get("participants_passed") != 15
            or release.get("functional_runs_audited") != 60):
        raise ValueError("Incomplete corrected-cohort release")
    if (release.get("functional_outcomes_read") is not False
            or release.get("temporal_signal_quality_assessed") is not False
            or release.get("source_mask_definitions_unchanged") is not True):
        raise ValueError("Release scope is inconsistent with the anatomical audit")

    pending = [RECEIPT]
    inspected_json = set()
    verified = {}

    def check(path, expected):
        path = Path(path)
        if not path.is_absolute():
            path = ROOT / path
        path = path.resolve()
        if path in verified:
            if verified[path] != expected:
                raise ValueError(f"Conflicting recorded hashes: {path}")
            return
        if digest(path) != expected:
            raise ValueError(f"Release evidence changed: {path}")
        verified[path] = expected
        if path.suffix.lower() == ".json":
            pending.append(path)

    def visit(value):
        if isinstance(value, dict):
            if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str):
                check(value["path"], value["sha256"])
            for key in ["input_sha256", "output_sha256"]:
                for path, expected in value.get(key, {}).items():
                    check(path, expected)
            for nested in value.values():
                visit(nested)
        elif isinstance(value, list):
            for nested in value:
                visit(nested)

    while pending:
        path = pending.pop()
        if path in inspected_json:
            continue
        inspected_json.add(path)
        document = json.loads(path.read_text())
        visit(document)
        if "audit_script_sha256" in document:
            check(ROOT / "setup/aachen_alignment_anatomical_review.py", document["audit_script_sha256"])
        if "script_sha256" in document:
            if path.name == "bilateral_target_planes_provenance.json":
                script = "aachen_alignment_target_planes.py"
            elif path.name == "reference_spotcheck_provenance.json":
                script = "aachen_alignment_reference_spotcheck.py"
            elif path.parent.name == "reference_intensity" and path.name == "sampled_bold_provenance.json":
                script = "aachen_alignment_sampled_bold.py"
            elif path.parent.name == "reference_intensity" and path.name == "provenance.json":
                script = "aachen_alignment_reference_intensity.py"
            else:
                raise ValueError(f"Unrecognized script provenance: {path}")
            check(ROOT / "setup" / script, document["script_sha256"])
        if "audit_provenance_sha256" in document:
            check(Path(document["independent_audit"]) / "provenance.json", document["audit_provenance_sha256"])

    statuses = pd.read_csv(COHORT / "subject_status.csv", dtype={"subject": str})
    if set(statuses.subject) != set(SUBJECTS) or not statuses.status.eq("pass").all():
        raise ValueError("Participant review status is incomplete")
    queue = json.loads((COHORT / "queue_final_snapshot.json").read_text())
    for subject in SUBJECTS:
        item = queue["finished"][subject]
        if not (item["success"] is True and item["state"]["ExitCode"] == 0
                and item["inventory"]["complete"] is True):
            raise ValueError(f"Processing success not established: {subject}")

    weights = pd.read_csv(COHORT / "target_mapping/subject_region_weights.csv")
    coverage = pd.read_csv(COHORT / "target_mapping/parcel_coverage.csv")
    groups = weights.groupby(["subject_id", "side"])
    if len(groups) != int(coverage.eligible.sum()) or not groups.size().eq(454).all():
        raise ValueError("Incomplete deposited-ROI weights")
    if (not np.isfinite(weights.weight).all() or (weights.weight < 0).any()
            or weights.loc[weights.roi > 400, "weight"].ne(0).any()):
        raise ValueError("Invalid cortical target weights")
    np.testing.assert_allclose(groups.weight.sum(), 1, atol=1e-10)

    result = {"created_utc": datetime.now(timezone.utc).isoformat(), "passed": True,
              "release_path": str(RECEIPT), "release_sha256": digest(RECEIPT),
              "verified_files": len(verified), "verified_bytes": sum(p.stat().st_size for p in verified),
              "inspected_json_records": len(inspected_json), "subjects": SUBJECTS,
              "target_weight_groups": len(groups), "target_weight_rows": len(weights),
              "functional_statistics_recomputed": False,
              "scope": "Recursive recorded-hash integrity, processing completion, review status and cortical weight invariants."}
    destination = COHORT / "release_integrity_verification.json"
    destination.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
