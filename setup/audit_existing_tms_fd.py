"""Audit registered FD rules on the exact, already-analysed TMS time windows.

Read-only provenance audit: no inference, censoring, reprocessing or new salvage.
Aachen inputs are confined to the corrected skull-strip-skip branch.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from io import BytesIO
import hashlib
import json
from pathlib import Path
import re

import nibabel as nib
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
CORRECTED = ROOT / "reruns/aachen_skullstrip_skip_20261125/pipeline"
BRANCHES = ("registered", "strict")
MOTION = [f"{kind}_{axis}" for kind in ("trans", "rot") for axis in "xyz"]
ROUTES = {
    "aachen_itbs": {
        "pipeline": CORRECTED,
        "derivatives": ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125",
        "manifests": [ROOT / "logs/aachen_corrected_20261125/downstream/complete_manifest.json"],
        "logs": sorted((ROOT / "logs/aachen_corrected_20261125/downstream").glob("batch*_prepare.log")),
        "normalization": CORE / "normalization",
    },
    "sudmex_tms": {
        "pipeline": ROOT / "pipeline",
        "derivatives": ROOT / "data/sudmex_tms/derivatives/fmriprep",
        "manifests": sorted((ROOT / "logs/tms_pipeline_20261125/bulk").glob("sudmex_tms_*/manifest.json")),
        "logs": sorted((ROOT / "logs/tms_pipeline_20261125/bulk").glob("sudmex_tms_*/prepare.log")) +
                [ROOT / "logs/tms_pipeline_20261125/status_lock_recovery/prepare_before_resume.log"],
        "normalization": ROOT / "output/tms_numit/20261125T155252_384543Z",
    },
}


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def read_tracked(path, hashes):
    path = Path(path)
    blob = path.read_bytes()
    key = path.relative_to(ROOT).as_posix()
    value = digest(blob)
    if key in hashes and hashes[key] != value:
        raise ValueError(f"Source changed during audit: {path}")
    hashes[key] = value
    return blob


def read_json(path, hashes):
    return json.loads(read_tracked(path, hashes))


def read_yaml(path, hashes):
    return yaml.safe_load(read_tracked(path, hashes))


def thresholds_from_config(config):
    q = config["qc"]
    values = {
        "mean_fd": q["max_fd_mean"], "peak_fd": q["max_fd_peak"],
        "min_timepoints": q["min_timepoints"],
        "translation": q["motion_thresholds"]["max_translation_mm"],
        "rotation": q["motion_thresholds"]["max_rotation_deg"],
        "first": config["preprocessing"]["volumes_to_remove"],
        "last": config["preprocessing"]["volumes_to_remove_last"],
        "tr": config["preprocessing"]["tr"],
    }
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or
           not np.isfinite(x) or x < 0 for x in values.values()):
        raise ValueError("Invalid or missing numeric registered QC setting")
    if any(int(values[k]) != values[k] for k in ("min_timepoints", "first", "last")):
        raise ValueError("Timepoint and boundary counts must be integers")
    return values


def retained_bounds(qc_row, original_count, thresholds):
    """Mirror the documented production window semantics, not length inference."""
    if qc_row["status"] == "INCLUDE":
        start = int(thresholds["first"])
        end = original_count - int(thresholds["last"])
    elif qc_row["status"] == "TRUNCATE":
        raw = [qc_row.get("truncate_start"), qc_row.get("truncate_end")]
        if any(pd.isna(x) or not np.isfinite(float(x)) or float(x) != int(float(x)) for x in raw):
            raise ValueError("TRUNCATE lacks explicit integer source bounds")
        start = max(int(float(raw[0])), int(thresholds["first"]))
        end = min(int(float(raw[1])), original_count)
    else:
        raise ValueError(f"Previously ineligible QC status: {qc_row['status']}")
    if not 0 <= start < end <= original_count:
        raise ValueError(f"Invalid retained interval [{start},{end}) of {original_count}")
    return start, end


def parse_prepare_log(text):
    # Tie the trim instruction directly to the adjacent source filename. Scan
    # headers alone are unsafe in the original parallel preparation logs.
    records = []
    for match in re.finditer(r"func source:\s*([^\r\n]+)\r?\n\s*trim:\s*([^\r\n]+)", text):
        path, instruction = match.groups()
        if instruction.startswith("first "):
            numbers = re.findall(r"\d+", instruction)
            if len(numbers) != 2:
                continue
            rule = {"kind": "boundary", "first": int(numbers[0]), "last": int(numbers[1])}
        elif instruction.startswith("QC truncation to volumes "):
            numbers = re.findall(r"\d+", instruction)
            if len(numbers) != 2:
                continue
            rule = {"kind": "truncate", "start": int(numbers[0]), "end": int(numbers[1])}
        else:
            continue
        records.append({"path": Path(path.strip()), "rule": rule, "instruction": instruction})
    return records


def verify_log_rule(rule, status, start, end, original_count, thresholds):
    if status == "INCLUDE":
        expected = {"kind": "boundary", "first": int(thresholds["first"]), "last": int(thresholds["last"])}
        return rule == expected and (start, end) == (expected["first"], original_count-expected["last"])
    return (rule["kind"] == "truncate" and
            (max(rule["start"], int(thresholds["first"])), min(rule["end"], original_count)) == (start, end))


def motion_metrics(confounds, start, end, status, thresholds):
    required = MOTION + ["framewise_displacement"]
    missing = [name for name in required if name not in confounds]
    if missing:
        raise ValueError(f"Missing required confounds: {missing}")
    retained = confounds.iloc[start:end][required].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if not len(retained) or not np.isfinite(retained).all() or (retained[:, -1] < 0).any():
        raise ValueError("Missing, nonfinite or negative FD/motion in retained rows")
    motion, fd = retained[:, :6].copy(), retained[:, -1]
    # Preserve the legacy screen's reference: INCLUDE uses absolute fMRIPrep
    # parameters; already-salvaged TRUNCATE uses its first retained frame.
    if status == "TRUNCATE":
        motion -= motion[0]
    translation = float(np.abs(motion[:, :3]).max())
    rotation = float(np.degrees(np.abs(motion[:, 3:])).max())
    mean_fd, peak_fd = float(fd.mean()), float(fd.max())
    fraction = float(np.mean(fd > .5))
    displacement_pass = (translation <= thresholds["translation"] + 1e-12 and
                         rotation <= thresholds["rotation"] + 1e-12)
    reasons = []
    if not displacement_pass:
        reasons.append("legacy_displacement_limit_exceeded")
    if len(fd) < thresholds["min_timepoints"]:
        reasons.append("registered_min_timepoints_not_met")
    if mean_fd > thresholds["mean_fd"] + 1e-12:
        reasons.append("registered_mean_fd_exceeded")
    if peak_fd > thresholds["peak_fd"] + 1e-12:
        reasons.append("registered_peak_fd_exceeded")
    strict_reasons = reasons.copy()
    if mean_fd > .30 + 1e-12:
        strict_reasons.append("strict_mean_fd_exceeded")
    if fraction > .20 + 1e-12:
        strict_reasons.append("strict_high_fd_fraction_exceeded")
    return {
        "n_retained": len(fd), "retained_mean_fd": mean_fd, "retained_peak_fd": peak_fd,
        "retained_high_fd_count": int(np.sum(fd > .5)), "retained_high_fd_fraction": fraction,
        "retained_max_translation_mm": translation, "retained_max_rotation_deg": rotation,
        "legacy_displacement_reference": "absolute" if status == "INCLUDE" else "first_retained",
        "legacy_displacement_pass": bool(displacement_pass),
        "registered_pass": not reasons, "strict_pass": not strict_reasons,
        "registered_fail_reasons": ";".join(reasons), "strict_fail_reasons": ";".join(strict_reasons),
    }


def load_evidence(dataset, route, hashes):
    config_path = ROOT / "pipeline/configs/datasets" / f"{dataset}.yaml"
    config = read_yaml(config_path, hashes)
    thresholds = thresholds_from_config(config)
    if dataset == "aachen_itbs":
        staged = read_yaml(route["pipeline"] / "configs/datasets/aachen_itbs.yaml", hashes)
        if thresholds != thresholds_from_config(staged):
            raise ValueError("Corrected Aachen run configuration differs from registered QC settings")
    qc_path = route["pipeline"] / "results/qc" / dataset / "exclusions.csv"
    qc = pd.read_csv(BytesIO(read_tracked(qc_path, hashes)), dtype={"subject": str})
    manifests = {}
    for path in route["manifests"]:
        manifest = read_json(path, hashes)
        if manifest["dataset"] != dataset:
            raise ValueError(f"Wrong dataset manifest: {path}")
        for scan in manifest["scans"]:
            key = scan["scan_id"]
            if key in manifests and manifests[key]["entry"] != scan:
                raise ValueError(f"Conflicting archived CONN manifests: {key}")
            manifests[key] = {"entry": scan, "source": str(path)}
    current_path = route["pipeline"] / "data/conn_input" / dataset / "manifest.json"
    current = read_json(current_path, hashes)
    for scan in current["scans"]:
        if scan["scan_id"] not in manifests or manifests[scan["scan_id"]]["entry"] != scan:
            raise ValueError(f"Current CONN manifest contradicts archived execution: {scan['scan_id']}")
    logs = {}
    for path in route["logs"]:
        for item in parse_prepare_log(read_tracked(path, hashes).decode("utf-8", errors="replace")):
            item["log"] = str(path)
            key = str(item["path"].resolve()).lower()
            if key in logs and logs[key]["rule"] != item["rule"]:
                raise ValueError(f"Conflicting source-specific trim instructions: {item['path']}")
            logs[key] = item
    prepared_path = route["normalization"] / "cache/prepared.json"
    prepared = read_json(prepared_path, hashes)
    normalized_hashes = {Path(k).as_posix(): value for k, value in prepared["input_sha256"].items()}
    pool_path = route["normalization"] / "cache" / f"{dataset}_baseline_A_pool.npy"
    pool_blob = read_tracked(pool_path, hashes)
    pool_shape = list(np.load(BytesIO(pool_blob)).shape)
    if pool_shape[0] != prepared["baseline_pool_sizes"][dataset]:
        raise ValueError("Normalization pool receipt disagrees with saved pool")
    return {
        "config": config, "thresholds": thresholds, "qc": qc, "qc_path": qc_path,
        "manifests": manifests, "logs": logs, "prepared_path": prepared_path,
        "normalization_hashes": normalized_hashes, "pool_path": pool_path, "pool_shape": pool_shape,
    }


def audit_scan(row, evidence, route, hashes):
    base = {name: getattr(row, name) for name in ("dataset", "subject_id", "session", "arm", "prefix")}
    base.update(bounds_verified=False, legacy_displacement_pass=False, registered_pass=False, strict_pass=False,
                registered_fail_reasons="", strict_fail_reasons="", evidence_error="")
    try:
        subject = row.subject_id.rsplit("_", 1)[1]
        scan_id = row.prefix.removeprefix(f"dataset-{row.dataset}_")
        expected_id = f"sub-{subject}_ses-{row.session}_task-rest_run-1"
        if scan_id != expected_id:
            raise ValueError("Scan prefix does not match its dataset/subject/session")
        qc_rows = evidence["qc"].loc[evidence["qc"].subject.eq(subject) &
                                         evidence["qc"].session.eq(row.session) & evidence["qc"].task.eq("rest")]
        if len(qc_rows) != 1:
            raise ValueError(f"Expected one QC decision, found {len(qc_rows)}")
        qc_row = qc_rows.iloc[0]
        base["legacy_qc_status"] = qc_row.status
        folder = route["derivatives"] / f"sub-{subject}" / f"ses-{row.session}" / "func"
        confound_paths = list(folder.glob(f"sub-{subject}_ses-{row.session}_task-rest*_desc-confounds_timeseries.tsv"))
        if len(confound_paths) != 1:
            raise ValueError(f"Expected one exact confounds file, found {len(confound_paths)}")
        confound_path = confound_paths[0]
        blob = read_tracked(confound_path, hashes)
        confounds = pd.read_csv(BytesIO(blob), sep="\t", usecols=lambda name: name in MOTION + ["framewise_displacement"])
        n_original = len(confounds)
        start, end = retained_bounds(qc_row, n_original, evidence["thresholds"])
        base.update(confounds_path=str(confound_path), confounds_sha256=digest(blob),
                    confounds_start=start, confounds_end=end, original_confounds_rows=n_original)
        if scan_id not in evidence["manifests"]:
            raise ValueError("No archived execution manifest for this scan")
        manifest = evidence["manifests"][scan_id]
        entry = manifest["entry"]
        expected_input_dir = route["pipeline"] / "data/conn_input" / row.dataset / scan_id
        if Path(entry["bold_trimmed"]).parent.resolve() != expected_input_dir.resolve():
            raise ValueError("Manifest points outside the current/corrected CONN branch")
        bold_path = folder / f"sub-{subject}_ses-{row.session}_task-rest_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz"
        image = nib.load(bold_path, mmap=True)
        if len(image.shape) != 4 or image.shape[3] != n_original:
            raise ValueError("BOLD header and original confounds length disagree")
        source_key = str(bold_path.resolve()).lower()
        if source_key not in evidence["logs"]:
            raise ValueError("No source-specific preparation log proving trim instruction")
        log = evidence["logs"][source_key]
        if not verify_log_rule(log["rule"], qc_row.status, start, end, n_original, evidence["thresholds"]):
            raise ValueError("Preparation log contradicts reconstructed QC/config window")
        ts_base = route["pipeline"] / "results/timeseries"
        meta_path = ts_base / f"{row.prefix}_timeseries_metadata.yaml"
        hrf_path = ts_base / f"{row.prefix}_timeseries_deconv_metadata.yaml"
        meta, hrf = read_yaml(meta_path, hashes), read_yaml(hrf_path, hashes)
        for key, expected in {"subject": subject, "session": row.session, "task": "rest", "dataset": row.dataset}.items():
            if str(meta[key]) != expected:
                raise ValueError(f"Time-series metadata identity mismatch: {key}")
        shapes = {}
        for suffix in ("_timeseries.npy", "_timeseries_deconv.npy"):
            path = ts_base / f"{row.prefix}{suffix}"
            raw = read_tracked(path, hashes)
            shapes[suffix] = list(np.load(BytesIO(raw)).shape)
        length = end - start
        lengths = [entry["n_volumes"], meta["n_timepoints"], hrf["original_shape"][0], hrf["deconv_shape"][0],
                   row.n_timepoints, shapes["_timeseries.npy"][0], shapes["_timeseries_deconv.npy"][0]]
        if any(int(value) != length for value in lengths):
            raise ValueError(f"Retained bounds disagree with CONN/time-series/HRF lengths: {length}, {lengths}")
        if any(shapes[suffix][1:] != [454] for suffix in shapes):
            raise ValueError("Unexpected parcellated time-series shape")
        tr = evidence["thresholds"]["tr"]
        if any(abs(float(value)-tr) > 1e-7 for value in (entry["tr"], meta["tr"], hrf["tr"], image.header.get_zooms()[3])):
            raise ValueError("TR differs across source header/CONN/time-series/HRF metadata")
        expected_original = ts_base / f"{row.prefix}_timeseries.npy"
        expected_deconv = ts_base / f"{row.prefix}_timeseries_deconv.npy"
        if Path(hrf["input_file"]).resolve() != expected_original.resolve() or Path(hrf["output_file"]).resolve() != expected_deconv.resolve():
            raise ValueError("HRF sidecar points outside exact audited signal branch")
        metrics = motion_metrics(confounds, start, end, qc_row.status, evidence["thresholds"])
        if abs(metrics["retained_mean_fd"] - row.mean_fd) > 1e-10:
            raise ValueError("Retained FD mean differs from current analysis covariate")
        for key, column in (("retained_max_translation_mm", "usable_max_trans_mm"),
                            ("retained_max_rotation_deg", "usable_max_rot_deg"),
                            ("retained_mean_fd", "usable_mean_fd_mm")):
            if pd.isna(qc_row[column]) or abs(metrics[key] - float(qc_row[column])) > 1e-8:
                raise ValueError(f"Recomputed retained metric disagrees with legacy QC: {column}")
        base.update(metrics)
        deconv_rel = expected_deconv.relative_to(ROOT).as_posix()
        is_baseline = row.session in ("pre90", "pre120", "t0")
        pool_member = is_baseline and deconv_rel in evidence["normalization_hashes"]
        if pool_member and hashes[deconv_rel] != evidence["normalization_hashes"][deconv_rel]:
            raise ValueError("Normalization baseline input hash differs from audited signal")
        base.update(
            bounds_verified=True, bounds_basis="QC/config rule + source-specific preparation log + archived CONN manifest + original BOLD header + both signal arrays and metadata + HRF provenance",
            input_window_semantics="zero-based, end-exclusive, uninterrupted; no new salvage",
            manifest_path=manifest["source"], preparation_log=log["log"], logged_trim_instruction=log["instruction"],
            bold_path=str(bold_path), bold_header_sha256=digest(image.header.binaryblock),
            bold_file_size_bytes=bold_path.stat().st_size, original_timepoints=n_original,
            timeseries_metadata_path=str(meta_path), hrf_metadata_path=str(hrf_path),
            timeseries_path=str(expected_original), deconv_path=str(expected_deconv),
            current_analysis_mean_fd=float(row.mean_fd),
            normalization_baseline_pool_member=bool(pool_member),
            normalization_pool_source=str(evidence["pool_path"]), normalization_prepared_source=str(evidence["prepared_path"]),
        )
    except (OSError, ValueError, KeyError, IndexError) as exc:
        base.update(bounds_verified=False, registered_pass=False, strict_pass=False,
                    registered_fail_reasons="input_evidence_unverified", strict_fail_reasons="input_evidence_unverified",
                    evidence_error=f"{type(exc).__name__}: {exc}")
    return base


def build_pair_audit(scans, original_pairs):
    lookup = {(row.dataset, row.subject_id, row.session): row for row in scans.itertuples()}
    result, cohorts = [], {branch: {} for branch in BRANCHES}
    for branch in BRANCHES:
        for dataset, block in scans.groupby("dataset"):
            sessions = {"pre90", "post90", "pre120", "post120"} if dataset == "aachen_itbs" else {"t0", "t1"}
            subjects = []
            for subject, group in block.loc[block.session.isin(sessions)].groupby("subject_id"):
                if set(group.session) == sessions and group[f"{branch}_pass"].all():
                    subjects.append(subject)
            cohorts[branch][dataset] = sorted(subjects)
        for pair in original_pairs.itertuples():
            records = [lookup.get((pair.dataset, pair.subject_id, session)) for session in
                       (pair.baseline_session, pair.followup_session)]
            eligible = all(row is not None and getattr(row, f"{branch}_pass") for row in records)
            result.append({
                "branch": branch, "dataset": pair.dataset, "subject_id": pair.subject_id,
                "comparison": pair.comparison, "arm": pair.arm,
                "baseline_session": pair.baseline_session, "followup_session": pair.followup_session,
                "baseline_prefix": records[0].prefix if records[0] else None,
                "followup_prefix": records[1].prefix if records[1] else None,
                "eligible": bool(eligible),
                "eligible_primary_cohort": pair.subject_id in cohorts[branch][pair.dataset],
            })
    return pd.DataFrame(result), cohorts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "output/existing_tms_fd_audit_20261126")
    args = parser.parse_args()
    out = args.out.resolve(); out.relative_to(ROOT / "output"); out.mkdir(parents=True, exist_ok=True)
    hashes = {}
    scans = pd.read_csv(BytesIO(read_tracked(CORE / "raw_scan_metrics.csv", hashes)))
    pairs = pd.read_csv(BytesIO(read_tracked(CORE / "raw_paired_changes.csv", hashes)))
    if scans.groupby("dataset").size().to_dict() != {"aachen_itbs": 60, "sudmex_tms": 133} or scans.prefix.duplicated().any():
        raise ValueError("Unexpected scan inventory; this audit targets the released 193 scans")
    evidence = {ds: load_evidence(ds, route, hashes) for ds, route in ROUTES.items()}
    for pipeline in (ROOT / "pipeline", CORRECTED):
        for relative in ("scripts/preprocessing/01a_prepare_conn_inputs.py", "scripts/archive/01_extract_timeseries_luppi_method.py",
                         "scripts/qc/fmriprep_motion_qc.py"):
            read_tracked(pipeline / relative, hashes)
    for relative in ("setup/run_tms_numit_sensitivity.py", "setup/normalize_corrected_aachen.py", "setup/normalize_pid_processing_modes.py",
                     "output/pid_modes_normalized/20261126T131608_596458Z/input_verification.json"):
        read_tracked(ROOT / relative, hashes)
    rows = []
    for index, scan in enumerate(scans.itertuples(), 1):
        rows.append(audit_scan(scan, evidence[scan.dataset], ROUTES[scan.dataset], hashes))
        if index % 40 == 0:
            print(f"Exact-window evidence checked: {index}/193", flush=True)
    frame = pd.DataFrame(rows)
    pair_audit, cohorts = build_pair_audit(frame, pairs)
    frame.to_csv(out / "scan_audit.csv", index=False)
    pair_audit.to_csv(out / "pair_audit.csv", index=False)
    (out / "cohort_subjects.json").write_text(json.dumps(cohorts, indent=2) + "\n", encoding="utf-8")
    baselines = frame.loc[frame.session.isin(["pre90", "pre120", "t0"])].copy()
    baselines.to_csv(out / "baseline_pool_audit.csv", index=False)
    baselines.loc[~baselines.registered_pass | ~baselines.strict_pass].to_csv(out / "baseline_failures.csv", index=False)
    counts = []
    for branch in BRANCHES:
        for dataset, block in frame.groupby("dataset"):
            selected = cohorts[branch][dataset]
            arms = block.loc[block.subject_id.isin(selected)].drop_duplicates("subject_id").arm.value_counts().to_dict()
            base = baselines.loc[baselines.dataset.eq(dataset)]
            member = base.normalization_baseline_pool_member.fillna(False)
            counts.append({"branch": branch, "dataset": dataset, "scans_audited": len(block),
                           "scans_retained": int(block[f"{branch}_pass"].sum()), "primary_subjects": len(selected),
                           "primary_arm_counts": arms,
                           "baseline_pool_members": int(member.sum()),
                           "baseline_pool_members_failing_scan_qc": int((member & ~base[f"{branch}_pass"]).sum()),
                           "baseline_pool_members_outside_primary_cohort": int((member & ~base.subject_id.isin(selected)).sum())})
    changed = [path for path, value in hashes.items() if digest((ROOT / path).read_bytes()) != value]
    receipt = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "state": "complete" if frame.bounds_verified.all() and not changed else "unverified_inputs_fail_closed",
        "scope": "Exact-input audit only; no inference or preprocessing changes",
        "all_193_windows_verified": bool(frame.bounds_verified.all()), "changed_source_files": changed,
        "registered_thresholds": {ds: ev["thresholds"] for ds, ev in evidence.items()},
        "strict_sensitivity": {"mean_fd_max": .30, "fraction_fd_above_0_5_max": .20,
                               "retains_registered_peak_and_displacement_and_length_rules": True},
        "legacy_displacement_failures": int((~frame.legacy_displacement_pass).sum()),
        "new_salvage_or_censoring": False, "counts": counts,
        "bounds_index_convention": "zero-based, start inclusive, end exclusive",
        "normalization_pool_membership_evidence": "Prepared input SHA256 plus baseline-only selection code and previously independently reconstructed identical pools; current baseline hashes matched. Pools are not regenerated by this audit.",
        "normalization_pools": {ds: {"path": str(ev["pool_path"]), "shape": ev["pool_shape"]} for ds, ev in evidence.items()},
        "source_sha256": hashes,
        "output_sha256": {name: digest((out / name).read_bytes()) for name in
                          ("scan_audit.csv", "pair_audit.csv", "cohort_subjects.json", "baseline_pool_audit.csv", "baseline_failures.csv")},
        "script_sha256": digest(Path(__file__).read_bytes()),
    }
    (out / "audit_receipt.json").write_text(json.dumps(receipt, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"state": receipt["state"], "all_193_windows_verified": receipt["all_193_windows_verified"],
                      "legacy_displacement_failures": receipt["legacy_displacement_failures"], "counts": counts}, indent=2))
    if receipt["state"] != "complete":
        print(frame.loc[~frame.bounds_verified, ["prefix", "evidence_error"]].to_string(index=False))
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
