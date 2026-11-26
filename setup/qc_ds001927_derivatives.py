"""Read-only, header-only completeness and declared motion gate for ds001927.

This is not an anatomical release, a full NIfTI integrity check, or denoising.
No volumes are censored, concatenated, rewritten or selected to salvage a run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DERIVATIVES = ROOT / "data/ds001927/derivatives/fmriprep_23_2_1"
DEFAULT_POLICY = ROOT / "PROVENANCE/ds001927_preprocessing_20261126/qc_policy.json"
STANDARD_SPACE = "MNI152NLin2009cAsym"
STANDARD_RESOLUTION_MM = 2.0
MOTION_COLUMNS = [f"{kind}_{axis}" for kind in ("trans", "rot") for axis in "xyz"]
NUMERIC_KEYS = (
    "tr_seconds", "expected_volumes", "initial_volumes_to_remove",
    "minimum_remaining_seconds", "max_translation_mm", "max_rotation_degrees",
    "max_mean_fd_mm", "high_fd_threshold_mm", "max_high_fd_fraction",
    "strict_sensitivity_max_mean_fd_mm", "standardized_dvars_flag",
)


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def json_safe(value):
    """Invalid headers still produce a receipt, with unavailable numbers null."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def validate_policy(policy):
    """Reject incomplete policies; thresholds never silently fall back."""
    for key in NUMERIC_KEYS:
        value = policy.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"Policy requires numeric {key}")
        if not np.isfinite(value) or value < 0:
            raise ValueError(f"Policy {key} must be finite and nonnegative")
    for key in ("expected_volumes", "initial_volumes_to_remove"):
        if int(policy[key]) != policy[key]:
            raise ValueError(f"Policy {key} must be an integer")
    if policy.get("dataset") != "ds001927":
        raise ValueError("Policy dataset must be ds001927")
    if policy.get("sessions") != ["ctr", "fro", "occ"]:
        raise ValueError("Policy sessions must be ctr, fro, occ in that order")
    if policy.get("tasks") != ["restpre", "restpost"]:
        raise ValueError("Policy tasks must be restpre, restpost in that order")
    if policy["expected_volumes"] != 600 or policy["tr_seconds"] != 1.25:
        raise ValueError("ds001927 acquisition requires 600 volumes and TR 1.25 s")
    if not 0 <= policy["initial_volumes_to_remove"] < policy["expected_volumes"]:
        raise ValueError("Initial trim must leave at least one original volume")
    if policy["minimum_remaining_seconds"] > 600 * 1.25:
        raise ValueError("Minimum retained duration exceeds the acquired record")
    if policy["max_high_fd_fraction"] > 1:
        raise ValueError("High-FD fraction must be between zero and one")
    if policy["strict_sensitivity_max_mean_fd_mm"] > policy["max_mean_fd_mm"]:
        raise ValueError("Sensitivity mean-FD threshold must be at least as strict")
    if policy.get("auto_salvage") is not False:
        raise ValueError("Policy must explicitly disable auto_salvage")
    if policy.get("delete_and_rejoin_scattered_frames") is not False:
        raise ValueError("Policy must explicitly disable deleting/rejoining frames")
    return policy


def load_policy(path):
    blob = Path(path).read_bytes()
    return validate_policy(json.loads(blob)), sha256_bytes(blob)


def normalize_subject(subject):
    subject = subject.removeprefix("sub-")
    if not re.fullmatch(r"s\d{3}", subject):
        raise ValueError("Subject must be sNNN or sub-sNNN")
    return f"sub-{subject}"


def assess_confounds(frame, policy):
    """Assess the uninterrupted suffix; return full original-index diagnostics."""
    validate_policy(policy)
    n = len(frame)
    errors, review, warnings = [], [], []
    if n != policy["expected_volumes"]:
        errors.append(f"confounds_rows_{n}_expected_{policy['expected_volumes']}")
    if frame.columns.duplicated().any():
        errors.append("duplicate_confounds_columns")
    required = MOTION_COLUMNS + ["framewise_displacement"]
    missing = [name for name in required if name not in frame]
    errors.extend(f"missing_required_column:{name}" for name in missing)

    def numeric(name):
        if name not in frame or name in frame.columns[frame.columns.duplicated()]:
            return np.full(n, np.nan)
        return pd.to_numeric(frame[name], errors="coerce").to_numpy(dtype=float)

    nss_columns = [name for name in frame if name.startswith("non_steady_state_outlier")]
    nss = np.zeros(n, dtype=bool)
    for name in nss_columns:
        values = numeric(name)
        if not (np.isfinite(values).all() and np.isin(values, [0, 1]).all()):
            errors.append(f"invalid_nonsteady_state_flags:{name}")
        nss |= np.isfinite(values) & (values != 0)
    prefix = 0
    while prefix < n and nss[prefix]:
        prefix += 1
    if nss[prefix:].any():
        review.append("nonsteady_state_flags_outside_leading_contiguous_prefix")
    trim = max(int(policy["initial_volumes_to_remove"]), prefix)
    indices = np.arange(n, dtype=int)
    retained = indices >= trim
    count = int(retained.sum())
    duration = count * policy["tr_seconds"]
    fd, dvars = numeric("framewise_displacement"), numeric("std_dvars")
    motion = np.column_stack([numeric(name) for name in MOTION_COLUMNS])
    finite_motion = np.isfinite(motion).all(axis=1) & np.isfinite(fd)
    if count and not finite_motion[retained].all():
        errors.append("nonfinite_required_motion_or_fd_after_trim")
    if count and (fd[retained] < 0).any():
        errors.append("negative_fd_after_trim")
    if "std_dvars" not in frame:
        warnings.append("std_dvars_missing_diagnostic_only")
    elif count and not np.isfinite(dvars[retained]).all():
        warnings.append("nonfinite_std_dvars_after_trim_diagnostic_only")

    relative = np.full((n, 6), np.nan)
    if count:
        relative[retained] = motion[retained] - motion[trim]
        relative[:, 3:] = np.degrees(relative[:, 3:])
    translations = np.full(n, np.nan)
    rotations = np.full(n, np.nan)
    if count:
        translations[retained] = np.max(np.abs(relative[retained, :3]), axis=1)
        rotations[retained] = np.max(np.abs(relative[retained, 3:]), axis=1)
    evaluable = count > 0 and finite_motion[retained].all() and not (fd[retained] < 0).any()
    high_fd = np.isfinite(fd) & (fd > policy["high_fd_threshold_mm"])
    high_dvars = np.isfinite(dvars) & (dvars > policy["standardized_dvars_flag"])
    maximum_translation = float(translations[retained].max()) if evaluable else None
    maximum_rotation = float(rotations[retained].max()) if evaluable else None
    mean_fd = float(fd[retained].mean()) if evaluable else None
    high_fraction = float(high_fd[retained].mean()) if evaluable else None

    # Tiny arithmetic tolerance only for inclusive aggregate boundaries, e.g.
    # radians -> degrees. Frame flags use the declared strict 'greater than'.
    def passes(value, limit):
        return None if value is None else bool(value <= limit + 1e-12)

    checks = {
        "duration": bool(duration + 1e-12 >= policy["minimum_remaining_seconds"]),
        "translation": passes(maximum_translation, policy["max_translation_mm"]),
        "rotation": passes(maximum_rotation, policy["max_rotation_degrees"]),
        "mean_fd": passes(mean_fd, policy["max_mean_fd_mm"]),
        "high_fd_fraction": passes(high_fraction, policy["max_high_fd_fraction"]),
    }
    motion_pass = not errors and not review and all(value is True for value in checks.values())
    strict_fd = passes(mean_fd, policy["strict_sensitivity_max_mean_fd_mm"])
    metrics = {
        "n_original_volumes": n, "leading_nonsteady_state_volumes": prefix,
        "nonsteady_state_flagged_indices": indices[nss].tolist(),
        "trim_first_volumes": trim, "first_retained_original_index": trim if count else None,
        "last_retained_original_index": n - 1 if count else None,
        "n_retained_volumes": count, "retained_duration_seconds": duration,
        "max_axis_translation_mm": maximum_translation,
        "max_axis_rotation_degrees": maximum_rotation,
        "mean_fd_mm": mean_fd, "fraction_fd_above_threshold_retained": high_fraction,
        "fd_above_threshold_count_retained": int(high_fd[retained].sum()),
        "fd_above_threshold_count_all": int(high_fd.sum()),
        "std_dvars_above_threshold_count_retained": int(high_dvars[retained].sum()),
        "std_dvars_above_threshold_count_all": int(high_dvars.sum()),
        "std_dvars_finite_count_retained": int(np.isfinite(dvars[retained]).sum()),
        "mean_std_dvars_retained": float(dvars[retained][np.isfinite(dvars[retained])].mean())
        if np.isfinite(dvars[retained]).any() else None,
        "motion_checks": checks, "motion_gate_pass": bool(motion_pass),
        "strict_mean_fd_sensitivity_pass": strict_fd,
        "strict_sensitivity_motion_gate_pass": bool(motion_pass and strict_fd),
        "input_errors": errors, "review_reasons": review, "warnings": warnings,
        "automatic_censoring_applied": False, "interior_frames_removed": 0,
    }
    flags = pd.DataFrame({
        "original_volume_index": indices, "original_volume_number_1based": indices + 1,
        "time_seconds": indices * policy["tr_seconds"],
        "retained_for_motion_assessment": retained,
        "in_declared_initial_trim": indices < policy["initial_volumes_to_remove"],
        "nonsteady_state_flag": nss, "framewise_displacement_mm": fd,
        "fd_above_declared_threshold": high_fd, "std_dvars": dvars,
        "std_dvars_above_declared_threshold": high_dvars,
        "max_axis_translation_relative_mm": translations,
        "max_axis_rotation_relative_degrees": rotations,
        "required_motion_and_fd_finite": finite_motion,
    })
    return metrics, flags


def entities(path):
    return dict(part.split("-", 1) for part in path.name.split("_") if "-" in part)


def read_header(path, kind, policy):
    """Do not access dataobj/get_fdata: reading 600-volume payloads is unnecessary."""
    before = path.stat()
    image = nib.load(path, mmap=True)
    shape = tuple(image.shape)
    affine = np.asarray(image.affine, dtype=float)
    zooms = tuple(map(float, image.header.get_zooms()))
    errors = []
    if len(shape) != (4 if kind == "bold" else 3):
        errors.append(f"{kind}_wrong_dimensions:{shape}")
    if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3, :3])) < 1e-8:
        errors.append(f"{kind}_invalid_affine")
    if not np.allclose(zooms[:3], STANDARD_RESOLUTION_MM, atol=1e-5, rtol=0):
        errors.append(f"{kind}_not_2mm_resolution")
    if not np.allclose(np.linalg.norm(affine[:3, :3], axis=0), 2, atol=1e-5, rtol=0):
        errors.append(f"{kind}_affine_not_2mm_resolution")
    units = image.header.get_xyzt_units()
    if units[0] != "mm":
        errors.append(f"{kind}_spatial_units_not_mm:{units[0]}")
    tr_seconds = None
    if kind == "bold" and len(shape) == 4:
        if shape[3] != policy["expected_volumes"]:
            errors.append(f"bold_volumes_{shape[3]}_expected_{policy['expected_volumes']}")
        multiplier = {"sec": 1, "msec": .001, "usec": .000001}.get(units[1])
        if multiplier is None:
            errors.append(f"bold_unknown_time_units:{units[1]}")
        else:
            tr_seconds = zooms[3] * multiplier
            if not np.isclose(tr_seconds, policy["tr_seconds"], atol=1e-6, rtol=0):
                errors.append(f"bold_wrong_tr_seconds:{tr_seconds}")
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        errors.append(f"{kind}_changed_during_header_read")
    return {
        "path": str(path), "shape": list(shape), "affine": affine.tolist(),
        "voxel_sizes_header": list(zooms), "units": list(units), "tr_seconds": tr_seconds,
        "file_size_bytes": after.st_size, "mtime_ns": after.st_mtime_ns,
        "header_sha256": sha256_bytes(image.header.binaryblock),
        "payload_read": False, "errors": errors,
    }


def same_grid(first, second):
    return (first["shape"][:3] == second["shape"][:3]
            and np.allclose(first["affine"], second["affine"], atol=1e-5, rtol=0))


def audit_subject(derivatives, subject, policy, out, policy_path=None, policy_sha256=None):
    validate_policy(policy)
    subject = normalize_subject(subject)
    derivatives, out = Path(derivatives).resolve(), Path(out).resolve()
    if out == derivatives or derivatives in out.parents:
        raise ValueError("QC output must be separate from read-only derivatives")
    out.mkdir(parents=True, exist_ok=True)
    (out / "frame_flags").mkdir(exist_ok=True)
    expected = [(session, task) for session in policy["sessions"] for task in policy["tasks"]]
    inventory = {key: {kind: [] for kind in ("bold", "mask", "confounds", "boldref")} for key in expected}
    unexpected = []
    for path in sorted((derivatives / subject).glob("ses-*/func/*")):
        if not path.is_file():
            continue
        ent = entities(path)
        key = (ent.get("ses"), ent.get("task"))
        image = path.name.endswith((".nii", ".nii.gz"))
        standard = ent.get("space") == STANDARD_SPACE and ent.get("res") == "2"
        kind = None
        if path.name.endswith("_desc-confounds_timeseries.tsv"):
            kind = "confounds"
        elif image and standard:
            if "_desc-preproc_bold.nii" in path.name:
                kind = "bold"
            elif "_desc-brain_mask.nii" in path.name:
                kind = "mask"
            elif "_boldref.nii" in path.name:
                kind = "boldref"
        if kind is not None:
            if key not in inventory or ent.get("sub") != subject[4:]:
                unexpected.append(str(path))
            else:
                inventory[key][kind].append(path)

    records, reference_grid = [], None
    for (session, task), candidates in inventory.items():
        record = {
            "subject": subject, "session": session, "task": task,
            "files": {k: [str(p) for p in v] for k, v in candidates.items()},
            "file_errors": [], "headers": {}, "motion": None,
            "files_complete": False, "run_motion_gate_pass": False,
            "analysis_release": False,
        }
        for kind in ("bold", "mask", "confounds"):
            if len(candidates[kind]) != 1:
                record["file_errors"].append(f"{kind}_count_{len(candidates[kind])}_expected_1")
        if len(candidates["boldref"]) > 1:
            record["file_errors"].append(f"boldref_count_{len(candidates['boldref'])}_expected_at_most_1")
        for kind in ("bold", "mask", "boldref"):
            if len(candidates[kind]) == 1:
                try:
                    header = read_header(candidates[kind][0], kind, policy)
                    record["headers"][kind] = header
                    record["file_errors"].extend(header["errors"])
                except (OSError, ValueError, EOFError, nib.filebasedimages.ImageFileError) as exc:
                    record["file_errors"].append(f"{kind}_unreadable_header:{type(exc).__name__}:{exc}")
        bold = record["headers"].get("bold")
        if bold is not None:
            for kind in ("mask", "boldref"):
                if kind in record["headers"] and not same_grid(bold, record["headers"][kind]):
                    record["file_errors"].append(f"{kind}_does_not_match_bold_grid")
            if not bold["errors"]:
                if reference_grid is None:
                    reference_grid = bold
                elif not same_grid(reference_grid, bold):
                    record["file_errors"].append("bold_grid_differs_from_other_runs")
        if len(candidates["confounds"]) == 1:
            path = candidates["confounds"][0]
            try:
                blob = path.read_bytes()
                columns = next(csv.reader(blob.decode("utf-8-sig").splitlines(), delimiter="\t"))
                if len(columns) != len(set(columns)):
                    raise ValueError("duplicate columns in source TSV")
                from io import BytesIO
                frame = pd.read_csv(BytesIO(blob), sep="\t", na_values=["n/a"])
                metrics, flags = assess_confounds(frame, policy)
                record["motion"] = metrics
                record["confounds_sha256"] = sha256_bytes(blob)
                if bold is not None and len(bold["shape"]) == 4 and len(frame) != bold["shape"][3]:
                    record["file_errors"].append("confounds_length_does_not_match_bold")
                name = f"{subject}_ses-{session}_task-{task}_frame_flags.csv"
                flags.to_csv(out / "frame_flags" / name, index=False)
                record["frame_flags"] = str(out / "frame_flags" / name)
                record["file_errors"].extend(metrics["input_errors"])
            except (OSError, ValueError, UnicodeError, StopIteration, pd.errors.ParserError) as exc:
                record["file_errors"].append(f"confounds_unreadable:{type(exc).__name__}:{exc}")
        record["files_complete"] = not record["file_errors"]
        if not record["files_complete"]:
            record["state"] = "incomplete_or_invalid_files"
        elif record["motion"]["review_reasons"]:
            record["state"] = "needs_review"
        elif record["motion"]["motion_gate_pass"]:
            record["state"] = "motion_pass_anatomy_pending"
            record["run_motion_gate_pass"] = True
        else:
            record["state"] = "motion_excluded"
        records.append(record)

    complete = not unexpected and all(r["files_complete"] for r in records)
    passed = complete and all(r["run_motion_gate_pass"] for r in records)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "dataset": "ds001927",
        "subject": subject, "derivatives": str(derivatives), "policy": policy,
        "policy_path": str(policy_path) if policy_path is not None else None,
        "policy_sha256": policy_sha256, "script_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "standard_space": STANDARD_SPACE, "standard_resolution_mm": STANDARD_RESOLUTION_MM,
        "expected_run_count": len(expected), "files_complete_run_count": sum(r["files_complete"] for r in records),
        "motion_pass_run_count": sum(r["run_motion_gate_pass"] for r in records),
        "six_run_files_complete": bool(complete), "six_run_motion_gate_pass": bool(passed),
        "unexpected_run_files": unexpected, "analysis_release": False,
        "anatomical_and_coverage_release": "pending_separate_visual_review",
        "fmriprep_process_exit_verified": False,
        "scope": "Headers and confounds only; no 4D payload integrity or anatomical alignment approval. A successful process exit still needs separate verification.",
        "automatic_censoring_applied": False,
        "state": "incomplete_or_invalid_files" if not complete else
                 ("motion_pass_anatomy_pending" if passed else "motion_exclusion_or_review"),
        "runs": records,
    }
    json_path = out / "derivative_qc.json"
    temporary = json_path.with_suffix(".json.tmp")
    summary = json_safe(summary)
    temporary.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(json_path)
    flat = []
    for record in records:
        row = {key: record[key] for key in ("subject", "session", "task", "state", "files_complete", "run_motion_gate_pass", "analysis_release")}
        row["file_errors"] = "; ".join(record["file_errors"])
        if record["motion"] is not None:
            row.update({key: value for key, value in record["motion"].items()
                        if not isinstance(value, (dict, list))})
            row["review_reasons"] = "; ".join(record["motion"]["review_reasons"])
            row["warnings"] = "; ".join(record["motion"]["warnings"])
        flat.append(row)
    pd.DataFrame(flat).to_csv(out / "run_qc.csv", index=False)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--derivatives", default=DEFAULT_DERIVATIVES, type=Path)
    parser.add_argument("--policy", default=DEFAULT_POLICY, type=Path)
    args = parser.parse_args()
    try:
        out = args.out.resolve()
        out.relative_to(ROOT / "output")
        policy, policy_hash = load_policy(args.policy)
        result = audit_subject(args.derivatives, args.subject, policy, out, args.policy.resolve(), policy_hash)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({key: result[key] for key in (
        "subject", "state", "six_run_files_complete", "motion_pass_run_count",
        "six_run_motion_gate_pass", "analysis_release", "policy_sha256")}, indent=2))
    return 0 if result["six_run_motion_gate_pass"] else (2 if not result["six_run_files_complete"] else 3)


if __name__ == "__main__":
    raise SystemExit(main())
