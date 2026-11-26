"""Recompute the existing Aachen target-flow plan from corrected inputs only.

Numerical definitions and inference are imported unchanged. This adapter owns
only its new output run, builds metadata from corrected QC, and never loads
historical Aachen signals, matrices, weights, or any HCP/psychedelic outcomes.
"""
from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.io import savemat
from threadpoolctl import threadpool_limits
import yaml

sys.dont_write_bytecode = True
import analyse_targeted_information_flow as original

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "reruns/aachen_skullstrip_skip_20261125/pipeline"
TS = STAGE / "results/timeseries"
PHI = STAGE / "results/phiid_matrices"
DOWNSTREAM = ROOT / "logs/aachen_corrected_20261125/downstream"
ANATOMY = ROOT / "output/aachen_alignment_review/corrected_cohort"
WEIGHTS = ANATOMY / "target_mapping/subject_region_weights.csv"
RELEASE_SHA = "d02401c42be2ea20c9e6de008a30c9665d92e7151467bb0aa4c873a8dfd81895"
WEIGHTS_SHA = "bd204ebcbb3e1897a4ed7b2f08cdd8c064118048081d022835d22b386970103f"
SUBJECTS = ["01", "03", "04", "05", "06", "07", "11", "13"]
SESSIONS = ["pre90", "post90", "pre120", "post120"]
DEFINITION = "deposited_individual_corrected"
HASHES = {}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def track(path, expected=None):
    path = Path(path).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError(f"Input outside the authorized project: {path}")
    if path.is_relative_to(ROOT / "pipeline/results"):
        raise ValueError(f"Historical/frozen result input forbidden: {path}")
    key = path.relative_to(ROOT).as_posix()
    digest = sha(path)
    if expected is not None and digest != expected:
        raise ValueError(f"Released input changed: {path}")
    if key in HASHES and HASHES[key] != digest:
        raise ValueError(f"Input changed during the analysis: {path}")
    HASHES[key] = digest
    return path


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def input_metadata(out):
    release = json.loads(track(ANATOMY / "cohort_release_receipt.json", RELEASE_SHA).read_text())
    if not (release["passed"] is True and release["participants_passed"] == 15
            and release["functional_runs_audited"] == 60):
        raise ValueError("Corrected anatomical release is incomplete")
    integrity = json.loads(track(ANATOMY / "release_integrity_verification.json").read_text())
    if integrity.get("passed") is not True or integrity["release_sha256"] != RELEASE_SHA:
        raise ValueError("Corrected anatomical evidence was not verified")
    state = json.loads(track(DOWNSTREAM / "status.json").read_text())
    if (state.get("phase") != "complete" or Path(state["pipeline"]).resolve() != STAGE.resolve()
            or [state[k] for k in ["parcellated_scans", "hrf_scans", "phiid_scans"]] != [60, 60, 60]
            or state["anatomy_release_sha256"] != RELEASE_SHA):
        raise ValueError("Corrected downstream processing is incomplete or mismatched")
    hrf = pd.read_csv(track(DOWNSTREAM / "hrf_output_verification.csv")).set_index("file")
    matrix_verification = json.loads(track(DOWNSTREAM / "phiid_output_verification.json").read_text())
    matrix_hashes = {item["scan"]: item["files"] for item in matrix_verification}
    if len(hrf) != 60 or len(matrix_hashes) != 60:
        raise ValueError("Expected 60 verified corrected runs")
    weights = pd.read_csv(track(WEIGHTS, WEIGHTS_SHA))
    coverage = pd.read_csv(track(ANATOMY / "target_mapping/parcel_coverage.csv"))
    weights["definition"] = DEFINITION
    expected_subjects = {f"aachen_itbs_{s}" for s in SUBJECTS}
    if set(weights.subject_id) != expected_subjects or len(weights) != 16 * 454:
        raise ValueError("The anatomically fixed eight-subject weight set changed")
    if len(coverage) != 16 or not coverage.eligible.all() or not coverage.cortical_fraction.gt(.1).all():
        raise ValueError("The unchanged anatomical eligibility rule is not satisfied")
    for (subject, side), block in weights.groupby(["subject_id", "side"]):
        if side not in {"l", "r"} or not np.array_equal(np.sort(block.roi), np.arange(1, 455)):
            raise ValueError(f"Incomplete seed weights: {subject} {side}")
        if not np.isfinite(block.weight).all() or block.weight.lt(0).any() or block.loc[block.roi > 400, "weight"].ne(0).any():
            raise ValueError("Target weights must be finite, nonnegative and cortical only")
        np.testing.assert_allclose(block.weight.sum(), 1, atol=1e-10)
    qc = pd.read_csv(track(STAGE / "results/qc/aachen_itbs/exclusions.csv"), dtype={"subject": str})
    qc.subject = qc.subject.str.zfill(2)
    if qc.duplicated(["subject", "session", "task"]).any():
        raise ValueError("Run selection would be ambiguous")
    qc = qc.loc[qc.subject.isin(SUBJECTS) & qc.session.isin(SESSIONS)]
    if len(qc) != 32 or not qc.status.isin(["INCLUDE", "TRUNCATE"]).all() or not qc.task.eq("rest").all():
        raise ValueError("The corrected eight-participant, four-session inventory is incomplete")
    scans, pairs = [], []
    for subject in SUBJECTS:
        for session in SESSIONS:
            row = qc.loc[qc.subject.eq(subject) & qc.session.eq(session)].iloc[0]
            stem = f"dataset-aachen_itbs_sub-{subject}_ses-{session}_task-rest_run-1"
            meta = yaml.safe_load(track(TS / f"{stem}_timeseries_metadata.yaml").read_text())
            deconv = yaml.safe_load(track(TS / f"{stem}_timeseries_deconv_metadata.yaml").read_text())
            if not (str(meta["subject"]).zfill(2) == subject and meta["session"] == session
                    and meta["dataset"] == "aachen_itbs" and meta["n_rois"] == 454 and meta["tr"] == 2):
                raise ValueError(f"Corrected metadata mismatch: {stem}")
            if Path(deconv["input_file"]).resolve() != (TS / f"{stem}_timeseries.npy").resolve():
                raise ValueError("Deconvolution references an uncorrected input")
            if Path(deconv["output_file"]).resolve() != (TS / f"{stem}_timeseries_deconv.npy").resolve():
                raise ValueError("Deconvolution output points outside the corrected branch")
            validated = hrf.loc[f"{stem}_timeseries.npy"]
            track(TS / f"{stem}_timeseries.npy", validated.source_sha256)
            track(TS / f"{stem}_timeseries_deconv.npy", validated.deconv_sha256)
            scans.append({"dataset": "aachen_itbs", "subject_id": f"aachen_itbs_{subject}", "subject": subject,
                "session": session, "arm": "iTBS", "prefix": stem, "mean_fd": float(row.usable_mean_fd_mm),
                "qc_status": row.status, "expected_timepoints": int(meta["n_timepoints"]),
                "deconv_shape": json.dumps(deconv["deconv_shape"])})
        subject_scans = {r["session"]: r for r in scans if r["subject"] == subject}
        for pre, post in [("pre90", "post90"), ("pre120", "post120")]:
            pairs.append({"dataset": "aachen_itbs", "subject_id": f"aachen_itbs_{subject}", "arm": "iTBS",
                "comparison": f"{post}_minus_{pre}", "baseline_session": pre, "followup_session": post,
                "delta_fd": subject_scans[post]["mean_fd"] - subject_scans[pre]["mean_fd"]})
    scans, pairs = pd.DataFrame(scans), pd.DataFrame(pairs)
    if not np.isfinite(scans.mean_fd).all():
        raise ValueError("Missing corrected motion covariate")
    scans.to_csv(out / "input_scan_inventory.csv", index=False)
    pairs.to_csv(out / "input_visit_pairs.csv", index=False)
    weights.to_csv(out / "anatomical_weights_used.csv", index=False)
    return scans, pairs, weights, matrix_hashes


def extract(out, scans, pairs, weights, matrix_hashes):
    records, checks, vectors, expected, matlab_meta = [], [], [], [], []
    for subject, wb in weights.groupby("subject_id"):
        block = scans.loc[scans.subject_id.eq(subject)]
        series = {r.session: np.load(track(TS / f"{r.prefix}_timeseries_deconv.npy")) for r in block.itertuples()}
        equal = min(map(len, series.values()))
        for r in block.itertuples():
            full = series[r.session]
            if full.shape != (r.expected_timepoints, 454) or not np.isfinite(full).all() or (full.std(axis=0) <= 0).any():
                raise ValueError(f"Invalid corrected regional signals: {r.prefix}")
            for variant, (lag, trim) in original.VARIANTS.items():
                ts = full[:equal] if trim else full
                sideweights = {side: b.sort_values("roi").weight.to_numpy() for side, b in wb.groupby("side")}
                seeds = np.unique(np.concatenate([np.flatnonzero(w > 0) for w in sideweights.values()]))
                remote = np.setdiff1d(np.arange(454), seeds)
                edges = np.array([(i, j) for i in seeds for j in remote], int)
                sigma = original.tn.pg.lagged_sigma_from_timeseries(ts, lag, edges)
                atoms = original.tn.pg.phiid_atoms(sigma)
                values = original.directed(atoms)
                if not all(np.isfinite(v).all() for v in values.values()):
                    raise ValueError("Nonfinite target information values")
                error = max(float(np.max(abs(values["te_out"] - original.cmi(sigma, 0, 3, 1)))),
                            float(np.max(abs(values["te_in"] - original.cmi(sigma, 1, 2, 0)))))
                if error >= 2e-10:
                    raise ValueError("Independent TE/CMI identity failed")
                check = {"prefix": r.prefix, "definition": DEFINITION, "variant": variant,
                    "n_pairs": len(edges), "n_seed_union": len(seeds), "n_remote_parcels": len(remote), "max_cmi_error": error}
                if variant == "lag1_full":
                    for name in original.OLD:
                        saved = np.load(track(PHI / f"{r.prefix}_{name}.npy", matrix_hashes[r.prefix][name]))
                        err = float(np.max(abs(atoms[name] - saved[edges[:, 0], edges[:, 1]])))
                        if err >= 2e-8:
                            raise ValueError(f"Corrected production matrix mismatch: {r.prefix} {name} {err}")
                        check[f"max_saved_{name}_error"] = err
                checks.append(check)
                for side, w in sideweights.items():
                    pair_weights = w[edges[:, 0]] / len(remote)
                    np.testing.assert_allclose(pair_weights.sum(), 1, atol=1e-10)
                    record = {"dataset": "aachen_itbs", "subject_id": subject, "session": r.session,
                        "arm": "iTBS", "definition": DEFINITION, "variant": variant, "side": side,
                        "mean_fd": r.mean_fd, "n_timepoints": len(ts), "n_seed_parcels": int((w > 0).sum()),
                        "n_remote_parcels": len(remote)}
                    record.update({name: float(pair_weights @ value) for name, value in values.items()})
                    records.append(record)
                # Same preselected participant and three edge positions as the original plan.
                if subject == "aachen_itbs_01" and variant != "lag1_equal_length":
                    for index in np.linspace(0, len(edges)-1, 3, dtype=int):
                        i, j = edges[index]
                        vectors.append(np.array([ts[:-lag, i], ts[:-lag, j], ts[lag:, i], ts[lag:, j]]))
                        expected.append([atoms[name][index] for name in original.tn.pg.ATOM_NAMES])
                        matlab_meta.append({"prefix": r.prefix, "lag": lag, "x_roi": int(i+1), "y_roi": int(j+1)})
        print("Extracted corrected target signals:", subject, flush=True)
    scan_metrics = pd.DataFrame(records)
    scan_metrics.to_csv(out / "scan_metrics.csv", index=False)
    pd.DataFrame(checks).to_csv(out / "numerical_scan_checks.csv", index=False)
    cells = np.empty((len(vectors), 1), dtype=object)
    for i, value in enumerate(vectors):
        cells[i, 0] = value
    savemat(out / "matlab_validation_input.mat", {"vectors": cells, "expected": np.asarray(expected),
        "atom_names": np.array(original.tn.pg.ATOM_NAMES, dtype=object)})
    pd.DataFrame(matlab_meta).to_csv(out / "matlab_validation_pairs.csv", index=False)
    return original.make_pairs(scan_metrics, pairs)


def main():
    out = ROOT / "output/corrected_aachen_target_flow" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True, exist_ok=False)
    for path in [Path(__file__), ROOT / "setup/verify_corrected_aachen_target_flow.py", Path(original.__file__),
                 ROOT / "setup/analyse_tms_familywise.py", Path(original.tn.__file__), Path(original.tn.pg.__file__),
                 ROOT / "PROVENANCE/TARGET_FLOW_WORKLOG_20261125.md", ROOT / "setup/verify_target_flow_matlab.m"]:
        track(path)
    for path in (ROOT / "pipeline/docs/luppi_supps/software/MMI_PhiID").rglob("*.m"):
        track(path)
    reuse = {name: hashlib.sha256(inspect.getsource(getattr(original, name)).encode()).hexdigest()
             for name in ["directed", "cmi", "numerical_unit_checks", "make_pairs", "aachen_stats", "sensitivities"]}
    specification = {"created_utc": stamp(), "source_timeseries": str(TS), "source_matrices_for_checks": str(PHI),
        "source_weights": str(WEIGHTS), "subjects": SUBJECTS, "required_sessions": SESSIONS,
        "metrics": original.METRICS, "scopes": original.SCOPES, "variants": original.VARIANTS,
        "primary": "lag1_full; average both post-minus-own-pre changes within each participant",
        "families": "Exact max-|t| across five metrics per scope; Holm across ten tests within each variant. Sensitivity variants are reported in full, not selected.",
        "sensitivities": "Lag2, minimum common length across all four within-participant sessions, baseline plus uncentered change in usable mean FD, all leave-one-participant-out fits, and intensity-specific descriptive tests (Holm20).",
        "seed_definition": "Unchanged deposited native masks; corrected anatomical weights. Same remote pool for both sides excludes the union of left/right seed parcels.",
        "right_scope": "Right values are retained descriptively; the original inferential comparator is left-minus-right, not an additional standalone right hypothesis.",
        "normalization": "Raw Gaussian MMI information in nats; no NuMIT prior or new normalization.",
        "causal_limit": "Eight anatomically available participants, no sham; predictive dependence is not identified direct neural causation. Right cortex is a spatial comparator.",
        "source_mask_subject_selection": "Anatomical availability and the unchanged >10% cortical overlap rule; all eight corrected bilateral participants retained before outcomes.",
        "unchanged_function_source_sha256": reuse, "historical_Aachen_result_inputs_allowed": False,
        "production_and_frozen_cohorts_modified": False}
    write_json(out / "analysis_specification.json", specification)
    with threadpool_limits(limits=2):
        unit = original.numerical_unit_checks()
        scans, pair_plan, weights, matrix_hashes = input_metadata(out)
        pairs = extract(out, scans, pair_plan, weights, matrix_hashes)
        if len(pairs) != 96 or not pairs.groupby(["variant", "scope", "subject_id"]).size().eq(2).all():
            raise ValueError("Participant-level two-visit pairing is incomplete")
        pairs.to_csv(out / "participant_visit_changes.csv", index=False)
        results = []
        for variant, block in pairs.groupby("variant"):
            result = original.aachen_stats(block)
            result["model"], result["dataset"], result["definition"], result["variant"] = "paired_mean", "aachen_itbs", DEFINITION, variant
            results.append(result)
        inference = pd.concat(results, ignore_index=True)
        if len(inference) != 30 or not inference.n.eq(8).all():
            raise ValueError("All planned hypotheses must use eight independent participants")
        inference.to_csv(out / "inference.csv", index=False)
        original.sensitivities(pairs, out)
    verification_source = (ROOT / "setup/verify_target_flow_matlab.m").read_text()
    start = verification_source.index("addpath(")
    script = "root='" + ROOT.as_posix().replace("'", "''") + "';\nfolder='" + out.as_posix().replace("'", "''") + "';\n" + verification_source[start:]
    (out / "verify_matlab.m").write_text(script, encoding="utf-8")
    provenance = {"created_utc": stamp(), "status": "computed_pending_independent_verification", "seed": original.SEED,
        "analysis_specification_created_before_extraction": True, "units": "nats; marginal 95% CIs are not multiplicity adjusted",
        "participants": 8, "bilateral_seed_sides": 16, "scans": 32, "paired_intensity_visits": 16,
        "scan_metric_rows": 192, "participant_visit_scope_variant_rows": 96, "inference_rows": 30,
        "exact_sign_patterns_per_test_family": 256, "primary_family_tests": 10, "matlab_cases": 24,
        "sign_flip_assumption": "Symmetric participant-mean change distributions; no sham randomization test.",
        "unchanged_function_source_sha256": reuse, "synthetic_validation": unit,
        "input_sha256": HASHES, "script_sha256": sha(Path(__file__)),
        "output_sha256": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()},
        "historical_Aachen_results_read": False, "frozen_HCP_or_psychedelic_results_read": False,
        "production_inputs_modified": False}
    write_json(out / "provenance.json", provenance)
    print("COMPUTED", out, flush=True)
    print(inference.loc[inference.variant.eq("lag1_full"), ["scope", "metric", "n", "estimate", "ci_low", "ci_high", "p_sign_flip", "p_maxT5", "p_holm10"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
