"""Independent checks of corrected target extraction and all planned inference.

Only the specified new target-flow run is written. Numerical checks use fixed
sentinel scans, all participant contrasts, all 256 sign patterns, and a bounded
24-case call to the original MATLAB PhiIDFull implementation.
"""
from argparse import ArgumentParser
from datetime import datetime, timezone
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

sys.dont_write_bytecode = True
import analyse_corrected_aachen_target_flow as analysis

ROOT = analysis.ROOT
METRICS = analysis.original.METRICS


def independent_holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adjusted = np.minimum(1, np.maximum.accumulate((len(p)-np.arange(len(p))) * p[order]))
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def assert_close(actual, expected, tolerance=2e-12):
    np.testing.assert_allclose(actual, expected, atol=tolerance, rtol=2e-10)


def verify_extraction(out, scans, inventory, weights):
    # Fixed before reading effects: two-session mask averaging, unequal lengths,
    # and the superficially placed right ROI are represented by these cases.
    cases = [("01", "pre90", "lag1_full"), ("04", "pre120", "lag2_full"),
             ("04", "post90", "lag1_equal_length"), ("13", "post120", "lag1_full")]
    errors, regression_errors = [], []
    for sid, session, variant in cases:
        subject = f"aachen_itbs_{sid}"
        scan = inventory.loc[inventory.subject_id.eq(subject) & inventory.session.eq(session)].iloc[0]
        ts = np.load(analysis.TS / f"{scan.prefix}_timeseries_deconv.npy")
        lag, trim = analysis.original.VARIANTS[variant]
        if trim:
            length = int(inventory.loc[inventory.subject_id.eq(subject), "expected_timepoints"].min())
            ts = ts[:length]
        wb = weights.loc[weights.subject_id.eq(subject)]
        seed_union = sorted(set(wb.loc[wb.weight > 0, "roi"].astype(int) - 1))
        remote = np.array(sorted(set(range(454)) - set(seed_union)))
        by_seed = {}
        for seed in seed_union:
            edges = np.column_stack([np.full(len(remote), seed), remote])
            sigma = analysis.original.tn.pg.lagged_sigma_from_timeseries(ts, lag, edges)
            values = analysis.original.directed(analysis.original.tn.pg.phiid_atoms(sigma))
            by_seed[seed] = {name: value.mean() for name, value in values.items()}
            for j in [0, len(remote)-1]:
                destination = remote[j]
                xp, yp, xf, yf = ts[:-lag, seed], ts[:-lag, destination], ts[lag:, seed], ts[lag:, destination]
                for target, own_past, other_past, expected in [
                    (yf, yp, xp, values["te_out"][j]), (xf, xp, yp, values["te_in"][j])]:
                    reduced = np.column_stack([np.ones(len(target)), own_past])
                    full = np.column_stack([reduced, other_past])
                    r0 = target - reduced @ np.linalg.lstsq(reduced, target, rcond=None)[0]
                    r1 = target - full @ np.linalg.lstsq(full, target, rcond=None)[0]
                    actual = .5 * np.log(np.dot(r0, r0) / np.dot(r1, r1))
                    regression_errors.append(abs(actual-expected))
                    assert_close(actual, expected, 2e-10)
        for side in ["l", "r"]:
            selected = scans.loc[scans.subject_id.eq(subject) & scans.session.eq(session)
                                 & scans.variant.eq(variant) & scans.side.eq(side)].iloc[0]
            w = wb.loc[wb.side.eq(side)].set_index("roi").weight
            for name in by_seed[seed_union[0]]:
                # Explicit average over destinations, then over seed parcels.
                actual = sum(w.loc[seed+1] * by_seed[seed][name] for seed in seed_union)
                errors.append(abs(actual-selected[name]))
                assert_close(actual, selected[name])
    return {"sentinel_scan_variants": len(cases), "independent_weighting_max_error": max(errors),
            "TE_regression_checks": len(regression_errors), "TE_regression_max_error": max(regression_errors)}


def verify_pairs_and_inference(out, scans, pairs, inference):
    names = METRICS + ["te_out", "te_in"] + analysis.original.tn.pg.ATOM_NAMES
    pair_plan = pd.read_csv(out / "input_visit_pairs.csv").set_index(["subject_id", "comparison"])
    for row in pairs.itertuples():
        plan = pair_plan.loc[(row.subject_id, row.comparison)]
        block = scans.loc[scans.subject_id.eq(row.subject_id) & scans.variant.eq(row.variant)]
        pre = block.loc[block.session.eq(plan.baseline_session)].set_index("side")
        post = block.loc[block.session.eq(plan.followup_session)].set_index("side")
        assert_close(row.delta_fd, plan.delta_fd)
        for metric in names:
            baseline, followup = pre.loc["l", metric], post.loc["l", metric]
            if row.scope == "left_minus_right":
                baseline -= pre.loc["r", metric]
                followup -= post.loc["r", metric]
            assert_close([baseline, followup, followup-baseline],
                         [getattr(row, metric+suffix) for suffix in ["_baseline", "_followup", "_delta"]])
    p_error = 0.
    for variant in analysis.original.VARIANTS:
        observed_rows = inference.loc[inference.variant.eq(variant)]
        expected_p = {}
        for scope in analysis.original.SCOPES:
            d = pairs.loc[pairs.variant.eq(variant) & pairs.scope.eq(scope)].groupby("subject_id").mean(numeric_only=True)
            assert len(d) == 8
            values = d[[metric+"_delta" for metric in METRICS]].to_numpy()
            signs = np.array(list(itertools.product([-1, 1], repeat=8)))
            distribution = np.abs(stats.ttest_1samp(signs[:, :, None] * values[None], 0, axis=1).statistic)
            test = stats.ttest_1samp(values, 0, axis=0)
            raw = (distribution >= abs(test.statistic)[None]-1e-10).mean(0)
            adjusted = (distribution.max(1)[:, None] >= abs(test.statistic)[None]-1e-10).mean(0)
            for index, metric in enumerate(METRICS):
                result = observed_rows.loc[observed_rows.scope.eq(scope) & observed_rows.metric.eq(metric)].iloc[0]
                v = values[:, index]
                ci = stats.t.interval(.95, 7, loc=v.mean(), scale=stats.sem(v))
                assert_close([result.estimate, result.ci_low, result.ci_high, result.t, result.p_t],
                             [v.mean(), ci[0], ci[1], test.statistic[index], test.pvalue[index]])
                assert_close([result.p_sign_flip, result.p_maxT5], [raw[index], adjusted[index]], 1e-14)
                expected_p[(scope, metric)] = raw[index]
                p_error = max(p_error, abs(result.p_maxT5-adjusted[index]))
        raw = [expected_p[(r.scope, r.metric)] for r in observed_rows.itertuples()]
        assert_close(observed_rows.p_holm10, independent_holm(raw), 1e-14)
    return {"paired_rows_reconstructed": len(pairs), "inference_rows_independently_verified": len(inference),
            "sign_patterns_per_scope_variant": 256, "max_exact_maxT_p_error": p_error}


def verify_sensitivities(out, pairs):
    primary = pairs.loc[pairs.variant.eq("lag1_full")]
    motion = pd.read_csv(out / "aachen_baseline_motion_sensitivity.csv")
    intensity = pd.read_csv(out / "aachen_intensity_descriptive.csv")
    loo = pd.read_csv(out / "leave_one_out.csv")
    if (len(motion), len(intensity), len(loo)) != (10, 20, 80):
        raise ValueError("A planned sensitivity was dropped")
    diagnostics = []
    for row in motion.itertuples():
        d = primary.loc[primary.scope.eq(row.scope)].groupby("subject_id").mean(numeric_only=True)
        baseline = d[row.metric+"_baseline"].to_numpy()
        y = d[row.metric+"_delta"].to_numpy()
        x = np.column_stack([np.ones(len(d)), baseline-baseline.mean(), d.delta_fd.to_numpy()])
        inverse = np.linalg.inv(x.T @ x)
        pinv = inverse @ x.T
        beta = pinv @ y
        residual = y - x @ beta
        leverage = np.einsum("ij,ji->i", x, pinv)
        covariance = (pinv * (residual/(1-leverage))**2) @ pinv.T
        se = np.sqrt(covariance[0, 0])
        ci = beta[0] + np.array([-1, 1]) * stats.t.ppf(.975, 5) * se
        p = 2 * stats.t.sf(abs(beta[0]/se), 5)
        assert_close([row.estimate, row.ci_low, row.ci_high, row.p], [beta[0], *ci, p], 2e-10)
        diagnostics.append({"scope": row.scope, "metric": row.metric, "n": 8, "residual_df": 5,
            "design_rank": int(np.linalg.matrix_rank(x)), "design_condition_number": float(np.linalg.cond(x)),
            "maximum_leverage": float(leverage.max()), "FD_reference": "delta usable mean FD = 0; original uncentered FD model retained"})
    assert_close(motion.p_holm10, independent_holm(motion.p), 1e-14)
    for row in intensity.itertuples():
        v = primary.loc[primary.scope.eq(row.scope) & primary.comparison.eq(row.comparison), row.metric+"_delta"]
        ci = stats.t.interval(.95, 7, loc=v.mean(), scale=stats.sem(v))
        assert_close([row.estimate, row.ci_low, row.ci_high, row.p_nominal], [v.mean(), *ci, stats.ttest_1samp(v, 0).pvalue])
    assert_close(intensity.p_holm20, independent_holm(intensity.p_nominal), 1e-14)
    for row in loo.itertuples():
        d = primary.loc[primary.scope.eq(row.scope)].groupby("subject_id").mean(numeric_only=True)
        v = d.drop(index=row.omitted)[row.metric+"_delta"]
        assert_close([row.estimate, row.p_nominal], [v.mean(), stats.ttest_1samp(v, 0).pvalue])
    pd.DataFrame(diagnostics).to_csv(out / "motion_model_diagnostics.csv", index=False)
    return {"motion_models_verified": 10, "intensity_descriptive_tests_verified": 20, "leave_one_out_fits_verified": 80}


def verify_matlab(out):
    executable = Path("C:/Program Files/MATLAB/R2025a/bin/matlab.exe")
    preferences = out / "matlab_preferences"
    preferences.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment["MATLAB_PREFDIR"] = str(preferences)
    command = [str(executable), "-nosplash", "-nodesktop", "-batch",
               "run('" + (out / "verify_matlab.m").as_posix().replace("'", "''") + "');"]
    started = datetime.now(timezone.utc).isoformat()
    with (out / "matlab_verification.log").open("w", encoding="utf-8") as log:
        process = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
            env=environment, creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0, timeout=300)
    analysis.write_json(out / "matlab_process_receipt.json", {"started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(), "command": command, "exit_code": process.returncode,
        "preferences_directory": str(preferences),
        "scope": "24 preselected Aachen pair cases only; no cohort/preprocessing driver"})
    if process.returncode:
        log_text = (out / "matlab_verification.log").read_text(encoding="utf-8", errors="replace")
        if "MathWorks Licensing Error 96" in log_text:
            return {"status": "unavailable", "reason": "MathWorks license server communication error 96",
                    "cases_computed": 0, "cases_prepared": 24, "atoms_per_case": 16,
                    "scope": "Fresh direct MATLAB atom check was optional if feasible in the original plan; corrected saved-production matrix agreement and independent Python checks remain completed."}
        raise RuntimeError("Bounded original-MATLAB verification failed; inspect matlab_verification.log")
    result = json.loads((out / "matlab_verification.json").read_text())
    if result["status"] != "passed" or result["cases"] != 24 or result["atoms_per_case"] != 16:
        raise ValueError("Incomplete original-MATLAB atom check")
    return result


def main():
    parser = ArgumentParser()
    parser.add_argument("--run-dir", type=Path, required=True)
    out = parser.parse_args().run_dir.resolve()
    if out.parent != (ROOT / "output/corrected_aachen_target_flow").resolve():
        raise ValueError("Verification must target a new corrected Aachen target-flow run")
    provenance = json.loads((out / "provenance.json").read_text())
    for name, digest in provenance["input_sha256"].items():
        if analysis.sha(ROOT / name) != digest:
            raise ValueError(f"Source changed: {name}")
    for name, digest in provenance["output_sha256"].items():
        if analysis.sha(out / name) != digest:
            raise ValueError(f"Computed output changed: {name}")
    scans = pd.read_csv(out / "scan_metrics.csv")
    pairs = pd.read_csv(out / "participant_visit_changes.csv")
    inference = pd.read_csv(out / "inference.csv")
    inventory = pd.read_csv(out / "input_scan_inventory.csv")
    weights = pd.read_csv(out / "anatomical_weights_used.csv")
    if (len(scans), len(pairs), len(inference), len(inventory), len(weights)) != (192, 96, 30, 32, 7264):
        raise ValueError("Unexpected cohort, side, variant or hypothesis counts")
    assert set(inventory.subject_id) == {f"aachen_itbs_{s}" for s in analysis.SUBJECTS}
    assert set(scans.side) == {"l", "r"} and set(pairs.scope) == set(analysis.original.SCOPES)
    assert_close(scans.net_te, scans.net_unique+scans.net_synergistic, 2e-15)
    for row in scans.itertuples():
        expected = inventory.loc[inventory.subject_id.eq(row.subject_id)]
        length = int(expected.expected_timepoints.min()) if row.variant == "lag1_equal_length" else int(expected.loc[expected.session.eq(row.session), "expected_timepoints"].iloc[0])
        if row.n_timepoints != length:
            raise ValueError("Equal-length sensitivity used the wrong time window")
    with threadpool_limits(limits=2):
        extraction = verify_extraction(out, scans, inventory, weights)
        inference_checks = verify_pairs_and_inference(out, scans, pairs, inference)
        sensitivity_checks = verify_sensitivities(out, pairs)
    matlab = verify_matlab(out)
    numerical = pd.read_csv(out / "numerical_scan_checks.csv")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if matlab["status"] == "passed" else "python_passed_matlab_unavailable",
        "required_numerical_checks_passed": True, "fresh_matlab_check_completed": matlab["status"] == "passed",
        "input_hashes_verified": len(provenance["input_sha256"]), "computed_output_hashes_verified": len(provenance["output_sha256"]),
        "participants": 8, "bilateral_sides": 16, "scans": 32, "source_branch": str(analysis.STAGE),
        "unchanged_metrics_and_test_families": True, **extraction, **inference_checks, **sensitivity_checks,
        "matlab": matlab, "max_saved_corrected_matrix_error": float(numerical.filter(regex="max_saved_").max().max()),
        "max_independent_CMI_error": float(numerical.max_cmi_error.max()),
        "equal_length_subjects_with_changed_length": int(inventory.groupby("subject_id").expected_timepoints.nunique().gt(1).sum()),
        "inference_limit": "Small anatomically available subset; no sham causal identification; complete planned sensitivity set retained.",
        "verification_script_sha256": analysis.sha(Path(__file__)), "analysis_provenance_sha256": analysis.sha(out / "provenance.json")}
    report["output_sha256"] = {p.name: analysis.sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name != "independent_verification.json"}
    analysis.write_json(out / "independent_verification.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
