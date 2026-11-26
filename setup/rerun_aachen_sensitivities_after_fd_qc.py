"""Refit cached corrected-Aachen representation/network checks after FD QC.

All four scans are required before averaging both intensity visits. No imaging
is reprocessed, normalization reference recalibrated, or HRF association tested.
These are repeated method/robustness checks, not independent treatment evidence.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

sys.dont_write_bytecode = True
import analyse_tms_workspace_scope as base

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z/sensitivities"
AUDIT = ROOT / "output/existing_tms_fd_audit_20261126"
OUT = ROOT / "output/tms_fd_reanalysis_20261126/sensitivities"
HASHES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def track(path, expected=None):
    path = Path(path).resolve()
    key = path.relative_to(ROOT).as_posix()
    value = sha(path)
    if expected is not None and value != expected:
        raise ValueError(f"Released input changed: {key}")
    if key in HASHES and HASHES[key] != value:
        raise ValueError(f"Input changed during this rerun: {key}")
    HASHES[key] = value
    return path


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def compare(actual, expected, keys, label):
    a = actual.set_index(keys).sort_index()
    b = expected.set_index(keys).sort_index()
    if not a.index.equals(b.index) or not a.index.is_unique or not b.index.is_unique:
        raise ValueError(f"Original-result keys differ: {label}")
    errors = []
    for column in a.columns.intersection(b.columns):
        if pd.api.types.is_numeric_dtype(a[column]) and pd.api.types.is_numeric_dtype(b[column]):
            x, y = a[column].to_numpy(float), b[column].to_numpy(float)
            np.testing.assert_allclose(x, y, rtol=2e-7, atol=2e-10, equal_nan=True,
                                       err_msg=f"Original-result parity: {label}/{column}")
            good = np.isfinite(x) & np.isfinite(y)
            errors.append(float(np.max(abs(x[good] - y[good]))) if good.any() else 0.0)
        elif not a[column].fillna("<NA>").equals(b[column].fillna("<NA>")):
            raise ValueError(f"Original-result column differs: {label}/{column}")
    return {"table": label, "rows": len(a), "maximum_absolute_error": max(errors, default=0.0)}


def pair(scans, metrics):
    rows = []
    for (branch, subject), block in scans.groupby(["branch", "subject_id"]):
        if len(block) != 4 or set(block.session) != {"pre90", "post90", "pre120", "post120"}:
            raise ValueError("Incomplete Aachen participant in cached sensitivity rerun")
        b = block.set_index("session")
        for intensity in (90, 120):
            pre, post = b.loc[f"pre{intensity}"], b.loc[f"post{intensity}"]
            row = {"branch": branch, "subject_id": subject, "intensity": intensity,
                   "delta_fd": post.mean_fd - pre.mean_fd}
            for metric in metrics:
                row.update({metric: post[metric] - pre[metric], metric + "_baseline": pre[metric]})
            rows.append(row)
    return pd.DataFrame(rows)


def summarize(data, family):
    if not np.isfinite(data).all().all():
        raise ValueError(f"Non-finite paired endpoint in {family}")
    p, corrected = base.exact_max_t(data.to_numpy())
    rows = [{"metric": metric, "family": family, **base.tms.paired_summary(data[metric]),
             "p_sign_flip": p[index], "p_maxT_family": corrected[index]} for index, metric in enumerate(data)]
    result = pd.DataFrame(rows)
    result["p_holm_family"] = base.holm(result.p)
    return result


def independent_exact_check(values, result):
    """Small chunks enumerate signed samples using SciPy's independent t test."""
    values = np.asarray(values, float)
    n = len(values)
    observed = stats.ttest_1samp(values, 0, axis=0).statistic
    direct = np.zeros(values.shape[1], int)
    corrected = np.zeros(values.shape[1], int)
    for start in range(0, 2 ** n, 256):
        signs = 2 * ((np.arange(start, min(start + 256, 2 ** n))[:, None] >> np.arange(n)) & 1) - 1
        t = abs(stats.ttest_1samp(signs[:, :, None] * values[None, :, :], 0, axis=1).statistic)
        direct += (t >= abs(observed)[None, :] - 1e-10).sum(axis=0)
        corrected += (t.max(axis=1)[:, None] >= abs(observed)[None, :] - 1e-10).sum(axis=0)
    p, adjusted = direct / (2 ** n), corrected / (2 ** n)
    np.testing.assert_allclose(result.p_sign_flip, p, atol=1e-12)
    np.testing.assert_allclose(result.p_maxT_family, adjusted, atol=1e-12)
    return {"participants": n, "sign_patterns": 2 ** n, "endpoints": values.shape[1],
            "maximum_adjusted_p_error": float(np.max(abs(result.p_maxT_family - adjusted)))}


def run(signals, profiles, order, out, independently_check=False):
    out.mkdir(parents=True, exist_ok=True)
    signal_metrics = [column for column in signals if "__" in column]
    profile_metrics = [column for column in profiles if "__" in column]
    sp, pp = pair(signals, signal_metrics), pair(profiles, profile_metrics)
    sp.to_csv(out / "signal_participant_changes.csv", index=False)
    pp.to_csv(out / "network_participant_changes.csv", index=False)
    signal_tests, network_tests, visit_tests, adjusted, influence, nonparam, diagnostics = [], [], [], [], [], [], []
    for branch, block in sp.groupby("branch"):
        pooled = block.groupby("subject_id").mean(numeric_only=True)
        columns = [scope + "__" + metric for scope in ("outside_workspace", "within_workspace", "difference")
                   for metric in ("tdmi", "self_sum", "redundancy")]
        result = summarize(pooled[columns], "temporal_components9")
        result["branch"] = branch
        signal_tests.append(result)
    variants = {}
    for branch, block in sp.groupby("branch"):
        pooled = block.groupby("subject_id").mean(numeric_only=True)
        for metric in ("whole_brain__phi", "within_workspace__phi", "outside_workspace__phi", "difference__phi"):
            variants[branch + "::" + metric] = pooled[metric]
    variant_frame = pd.DataFrame(variants)
    estimator = summarize(variant_frame, "estimator_variants12")
    for branch, block in pp.groupby("branch"):
        pooled = block.groupby("subject_id").mean(numeric_only=True)
        n = len(pooled)
        for family, columns in [("network_blocks", [metric for metric in profile_metrics if metric.startswith("block__")]),
                                ("network_omissions8", [metric for metric in profile_metrics if metric.startswith("without__")])]:
            result = summarize(pooled[columns], family)
            result["branch"] = branch
            network_tests.append(result)
        visit = {}
        for metric in ("outside_workspace__phi", "difference__phi"):
            wide = block.pivot(index="subject_id", columns="intensity", values=metric).sort_index()
            first = np.array([wide.loc[subject, order[subject]] for subject in wide.index])
            second = np.array([wide.loc[subject, 210 - order[subject]] for subject in wide.index])
            for name, values in [("90", wide[90]), ("120", wide[120]), ("120_minus_90", wide[120] - wide[90]),
                                 ("first_visit", first), ("second_visit", second), ("second_minus_first", second - first)]:
                visit[metric + "::" + name] = np.asarray(values)
            values = pooled[metric].to_numpy()
            for index, subject in enumerate(pooled.index):
                influence.append({"branch": branch, "metric": metric, "omitted_subject": subject,
                                  **base.tms.paired_summary(np.delete(values, index))})
            rng = np.random.default_rng(2026112528)
            bootstrap = values[rng.integers(n, size=(10000, n))].mean(axis=1)
            low, high = np.quantile(bootstrap, [.025, .975])
            nonparam.append({"branch": branch, "metric": metric, "positive_changes": int((values > 0).sum()), "n": n,
                             "wilcoxon_p": stats.wilcoxon(values, method="exact").pvalue,
                             "sign_test_p": stats.binomtest(int((values > 0).sum()), n, .5).pvalue,
                             "bootstrap_low": low, "bootstrap_high": high,
                             "visit_correlation": wide[90].corr(wide[120])})
            baseline = pooled[metric + "_baseline"]
            frame = pd.DataFrame({"y": values, "intercept": 1.0, "fd": pooled.delta_fd,
                                  "baseline_centered": baseline - baseline.mean()})
            for name, predictors in [("motion", ["intercept", "fd"]),
                                     ("baseline", ["intercept", "baseline_centered"]),
                                     ("baseline_motion", ["intercept", "baseline_centered", "fd"])]:
                design = frame[predictors].to_numpy()
                rank = np.linalg.matrix_rank(design)
                if rank != len(predictors) or n - rank < 2:
                    raise ValueError("Insufficient design rank or residual degrees of freedom")
                leverage = np.diag(design @ np.linalg.pinv(design))
                if leverage.max() >= 1 - 1e-9:
                    raise ValueError("Unit leverage in a sensitivity model")
                diagnostics.append({"branch": branch, "metric": metric, "model": name, "n": n,
                                    "rank": int(rank), "residual_df": int(n - rank), "max_leverage": float(leverage.max())})
                adjusted.append({"branch": branch, "metric": metric, "model": name,
                                 **base.hc3(frame, "y", predictors, "intercept")})
        result = summarize(pd.DataFrame(visit), "visit_and_intensity12")
        result["branch"] = branch
        visit_tests.append(result)
    means = sp.groupby(["branch", "subject_id"]).mean(numeric_only=True)
    method = {}
    for scope in ("whole_brain", "within_workspace", "outside_workspace", "difference"):
        for measure in ("phi", "phi_r"):
            metric = scope + "__" + measure
            method[metric] = means.xs("deconvolved")[metric] - means.xs("no_deconvolution")[metric]
    method_frame = pd.DataFrame(method)
    method_results = summarize(method_frame, "deconvolution_interactions8")
    adjusted = pd.DataFrame(adjusted)
    adjusted["p_holm6"] = adjusted.groupby("branch").p.transform(base.holm)
    tables = {"estimator_sensitivity_results.csv": estimator, "temporal_component_results.csv": pd.concat(signal_tests, ignore_index=True),
              "network_results.csv": pd.concat(network_tests, ignore_index=True), "visit_consistency_results.csv": pd.concat(visit_tests, ignore_index=True),
              "baseline_motion_sensitivity.csv": adjusted, "participant_influence.csv": pd.DataFrame(influence),
              "nonparametric_checks.csv": pd.DataFrame(nonparam), "method_interactions.csv": method_results,
              "model_diagnostics.csv": pd.DataFrame(diagnostics)}
    for filename, table in tables.items():
        table.to_csv(out / filename, index=False)
    if independently_check:
        save_json(out / "independent_signflip_verification.json", {
            "passed": True, "estimator12": independent_exact_check(variant_frame, estimator),
            "method8": independent_exact_check(method_frame, method_results)})
    return tables, sp, pp


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "completion_receipt.json").exists():
        raise FileExistsError("This completed sensitivity rerun is immutable")
    original_verification = json.loads(track(SOURCE / "followup_verification.json").read_text())
    if original_verification["status"] != "passed":
        raise ValueError("Original cached follow-up is unverified")
    for path in (SOURCE / "corrected_follow.py", SOURCE / "corrected_follow_verify.py",
                 ROOT / "setup/deepen_aachen_spatial_effect.py", ROOT / "setup/verify_aachen_spatial_followup.py",
                 ROOT / "setup/analyse_tms_familywise.py", ROOT / "setup/analyse_tms_all_post_average.py"):
        track(path)
    signals = pd.read_csv(track(SOURCE / "signal_decomposition_scans.csv"))
    profiles = pd.read_csv(track(SOURCE / "network_profile_scans.csv"))
    if len(signals) != 180 or len(profiles) != 120 or signals.subject_id.nunique() != 15:
        raise ValueError("Unexpected original corrected sensitivity inventory")
    weights = pd.read_csv(track(SOURCE / "network_pair_counts.csv")).query("edges > 0")
    computed = sum(profiles["block__" + row.network_pair] * row.weight_of_outside_mean for row in weights.itertuples())
    np.testing.assert_allclose(computed, profiles.outside_workspace__phi, rtol=1e-11, atol=1e-12)
    source_receipt = json.loads(track(SOURCE / "provenance.json").read_text())
    # Denoised-BOLD cached means were calculated in float32. Preserve their
    # previously recorded round-off bound; deconvolved branches are float64.
    for representation, block in signals.groupby("branch"):
        tolerance = source_receipt["max_temporal_identity_error"] + 1e-12 if representation == "no_deconvolution" else 1e-12
        for scope in ("whole_brain", "within_workspace", "outside_workspace", "between"):
            np.testing.assert_allclose(block[scope + "__phi"],
                block[scope + "__tdmi"] - block[scope + "__self_sum"], rtol=0, atol=tolerance)
    order = {}
    for line in track(ROOT / "data/aachen_itbs/derivatives/README.txt").read_text().splitlines():
        match = re.match(r"^\s*(\d+)\s+([1342/]+)\s+", line)
        if match:
            order["aachen_itbs_" + match[1].zfill(2)] = 120 if match[2].startswith("3") else 90
    keys = {"estimator_sensitivity_results.csv": ["metric", "family"], "temporal_component_results.csv": ["branch", "metric", "family"],
            "network_results.csv": ["branch", "metric", "family"], "visit_consistency_results.csv": ["branch", "metric", "family"],
            "baseline_motion_sensitivity.csv": ["branch", "metric", "model"], "participant_influence.csv": ["branch", "metric", "omitted_subject"],
            "nonparametric_checks.csv": ["branch", "metric"], "method_interactions.csv": ["metric", "family"]}
    with threadpool_limits(limits=2):
        tables, sp, pp = run(signals, profiles, order, OUT / "full_sample_reproduction")
        parity = [compare(table, pd.read_csv(track(SOURCE / filename)), keys[filename], filename)
                  for filename, table in tables.items() if filename in keys]
        for filename, table in (("signal_participant_changes.csv", sp), ("network_participant_changes.csv", pp)):
            parity.append(compare(table, pd.read_csv(track(SOURCE / filename)), ["branch", "subject_id", "intensity"], filename))
        save_json(OUT / "full_sample_reproduction/parity.json", {"passed": True, "checks": parity})
        receipt = json.loads(track(AUDIT / "audit_receipt.json").read_text())
        if receipt["state"] != "complete" or not receipt["all_193_windows_verified"]:
            raise ValueError("Exact-window QC audit is not complete")
        for filename in ("scan_audit.csv", "cohort_subjects.json"):
            track(AUDIT / filename, receipt["output_sha256"][filename])
        audit = pd.read_csv(AUDIT / "scan_audit.csv").query("dataset == 'aachen_itbs'")
        cohorts = json.loads((AUDIT / "cohort_subjects.json").read_text())
        summary = []
        for branch in ("registered", "strict"):
            subjects = cohorts[branch]["aachen_itbs"]
            selected = audit.loc[audit.subject_id.isin(subjects)]
            if len(selected) != len(subjects) * 4 or not selected[f"{branch}_pass"].all():
                raise ValueError("Aachen complete-cohort whitelist does not match scan QC")
            fd = selected[["subject_id", "session", "retained_mean_fd"]]
            frames = []
            for frame in (signals, profiles):
                subset = frame.loc[frame.subject_id.isin(subjects)].merge(fd, on=["subject_id", "session"], validate="many_to_one")
                subset["mean_fd"] = subset.retained_mean_fd
                frames.append(subset.drop(columns="retained_mean_fd"))
            tables, _, _ = run(*frames, order, OUT / branch, independently_check=True)
            for filename, table in tables.items():
                if "p_maxT_family" in table:
                    group = [column for column in ("branch", "family") if column in table]
                    for name, block in table.groupby(group):
                        summary.append({"qc_branch": branch, "participants": len(subjects), "table": filename,
                                        "analysis_family": str(name), "correction": "exact maxT", "tests": len(block),
                                        "minimum_adjusted_p": float(block.p_maxT_family.min()),
                                        "survivors": int(block.p_maxT_family.lt(.05).sum())})
                elif "p_holm6" in table:
                    for representation, block in table.groupby("branch"):
                        summary.append({"qc_branch": branch, "participants": len(subjects), "table": filename,
                                        "analysis_family": representation + "/baseline_motion6", "correction": "Holm6",
                                        "tests": len(block), "minimum_adjusted_p": float(block.p_holm6.min()),
                                        "survivors": int(block.p_holm6.lt(.05).sum())})
            save_json(OUT / branch / "cohort.json", {"qc_branch": branch, "participants": len(subjects), "subjects": subjects,
                                                      "four_scans_required": True, "normalization": "Cached corrected reference held fixed"})
        pd.DataFrame(summary).to_csv(OUT / "family_summary.csv", index=False)
        if any(sha(ROOT / path) != value for path, value in HASHES.items()):
            raise ValueError("A source changed during rerun")
        save_json(OUT / "completion_receipt.json", {
            "created_utc": datetime.now(timezone.utc).isoformat(), "completed": True,
            "full_sample_original_result_parity_passed": True, "exact_signflip_independent_checks_passed": True,
            "normalization": "Cached original corrected-baseline Phi quantiles held fixed; normalization references are not rebuilt by this rerun",
            "data": "Corrected Aachen only; cached regional/network scalar scores; complete four-scan participants",
            "motion_covariates": "Exact retained-window mean FD from verified scan audit",
            "families_unchanged": ["estimator12", "method8", "temporal9 per representation", "network blocks per representation",
                                    "network omissions8 per representation", "visit12 per representation", "baseline-motion6 per representation"],
            "not_run": "HRF parameter changes and HRF associations are excluded from this task",
            "interpretation": "Method/robustness diagnostics in one unshammed cohort; not independent intervention evidence or cross-dataset replication",
            "source_sha256": HASHES, "script_sha256": sha(Path(__file__)),
            "output_sha256": {path.relative_to(OUT).as_posix(): sha(path) for path in OUT.rglob("*") if path.is_file()}})
    print(pd.DataFrame(summary).to_string(index=False), flush=True)
    print("COMPLETE", OUT, flush=True)


if __name__ == "__main__":
    main()
