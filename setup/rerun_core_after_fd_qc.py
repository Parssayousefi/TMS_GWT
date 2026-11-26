"""Re-estimate existing TMS statistics after enforcing scan-level FD QC.

Read-only inputs: corrected Aachen scalars, SUDMEX scalars and frozen HCP mask.
The cached normalized scores use the original calibration; a separate option
accepts freshly recalibrated scores. No p-value or endpoint selection is used.
Group models follow the user's selected QC-only motion strategy: no continuous
motion covariates are fitted. Motion QC and denoising remain unchanged.
"""
from pathlib import Path
from datetime import datetime, timezone
from itertools import product
import argparse
import hashlib
import json
import re
import sys

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "setup"), str(ROOT / "pipeline/scripts/analysis")]
import analyse_tms_workspace_scope as ws
import test_corrected_gwt_expectations as gwt
from analyse_tms_all_post_average import POST

SOURCE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
OUT = ROOT / "output/tms_fd_reanalysis_20261126/core_qc_only"
KEY = ws.KEY
PLAN_COLUMNS = [*KEY[:2], "arm", "comparison", "baseline_session", "followup_session"]
FAMILIES = {m: m.split("__")[0] if "__" in m else ws.tms.METRIC_INFO[m][0] for m in ws.METRICS}
GROUP_MODEL_DESIGNS = set()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_verified_audit(path):
    """Fail closed unless the final receipt and its inputs still match."""
    path = Path(path).resolve()
    receipt_path = path.parent / "audit_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if (receipt.get("state") != "complete" or
            receipt.get("all_193_windows_verified") is not True or
            receipt.get("changed_source_files") != []):
        raise ValueError("FD audit has not released all 193 exact input windows")
    output_hashes = receipt.get("output_sha256", {})
    required_outputs = {"scan_audit.csv", "pair_audit.csv", "cohort_subjects.json",
                        "baseline_pool_audit.csv", "baseline_failures.csv"}
    if path.name != "scan_audit.csv" or not required_outputs.issubset(output_hashes):
        raise ValueError("FD audit receipt is incomplete or refers to another whitelist")
    for relative, expected in output_hashes.items():
        source = (path.parent / relative).resolve()
        if not source.is_relative_to(path.parent) or sha(source) != expected:
            raise ValueError(f"FD audit output changed after release: {relative}")
    input_hashes = receipt.get("source_sha256", {})
    expected_raw = str((SOURCE / "raw_scan_metrics.csv").relative_to(ROOT)).replace("\\", "/")
    if expected_raw not in input_hashes:
        raise ValueError("FD audit does not identify the corrected source inventory")
    for relative, expected in input_hashes.items():
        source = (ROOT / relative).resolve()
        if not source.is_relative_to(ROOT) or sha(source) != expected:
            raise ValueError(f"FD audit source changed after release: {relative}")
    audit = pd.read_csv(path)
    required = {*KEY, "prefix", "bounds_verified", "legacy_displacement_pass",
                "registered_pass", "strict_pass", "retained_mean_fd", "n_retained",
                "confounds_start", "confounds_end"}
    if not required.issubset(audit) or audit.duplicated(KEY).any() or audit.prefix.duplicated().any():
        raise ValueError("Incomplete or duplicate FD audit keys")
    if audit.groupby("dataset").size().to_dict() != {"aachen_itbs": 60, "sudmex_tms": 133}:
        raise ValueError("FD audit does not cover the original 193-scan inventory")
    for col in ("bounds_verified", "legacy_displacement_pass", "registered_pass", "strict_pass"):
        if audit[col].dtype != bool:
            raise ValueError(f"Invalid or missing audit boolean: {col}")
    if not audit.bounds_verified.all() or (audit.strict_pass & ~audit.registered_pass).any():
        raise ValueError("Unverified bounds or inconsistent nested inclusion flags")
    if ((audit.registered_pass | audit.strict_pass) & ~audit.legacy_displacement_pass).any():
        raise ValueError("Whitelisted scan failed its existing displacement gate")
    numeric = audit[["retained_mean_fd", "n_retained", "confounds_start", "confounds_end"]].to_numpy(float)
    if (not np.isfinite(numeric).all() or (audit.retained_mean_fd < 0).any() or
            not (audit.confounds_end - audit.confounds_start).eq(audit.n_retained).all()):
        raise ValueError("Invalid retained-window motion or boundaries")
    inventory = pd.read_csv(SOURCE / "raw_scan_metrics.csv")
    if set(map(tuple, audit[[*KEY, "prefix"]].to_numpy())) != set(map(tuple, inventory[[*KEY, "prefix"]].to_numpy())):
        raise ValueError("Audit scan identities differ from the corrected source")
    return audit, receipt_path


def verify_corrected_cache():
    """Check saved scalars against their earlier independent release receipt."""
    review_path = ROOT / "PROVENANCE/corrected_release_20261126/corrected_core_review.json"
    review = json.loads(review_path.read_text(encoding="utf-8"))
    if review.get("status") != "passed_bounded_review_no_material_discrepancies":
        raise ValueError("Corrected source cache lacks its independent release review")
    expected = {str(Path(k).resolve() if Path(k).is_absolute() else (ROOT / k).resolve()): v
                for k, v in review["input_and_artifact_sha256"].items()}
    required = [SOURCE / name for name in ("raw_scan_metrics.csv", "normalized_scan_metrics.csv",
                 "raw_paired_changes.csv", "normalized_paired_changes.csv", "raw_input_sha256.json",
                 "statistics_source_sha256.json", "normalization/verification.json",
                 "normalization/extended_sudmex_reuse_guard.json")]
    for path in required:
        if str(path.resolve()) not in expected or sha(path) != expected[str(path.resolve())]:
            raise ValueError(f"Corrected source artifact changed since independent review: {path.name}")
    raw_receipt = json.loads((SOURCE / "raw_input_sha256.json").read_text(encoding="utf-8"))
    references = [(ROOT / k, value) for k, value in raw_receipt.items()
                  if "group_analysis" in k or "atlases" in k]
    if not any(path.name == "hcp_rest_workspace_mask.npy" for path, _ in references):
        raise ValueError("Corrected source receipt lacks the frozen HCP workspace")
    for path, expected_hash in references:
        if sha(path) != expected_hash:
            raise ValueError(f"Fixed reference changed since corrected extraction: {path.name}")
    statistical_sources = json.loads((SOURCE / "statistics_source_sha256.json").read_text(encoding="utf-8"))
    for relative, expected_hash in statistical_sources.items():
        if sha(ROOT / relative) != expected_hash:
            raise ValueError(f"Corrected statistical input changed: {relative}")
    return [review_path, *required, *[path for path, _ in references],
            *[ROOT / relative for relative in statistical_sources]]


def verify_recalibrated_input(path, policy, audit, audit_path):
    """Accept a fresh reference only after both policies finish verification."""
    path = Path(path).resolve()
    if path.name != "normalized_scan_metrics.csv" or path.parent.name != policy:
        raise ValueError("Recalibrated input must identify its policy and combined core scan table")
    folder = path.parent.parent
    global_path, policy_path = folder / "verification.json", path.parent / "verification.json"
    manifest_path, sources_path = folder / "output_sha256.json", folder / "source_sha256.json"
    global_receipt = json.loads(global_path.read_text(encoding="utf-8"))
    policy_receipt = json.loads(policy_path.read_text(encoding="utf-8"))
    if (global_receipt.get("passed") is not True or global_receipt.get("state") != "complete" or
            global_receipt.get("sources_unchanged") is not True or policy_receipt.get("passed") is not True or
            policy_receipt.get("state") != "complete" or policy_receipt.get("qc_policy") != policy):
        raise ValueError("Fresh normalization has not completed its release checks")
    output_hashes = json.loads(manifest_path.read_text(encoding="utf-8"))
    key = str(path.relative_to(folder))
    expected = {str(Path(k)): v for k, v in output_hashes.items()}
    if key not in expected or sha(path) != expected[key]:
        raise ValueError("Fresh normalized scalar table differs from its release hash")
    sources = json.loads(sources_path.read_text(encoding="utf-8"))
    audit_key = str(Path(audit_path).resolve().relative_to(ROOT))
    source_hashes = {str(Path(k)): v for k, v in sources.items()}
    if audit_key not in source_hashes or sha(audit_path) != source_hashes[audit_key]:
        raise ValueError("Fresh normalization used a different scan-QC whitelist")
    scans = pd.read_csv(path)
    allowed = audit[audit[policy + "_pass"]]
    if (scans.duplicated(KEY).any() or
            set(map(tuple, scans[KEY].to_numpy())) != set(map(tuple, allowed[KEY].to_numpy())) or
            len(scans) != policy_receipt.get("passing_scans")):
        raise ValueError("Fresh normalization does not cover exactly the QC-passing scans")
    if not np.isfinite(scans[ws.METRICS].to_numpy(float)).all():
        raise ValueError("Fresh normalized core endpoints are incomplete")
    pool_expected = (audit[audit.normalization_baseline_pool_member & audit[policy + "_pass"]]
                     .groupby("dataset").size().mul(1024).to_dict())
    if pool_expected != policy_receipt.get("baseline_pool_sizes"):
        raise ValueError("Fresh baseline prior does not match the declared QC-clean baseline pool")
    return [global_path, policy_path, manifest_path, sources_path]


def conjunction_p(aachen_p, sudmex_p):
    """An unestimable cohort cannot supply evidence for a conjunction."""
    return float(np.maximum(aachen_p, sudmex_p))


def group_hc3(frame, outcome, predictors, term):
    """Enforce the selected group-level policy at every regression call."""
    forbidden = [name for name in predictors if re.search(r"(^|_)fd($|_)|motion|framewise|^trans_|^rot_", name.lower())]
    if forbidden:
        raise ValueError(f"Continuous motion covariates are excluded from group models: {forbidden}")
    GROUP_MODEL_DESIGNS.add((outcome, tuple(predictors), term))
    return {**ws.hc3(frame, outcome, predictors, term), "predictors": "|".join(predictors),
            "continuous_motion_covariate": False}


def self_check_no_group_motion():
    """Independent HC3 and sign-flip checks using no motion predictor."""
    rng = np.random.default_rng(2026112604)
    frame = pd.DataFrame({"y": rng.normal(size=40), "intercept": 1., "active": np.arange(40) % 2,
                          "baseline": rng.normal(size=40)})
    predictors = ["intercept", "active", "baseline"]
    fitted = group_hc3(frame, "y", predictors, "active")
    independent = gwt.direct_hc3(frame, "y", predictors, "active")
    for name in ("estimate", "se", "p", "ci_low", "ci_high"):
        np.testing.assert_allclose(fitted[name], independent[name], atol=1e-12)
    values = rng.normal(size=(6, 4))
    p, adjusted = gwt.exact_signflip(values)
    null = np.asarray([np.abs(ws.stats.ttest_1samp(values * np.asarray(signs)[:, None], 0, axis=0).statistic)
                       for signs in product((-1., 1.), repeat=len(values))])
    observed = abs(ws.stats.ttest_1samp(values, 0, axis=0).statistic)
    np.testing.assert_array_equal(p, (null >= observed - 1e-10).mean(0))
    np.testing.assert_array_equal(adjusted, (null.max(1)[:, None] >= observed - 1e-10).mean(0))


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")


def retain_plan(scans, plan):
    """Require complete pairs and both Aachen visits, never available-visit means."""
    if scans.duplicated(KEY).any() or plan.duplicated(["dataset", "subject_id", "comparison"]).any():
        raise ValueError("Duplicate scan or contrast keys")
    keys = set(map(tuple, scans[KEY].to_numpy()))
    keep = [(p.dataset, p.subject_id, p.baseline_session) in keys and
            (p.dataset, p.subject_id, p.followup_session) in keys for p in plan.itertuples()]
    result = plan.loc[keep].copy()
    a = result[result.dataset.eq("aachen_itbs")]
    eligible = a.groupby("subject_id").comparison.agg(set)
    eligible = eligible[eligible.map(lambda s: s == {"post90_minus_pre90", "post120_minus_pre120"})].index
    result = result[result.dataset.ne("aachen_itbs") | result.subject_id.isin(eligible)]
    sud = result[result.dataset.eq("sudmex_tms")]
    assert (sud.baseline_session.eq("t0") & sud.followup_session.eq("t1")).all()
    assert not sud.subject_id.duplicated().any()
    return result.reset_index(drop=True)


def arm_tests_qc_only(pairs, metrics, families):
    frame = pairs.assign(intercept=1., active=pairs.arm.eq("active").astype(float))
    rows = []
    for m in metrics:
        rows.append({"model": "unadjusted_change_difference", "metric": m, "family": families[m],
                     **ws.tms.welch_summary(frame.loc[frame.active.eq(1), m + "_delta"],
                                            frame.loc[frame.active.eq(0), m + "_delta"])})
        rows.append({"model": "baseline_adjusted", "metric": m, "family": families[m],
                     **group_hc3(frame, m + "_followup", ["intercept", "active", m + "_baseline"], "active")})
    result = pd.DataFrame(rows)
    result["p_holm_family"] = result.groupby(["model", "family"]).p.transform(ws.holm)
    return result


def sudmex_tests(pairs):
    sud = pairs[pairs.dataset.eq("sudmex_tms")]
    result = arm_tests_qc_only(sud, ws.METRICS, FAMILIES)
    result["p_holm_all40"] = result.groupby("model").p.transform(ws.holm)
    for label, metrics in [("main30", ws.MAIN_METRICS), ("information20", ws.INFO_METRICS),
                           ("whole_and_workspace10", ws.MAIN_METRICS[:10]), ("workspace25", ws.WS25)]:
        mask = result.metric.isin(metrics)
        result.loc[mask, "p_holm_" + label] = result.loc[mask].groupby("model").p.transform(ws.holm)
    result["n_active"] = int(sud.arm.eq("active").sum())
    result["n_sham"] = int(sud.arm.eq("sham").sum())
    return result


def gwt_tests(pairs, branch):
    pooled = ws.aggregate_aachen(pairs)
    pflip, pmax = gwt.exact_signflip(pooled[[m + "_delta" for m in gwt.METRICS]].to_numpy())
    rows = []
    for i, m in enumerate(gwt.METRICS):
        rows.append({"dataset": "aachen_itbs", "metric": m, "model": "mean of two visit changes",
                     **ws.tms.paired_summary(pooled[m + "_delta"]), "p_sign_flip": pflip[i],
                     "p_maxT4": pmax[i], "sign_vectors": 2 ** len(pooled)})
    sud = pairs[pairs.dataset.eq("sudmex_tms")].copy()
    sud["intercept"], sud["active"] = 1., sud.arm.eq("active").astype(float)
    for m in gwt.METRICS:
        rows.append({"dataset": "sudmex_tms", "metric": m,
                     "model": "followup ~ active + baseline outcome; HC3/t",
                     "n_active": int(sud.active.sum()), "n_sham": int((1 - sud.active).sum()),
                     **group_hc3(sud, m + "_followup", ["intercept", "active", m + "_baseline"], "active")})
    result = pd.DataFrame(rows)
    result["p_holm4"] = result.groupby("dataset").p.transform(ws.holm)
    result["branch"], result["expected_direction"], result["test_sidedness"] = branch, "increase", "two-sided"
    result["family_size"] = 4
    return result


def all_post(scans):
    """Rebuild each person's average from passing visits, with passing t0 required."""
    rows = []
    for (dataset, subject), group in scans[scans.dataset.eq("sudmex_tms")].groupby(["dataset", "subject_id"]):
        if "t0" not in set(group.session):
            continue
        pre = group[group.session.eq("t0")].iloc[0]
        post = group[group.session.isin(POST)]
        if post.empty:
            continue
        assert post.valid_roi_signature.eq(pre.valid_roi_signature).all()
        assert post.arm.eq(pre.arm).all()
        row = {"dataset": dataset, "subject_id": subject, "arm": pre.arm, "n_pre": 1, "n_post": len(post),
               "pre_sessions": "t0", "post_sessions": "|".join(sorted(post.session, key=POST.index)),
               "delta_fd": float(post.mean_fd.mean() - pre.mean_fd)}
        for m in ws.METRICS:
            row.update({m + "_baseline": pre[m], m + "_followup": post[m].mean(), m + "_delta": post[m].mean() - pre[m]})
        rows.append(row)
    return pd.DataFrame(rows)


def spatial_sudmex(pairs, branch):
    sud = pairs[pairs.dataset.eq("sudmex_tms")].copy()
    sud["intercept"], sud["active"] = 1., sud.arm.eq("active").astype(float)
    rows = []
    for m in ws.MEASURES:
        for time in ["baseline", "followup"]:
            sud["spatial_" + time] = sud["outside_workspace__" + m + "_" + time] - sud["within_workspace__" + m + "_" + time]
        rows.append({"branch": branch, "metric": m, "model": "baseline_adjusted", "contrast": "outside_minus_within_workspace",
                     **group_hc3(sud, "spatial_followup", ["intercept", "active", "spatial_baseline"], "active")})
    result = pd.DataFrame(rows)
    result["p_holm5"] = result.groupby("model").p.transform(ws.holm)
    return result


def information_leave_one_out(pooled, branch):
    """Influence diagnostic for every information endpoint, not a new search."""
    rows = []
    for subject in pooled.index:
        result = ws.test_aachen(pooled.drop(index=subject), ws.INFO_METRICS)
        result["branch"] = branch
        result["omitted_subject"] = subject
        result["diagnostic_only"] = True
        result["correction_scope"] = "scope5 and all20 within each omission; omissions are influence diagnostics"
        rows.append(result)
    return pd.concat(rows, ignore_index=True)


def nonmotion_robustness(out, pairs_by_branch):
    robust, loo, intensity, spatial = [], [], [], []
    order = {}
    for line in (ROOT / "data/aachen_itbs/derivatives/README.txt").read_text().splitlines():
        found = re.match(r"^\s*(\d+)\s+([1342/]+)\s+", line)
        if found:
            order["aachen_itbs_" + found[1].zfill(2)] = int(found[2].startswith("3"))
    for branch, pairs in pairs_by_branch.items():
        pooled = ws.aggregate_aachen(pairs)
        a90 = pairs[pairs.comparison.eq("post90_minus_pre90")].set_index("subject_id").sort_index()
        a120 = pairs[pairs.comparison.eq("post120_minus_pre120")].set_index("subject_id").sort_index()
        for m in ws.MAIN_METRICS[:10]:
            family, values = m.split("__")[0], pooled[m + "_delta"].to_numpy()
            d90, d120 = a90[m + "_delta"], a120[m + "_delta"]
            for contrast, change in [("90_change", d90), ("120_change", d120), ("120_minus_90", d120-d90)]:
                intensity.append({"branch": branch, "metric": m, "family": family, "contrast": contrast,
                                  **ws.tms.paired_summary(change)})
            frame = pd.DataFrame({"y": d120-d90, "intercept": 1.,
                                  "order": np.array([order[s] for s in a90.index])-.5,
                                  "baseline": a120[m + "_baseline"]-a90[m + "_baseline"]})
            intensity.append({"branch": branch, "metric": m, "family": family,
                              "contrast": "120_minus_90_order_baseline",
                              **group_hc3(frame, "y", ["intercept", "order", "baseline"], "intercept")})
            if m.endswith(("__phi", "__phi_r")):
                rng = np.random.default_rng(2026112526)
                means = values[rng.integers(len(values), size=(10000, len(values)))].mean(1)
                low, high = np.quantile(means, [.025, .975])
                robust.extend([{"branch": branch, "metric": m, "family": family, "check": "bootstrap10000",
                                "estimate": values.mean(), "ci_low": low, "ci_high": high},
                               {"branch": branch, "metric": m, "family": family, "check": "wilcoxon",
                                "estimate": np.median(values), "p": ws.stats.wilcoxon(values, method="exact").pvalue}])
                for k, subject in enumerate(pooled.index):
                    loo.append({"branch": branch, "metric": m, "omitted_subject": subject,
                                **ws.tms.paired_summary(np.delete(values, k))})
        differences = pd.DataFrame({m + "_delta": pooled["outside_workspace__" + m + "_delta"]-
                                     pooled["within_workspace__" + m + "_delta"] for m in ws.MEASURES})
        p, adjusted = ws.exact_max_t(differences)
        for k, m in enumerate(ws.MEASURES):
            spatial.append({"branch": branch, "metric": m, "contrast": "outside_minus_within_workspace",
                            **ws.tms.paired_summary(differences[m + "_delta"]), "p_sign_flip": p[k], "p_maxT_5": adjusted[k]})
    pd.DataFrame(robust).to_csv(out / "aachen_robustness.csv", index=False)
    pd.DataFrame(loo).to_csv(out / "aachen_leave_one_out.csv", index=False)
    intensity = pd.DataFrame(intensity)
    intensity["contrast_family"] = np.where(intensity.contrast.eq("120_minus_90_order_baseline"),
                                            "adjusted_intensity", "three_intensity_contrasts")
    intensity["p_holm_family"] = intensity.groupby(["branch", "family", "contrast_family"]).p.transform(ws.holm)
    intensity.to_csv(out / "aachen_intensity_results.csv", index=False)
    spatial = pd.DataFrame(spatial)
    spatial["p_holm5"] = spatial.groupby("branch").p.transform(ws.holm)
    spatial.to_csv(out / "aachen_spatial_interactions.csv", index=False)


def clinical_qc_only(out, pairs_by_branch):
    rows = []
    neural = ["gateways_FPN", "gateways_DAN", "mean_synergy_pc", "mean_redundancy_pc",
              "within_workspace__phi_r", "within_workspace__phi", "whole_brain__phi_r", "whole_brain__phi"]
    for branch, pairs in pairs_by_branch.items():
        sud = pairs[pairs.dataset.eq("sudmex_tms")].drop(columns=["baseline_session", "followup_session"])
        frame = ws.merge_clinical_neural(sud).assign(active=lambda d: d.arm.eq("active").astype(float), intercept=1.)
        for cm in ws.CLINICAL:
            for nm in neural:
                data = frame.copy()
                data["neural_change_z"] = (data[nm + "_delta"]-data[nm + "_delta"].mean())/data[nm + "_delta"].std(ddof=1)
                result = group_hc3(data, cm + "_delta", ["intercept", "neural_change_z", "active", cm + "_t0", nm + "_baseline"], "neural_change_z")
                rows.append({"branch": branch, "clinical_metric": cm, "neural_metric": nm, **result})
    result = pd.DataFrame(rows)
    result["p_holm56"] = result.groupby("branch").p.transform(ws.holm)
    selected = result.neural_metric.isin(neural[:5])
    result.loc[selected, "p_holm_workspace35"] = result[selected].groupby("branch").p.transform(ws.holm)
    result.to_csv(out / "clinical_neural_associations.csv", index=False)


def compare_numeric(actual, expected, keys, columns):
    a, b = actual.set_index(keys).sort_index(), expected.set_index(keys).sort_index()
    assert a.index.equals(b.index), (len(a), len(b))
    for c in columns:
        np.testing.assert_allclose(a[c].to_numpy(float), b[c].to_numpy(float), atol=2e-10, rtol=1e-8, equal_nan=True, err_msg=c)


def parity(out, plan):
    records = []
    for branch in ["raw", "normalized"]:
        scans = pd.read_csv(SOURCE / (branch + "_scan_metrics.csv"))
        pairs = ws.make_pairs(scans, plan)
        previous = pd.read_csv(SOURCE / (branch + "_paired_changes.csv"))
        compare_numeric(pairs, previous, ["dataset", "subject_id", "comparison"],
                        ["delta_fd", *[m + suffix for m in ws.METRICS for suffix in ["_baseline", "_followup", "_delta"]]])
        a = ws.test_aachen(ws.aggregate_aachen(pairs), ws.METRICS)
        old = pd.read_csv(SOURCE / (branch + "_aachen_all40.csv"))
        compare_numeric(a, old, ["metric"], ["estimate", "p", "p_sign_flip", "p_maxT_family", "p_holm_family", "p_holm_all", "p_maxT_all"])
        s = sudmex_tests(pairs)
        old = pd.read_csv(SOURCE / (branch + "_sudmex_acute.csv"))
        compare_numeric(s[s.metric.isin(ws.MAIN_METRICS)], old[old.model.isin(s.model.unique())],
                        ["model", "metric"], ["estimate", "p", "p_holm_family"])
        g = gwt_tests(pairs, branch)
        old = pd.read_csv(SOURCE / ("gwt_expectations_" + branch + ".csv"))
        compare_numeric(g[g.dataset.eq("aachen_itbs")], old[old.dataset.eq("aachen_itbs")],
                        ["dataset", "metric"], ["estimate", "p", "p_holm4", "p_maxT4"])
        compare_numeric(g[g.dataset.eq("sudmex_tms")],
                        s[s.model.eq("baseline_adjusted") & s.metric.isin(gwt.METRICS)],
                        ["metric"], ["estimate", "p"])
        ap = all_post(scans)
        old = pd.read_csv(SOURCE / (branch + "_all_post_all40_changes.csv"))
        compare_numeric(ap, old, ["dataset", "subject_id"], ["n_post", "delta_fd", *[m + "_delta" for m in ws.METRICS]])
        records.append({"branch": branch, "paired_values": "matched", "all40_Aachen": "matched", "SUDMEX_main30": "matched",
                        "GWT4_Aachen": "matched", "GWT4_SUDMEX": "matches prior baseline-only models; group FD removed by user choice",
                        "all_post_values": "matched", "group_motion_models_refit": False, "atol": 2e-10, "rtol": 1e-8})
    save_json(out / "full_sample_reproduction.json", records)


def analyze(out, scans_by_branch, plan):
    out.mkdir(parents=True, exist_ok=True)
    pairs_by_branch = {}
    allpost_results, allpost_arms, spatial, counts, conjunction, information_loo = [], [], [], [], [], []
    for branch, scans in scans_by_branch.items():
        selected_plan = retain_plan(scans, plan)
        pairs = ws.make_pairs(scans, selected_plan)
        pooled = ws.aggregate_aachen(pairs)
        assert len(pooled) >= 6
        pairs_by_branch[branch] = pairs
        scans.to_csv(out / (branch + "_scan_metrics.csv"), index=False)
        pairs.to_csv(out / (branch + "_paired_changes.csv"), index=False)
        pooled.to_csv(out / (branch + "_aachen_participant_changes.csv"))
        aa_tables = {}
        for label, metrics in [("main30", ws.MAIN_METRICS), ("information20", ws.INFO_METRICS), ("all40", ws.METRICS)]:
            table = ws.test_aachen(pooled, metrics)
            table["branch"] = branch
            table.to_csv(out / (branch + "_aachen_" + label + ".csv"), index=False)
            aa_tables[label] = table
        sud = sudmex_tests(pairs)
        sud["branch"] = branch
        sud.to_csv(out / (branch + "_sudmex_acute_all40.csv"), index=False)
        gwt_tests(pairs, branch).to_csv(out / ("gwt_expectations_" + branch + ".csv"), index=False)
        information_loo.append(information_leave_one_out(pooled, branch))
        spatial.append(spatial_sudmex(pairs, branch))
        ap = all_post(scans)
        ap.to_csv(out / (branch + "_all_post_all40_changes.csv"), index=False)
        test = ws.paired_family_tests(ap, ws.METRICS, FAMILIES, 2026112527)
        test["branch"] = branch
        allpost_results.append(test)
        arm = arm_tests_qc_only(ap, ws.METRICS, FAMILIES)
        arm["branch"] = branch
        allpost_arms.append(arm)
        counts.append({"branch": branch, "aachen_participants": len(pooled), "aachen_visit_pairs": len(pairs[pairs.dataset.eq("aachen_itbs")]),
                       "sudmex_acute_participants": int(pairs.dataset.eq("sudmex_tms").sum()),
                       "sudmex_acute_arms": pairs[pairs.dataset.eq("sudmex_tms")].arm.value_counts().to_dict(),
                       "sudmex_all_post_participants": len(ap), "all_post_n_visits": ap.n_post.value_counts().to_dict()})
        # Explicit conjunction requires evidence in both cohorts; same direction
        # is an additional condition, never a reason to select a one-sided test.
        aa = aa_tables["all40"].set_index("metric")
        ss = sud[sud.model.eq("baseline_adjusted")].set_index("metric")
        for m in ws.METRICS:
            conjunction.append({"branch": branch, "metric": m, "family": FAMILIES[m],
                                "aachen_estimate": aa.loc[m, "estimate"], "sudmex_estimate": ss.loc[m, "estimate"],
                                "same_direction": bool(aa.loc[m, "estimate"] * ss.loc[m, "estimate"] > 0),
                                "p_conjunction": conjunction_p(aa.loc[m, "p_sign_flip"], ss.loc[m, "p"])})
        print(branch, counts[-1], flush=True)
    nonmotion_robustness(out, pairs_by_branch)
    clinical_qc_only(out, pairs_by_branch)
    pd.concat(allpost_results).to_csv(out / "sudmex_all_post_all40_results.csv", index=False)
    pd.concat(allpost_arms).to_csv(out / "sudmex_all_post_arm_results.csv", index=False)
    pd.concat(spatial).to_csv(out / "sudmex_spatial_interactions.csv", index=False)
    conjunction = pd.DataFrame(conjunction)
    conjunction["p_holm_family"] = conjunction.groupby(["branch", "family"]).p_conjunction.transform(ws.holm)
    conjunction["p_holm_all40"] = conjunction.groupby("branch").p_conjunction.transform(ws.holm)
    conjunction.to_csv(out / "cross_dataset_conjunction40.csv", index=False)
    pd.concat(information_loo).to_csv(out / "aachen_information_leave_one_out.csv", index=False)
    save_json(out / "cohort_counts.json", counts)


def check_pair_filter():
    plan = pd.DataFrame([
        {"dataset": "aachen_itbs", "subject_id": "A", "arm": "active", "comparison": f"post{x}_minus_pre{x}",
         "baseline_session": f"pre{x}", "followup_session": f"post{x}"} for x in [90, 120]
    ] + [{"dataset": "sudmex_tms", "subject_id": "S", "arm": "sham", "comparison": "t1_minus_t0", "baseline_session": "t0", "followup_session": "t1"}])
    scans = pd.DataFrame([{"dataset": p.dataset, "subject_id": p.subject_id, "session": s}
                          for p in plan.itertuples() for s in [p.baseline_session, p.followup_session]])
    assert len(retain_plan(scans, plan)) == 3
    assert retain_plan(scans[scans.session.ne("post90")], plan).dataset.tolist() == ["sudmex_tms"]
    assert len(retain_plan(scans[scans.session.ne("t0")], plan)) == 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--registered-column", default="registered_pass")
    parser.add_argument("--strict-column", default="strict_pass")
    parser.add_argument("--recalibrated-registered", type=Path)
    parser.add_argument("--recalibrated-strict", type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    assert out.is_relative_to(ROOT)
    archive = ROOT / "output/tms_fd_reanalysis_20261126/core"
    if out == archive or out.is_relative_to(archive):
        raise ValueError("Earlier FD-covariate results are archived; use the separate QC-only output")
    out.mkdir(parents=True, exist_ok=True)
    self_check_no_group_motion()
    check_pair_filter()
    audit, audit_receipt_path = load_verified_audit(args.audit)
    cache_sources = verify_corrected_cache()
    plan = pd.read_csv(SOURCE / "raw_paired_changes.csv")[PLAN_COLUMNS]
    parity(out, plan)
    assert not audit.duplicated(KEY).any()
    sources = [args.audit, audit_receipt_path, Path(__file__), SOURCE / "raw_paired_changes.csv", Path(ws.__file__), Path(gwt.__file__),
               ROOT / "setup/analyse_tms_all_post_average.py", ROOT / "setup/analyse_tms_familywise.py",
               ROOT / "pipeline/scripts/utils/clinical_data.py", *cache_sources]
    specification = {"created_utc": datetime.now(timezone.utc).isoformat(), "primary_qc": "Existing displacement checks plus configured mean FD <=0.5 mm, peak FD <=3 mm and min100 retained volumes",
                     "strict_sensitivity": "Additionally mean FD <=0.30 mm and fraction FD >0.5 <=0.20; thresholds chosen before result review",
                     "Aachen": "Both complete intensity pairs required, averaged equally per participant",
                     "SUDMEX": "Acute t0/t1 randomized comparison: unadjusted change difference and baseline-outcome-adjusted HC3/t; baseline-adjusted model supplies GWT and cross-dataset comparisons",
                     "FD": "Motion enters scan QC and descriptive metadata only. No additional continuous group-level motion covariate is fitted, following explicit user choice. Preprocessing denoising remains unchanged.",
                     "group_motion_covariates": False,
                     "clinical": "Clinical change ~ neural change z + randomized arm + clinical baseline + neural baseline; complete-case HC3/t, no motion covariate",
                     "normalization": "normalized_frozen_reference reuses previous calibration; normalized_recalibrated only if supplied and explicitly verified",
                     "families": "Original five measures per scope; role subfamilies; main30, information20 and all40. Four GWT expectations. Five direct spatial contrasts. Clinical56. Cross-dataset conjunction40 additionally shown.",
                     "information_influence": "All20 information endpoints reestimated after each Aachen omission, with scope5 and all20 correction within omission; influence diagnostic, not another inferential family or participant-selection rule.",
                     "all_post": "Rebuild each person's average from all QC-passing original post sessions (t1,t14,t2,t3,t4); passing t0 required. Open-label descriptive sensitivity, not randomized efficacy.",
                     "scope": "No hemodynamic tests or stability claims; earlier invalid Aachen branch is never used."}
    save_json(out / "analysis_specification.json", specification)
    for policy, column, recalibrated in [("registered", args.registered_column, args.recalibrated_registered),
                                         ("strict", args.strict_column, args.recalibrated_strict)]:
        assert column in audit and "retained_mean_fd" in audit
        flags = audit[column]
        if flags.dtype != bool:
            assert flags.astype(str).str.lower().isin(["true", "false"]).all()
            flags = flags.astype(str).str.lower().eq("true")
        passing = audit.loc[flags, KEY + ["retained_mean_fd"]]
        branches = {}
        inputs = [("raw", SOURCE / "raw_scan_metrics.csv"), ("normalized_frozen_reference", SOURCE / "normalized_scan_metrics.csv")]
        if recalibrated:
            sources.extend(verify_recalibrated_input(recalibrated, policy, audit, args.audit))
            inputs.append(("normalized_recalibrated", recalibrated))
        for branch, path in inputs:
            sources.append(path)
            scans = pd.read_csv(path)
            assert not scans.duplicated(KEY).any()
            assert set(map(tuple, scans[KEY].to_numpy())) <= set(map(tuple, audit[KEY].to_numpy()))
            filtered = scans.merge(passing, on=KEY, how="inner", validate="one_to_one")
            filtered["legacy_mean_fd"] = filtered.mean_fd
            filtered["mean_fd"] = filtered.pop("retained_mean_fd")
            branches[branch] = filtered
        analyze(out / policy, branches, plan)
    save_json(out / "verification.json", {"completed_utc": datetime.now(timezone.utc).isoformat(), "full_sample_parity": "passed",
              "filtering_tests": "passed", "independent_HC3_and_exact_tests": "passed",
              "group_motion_covariates": False, "user_selected_motion_strategy": "motion QC and preprocessing only",
              "group_model_designs": [{"outcome": outcome, "predictors": list(predictors), "term": term}
                                       for outcome, predictors, term in sorted(GROUP_MODEL_DESIGNS)],
              "source_sha256": {str(p.resolve().relative_to(ROOT)): sha(p) for p in set(sources)},
              "outputs_sha256": {str(p.relative_to(out)): sha(p) for p in out.rglob("*.csv")}})
    print("Core QC reanalysis complete", flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
