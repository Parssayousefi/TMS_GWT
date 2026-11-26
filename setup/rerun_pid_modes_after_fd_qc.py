"""Re-estimate existing PID mode tests after the verified retained-frame FD audit.

This runner reuses validated scan-level outcomes and fixed normalized references.
It refits statistics and independent-cohort regional templates after QC. No new
endpoints, spatial masks, production preprocessing, or null references are fitted.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "setup"), str(ROOT / "pipeline/scripts/analysis")]
import explore_pid_processing_modes as original
import explore_pid_change_patterns as regional_original
from analyse_tms_familywise import exact_max_t
from test_corrected_gwt_expectations import direct_hc3

RAW = ROOT / "output/pid_processing_modes/20261126T130022_906563Z"
NORMALIZED = ROOT / "output/pid_modes_normalized/20261126T131608_596458Z"
REGIONAL = ROOT / "output/pid_change_patterns/20261126T130710_947110Z"
CORE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
DEFAULT_AUDIT = ROOT / "output/existing_tms_fd_audit_20261126"
DEFAULT_OUT = ROOT / "output/tms_fd_reanalysis_20261126/modes"
MODES = list(original.MODES)
SCOPES = ["whole_brain", "within_workspace", "outside_workspace", "workspace_boundary",
          "FPN_striatum", "FPN_DAN", "FPN_DMN"]
METRICS = [f"{s}__{m}" for s in SCOPES for m in MODES]
SEEDS = ["2026112511", "2026112512", "combined"]
SOURCES: dict[str, str] = {}
COMPARISONS = []


def track(path: Path) -> Path:
    path = path.resolve()
    SOURCES[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(track(path), **kwargs)


def read_json(path: Path):
    return json.loads(track(path).read_text(encoding="utf-8-sig"))


def save_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + "\n", encoding="utf-8")


def bool_column(series: pd.Series) -> pd.Series:
    normalized = series.astype(str).str.lower()
    if not normalized.isin(["true", "false", "1", "0"]).all():
        raise ValueError(f"Invalid boolean column: {series.name}")
    return normalized.isin(["true", "1"])


def compare_table(actual: pd.DataFrame, source: Path, keys: list[str], label: str):
    """Check all pre-existing columns, including numeric results and row identity."""
    expected = read_csv(source)
    a = actual.set_index(keys).sort_index()
    b = expected.set_index(keys).sort_index()
    if not a.index.equals(b.index):
        raise AssertionError(f"Reproduction row mismatch: {label}")
    maximum = 0.
    for col in b:
        if col not in a:
            raise AssertionError(f"Reproduction missing column {label}: {col}")
        if pd.api.types.is_numeric_dtype(b[col]) and not pd.api.types.is_bool_dtype(b[col]):
            x, y = a[col].to_numpy(float), b[col].to_numpy(float)
            np.testing.assert_allclose(x, y, rtol=1e-8, atol=2e-12, equal_nan=True,
                                       err_msg=f"{label}: {col}")
            good = np.isfinite(x) & np.isfinite(y)
            if good.any():
                maximum = max(maximum, float(abs(x[good] - y[good]).max()))
        else:
            if not a[col].fillna("<NA>").astype(str).equals(b[col].fillna("<NA>").astype(str)):
                raise AssertionError(f"Reproduction value mismatch {label}: {col}")
    COMPARISONS.append(dict(table=label, rows=len(a), columns=len(b.columns),
                            max_absolute_numeric_difference=maximum, passed=True,
                            source=str(source.relative_to(ROOT))))


def make_pairs(scans: pd.DataFrame, plan: pd.DataFrame) -> pd.DataFrame:
    if scans.duplicated(["dataset", "subject_id", "session"]).any():
        raise ValueError("Scan identity is not unique")
    si = scans.set_index(["dataset", "subject_id", "session"])
    records = []
    for row in plan.itertuples():
        before = si.loc[(row.dataset, row.subject_id, row.baseline_session)]
        after = si.loc[(row.dataset, row.subject_id, row.followup_session)]
        if before.arm != row.arm or after.arm != row.arm:
            raise AssertionError("Plan/scan arm mismatch")
        values = {key: getattr(row, key) for key in ["dataset", "subject_id", "arm", "comparison"]}
        values["delta_fd"] = after.mean_fd - before.mean_fd
        for metric in METRICS:
            values[metric + "_pre"] = before[metric]
            values[metric + "_post"] = after[metric]
            values[metric + "_delta"] = after[metric] - before[metric]
        records.append(values)
    pairs = pd.DataFrame(records)
    if not np.isfinite(pairs.select_dtypes(include=[np.number])).all().all():
        raise ValueError("Nonfinite paired values")
    counts = pairs.groupby(["dataset", "subject_id"]).size()
    assert counts.loc["aachen_itbs"].eq(2).all()
    assert counts.loc["sudmex_tms"].eq(1).all()
    return pairs


def mode_inference(scans: pd.DataFrame, plan: pd.DataFrame, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    pairs = make_pairs(scans, plan)
    pairs.to_csv(out / "participant_changes.csv", index=False)
    aa = pairs.loc[pairs.dataset.eq("aachen_itbs")].groupby("subject_id").mean(numeric_only=True).sort_index()
    sud = pairs.loc[pairs.dataset.eq("sudmex_tms")].copy()
    sud["active"] = sud.arm.eq("active").astype(float)
    sud["intercept"] = 1.
    assert 3 <= len(aa) <= 20 and sud.active.nunique() == 2
    values = aa[[m + "_delta" for m in METRICS]].to_numpy()
    pexact, pmax = exact_max_t(values)
    records = []
    for scope in SCOPES:
        columns = [f"{scope}__{m}" for m in MODES]
        pscope = exact_max_t(aa[[m + "_delta" for m in columns]].to_numpy())[1]
        for j, metric in enumerate(columns):
            k = METRICS.index(metric)
            change = aa[metric + "_delta"]
            estimate, se = change.mean(), stats.sem(change)
            ci = estimate + np.array([-1, 1]) * stats.t.ppf(.975, len(aa) - 1) * se
            records.append(dict(dataset="aachen_itbs", scope=scope, mode=metric.split("__")[1],
                                metric=metric, model="paired_change", n=len(aa), estimate=estimate,
                                ci_low=ci[0], ci_high=ci[1], p=pexact[k],
                                p_t=stats.ttest_1samp(change, 0).pvalue, p_maxT6=pscope[j],
                                p_maxT42=pmax[k], pre_mean=aa[metric + "_pre"].mean(),
                                post_mean=aa[metric + "_post"].mean()))
            for model, predictors in [
                ("baseline_motion_adjusted", ["intercept", "active", metric + "_pre", "delta_fd"]),
                ("baseline_adjusted", ["intercept", "active", metric + "_pre"]),
            ]:
                fit = direct_hc3(sud, metric + "_post", predictors, "active")
                if fit["status"] != "ok" or fit["n"] != len(sud):
                    raise AssertionError(f"Unestimable mode model: {metric}/{model}: {fit}")
                check = sm.OLS(sud[metric + "_post"], sud[predictors]).fit(cov_type="HC3", use_t=True)
                np.testing.assert_allclose([fit["estimate"], fit["se"], fit["p"]],
                                           [check.params.active, check.bse.active, check.pvalues.active],
                                           rtol=1e-7, atol=1e-10)
                records.append(dict(dataset="sudmex_tms", scope=scope, mode=metric.split("__")[1],
                                    metric=metric, model=model, n_active=int(sud.active.sum()),
                                    n_sham=int((sud.active == 0).sum()),
                                    active_change=sud.loc[sud.active.eq(1), metric + "_delta"].mean(),
                                    sham_change=sud.loc[sud.active.eq(0), metric + "_delta"].mean(), **fit))
    results = pd.DataFrame(records)
    results["p_holm6"] = results.groupby(["dataset", "model", "scope"]).p.transform(
        lambda p: multipletests(p, method="holm")[1])
    results["p_holm42"] = results.groupby(["dataset", "model"]).p.transform(
        lambda p: multipletests(p, method="holm")[1])
    results.to_csv(out / "mode_treatment_results.csv", index=False)
    primary = results.loc[~results.model.eq("baseline_adjusted")]
    records = []
    for metric in METRICS:
        a = primary.loc[primary.dataset.eq("aachen_itbs") & primary.metric.eq(metric)].iloc[0]
        b = primary.loc[primary.dataset.eq("sudmex_tms") & primary.metric.eq(metric)].iloc[0]
        same = np.sign(a.estimate) == np.sign(b.estimate)
        records.append(dict(metric=metric, scope=a.scope, mode=a["mode"], same_direction=bool(same),
                            aachen_estimate=a.estimate, sudmex_estimate=b.estimate,
                            p_aachen=a.p, p_sudmex=b.p, p=max(a.p, b.p) if same else 1.))
    cc = pd.DataFrame(records)
    cc["p_holm6"] = cc.groupby("scope").p.transform(lambda p: multipletests(p, method="holm")[1])
    cc["p_holm42"] = multipletests(cc.p, method="holm")[1]
    cc.to_csv(out / "cross_dataset_conjunction.csv", index=False)
    return results, cc, pairs


def normalized_inference(scans: pd.DataFrame, plan: pd.DataFrame, out: Path):
    treatments, conjunctions, changes = [], [], []
    for seed in SEEDS:
        destination = out if seed == "combined" else out / seed
        r, c, p = mode_inference(scans.loc[scans.reference_seed.eq(seed)], plan, destination)
        treatments.append(r.assign(reference_seed=seed))
        conjunctions.append(c.assign(reference_seed=seed))
        changes.append(p.assign(reference_seed=seed))
    results, cc, pairs = pd.concat(treatments), pd.concat(conjunctions), pd.concat(changes)
    results.to_csv(out / "mode_treatment_results_all_seeds.csv", index=False)
    cc.to_csv(out / "cross_dataset_conjunction_all_seeds.csv", index=False)
    pairs.to_csv(out / "participant_changes_all_seeds.csv", index=False)
    seed_rows = []
    for dataset in ["aachen_itbs", "sudmex_tms"]:
        for kind, data, keys, cols in [
            ("scans", scans, ["subject_id", "session"], METRICS),
            ("changes", pairs, ["subject_id", "comparison"], [m + "_delta" for m in METRICS]),
        ]:
            selected = data.loc[data.dataset.eq(dataset)]
            a, b = [selected.loc[selected.reference_seed.eq(seed)].set_index(keys).sort_index() for seed in SEEDS[:2]]
            assert a.index.equals(b.index)
            for metric in cols:
                difference = a[metric] - b[metric]
                seed_rows.append(dict(dataset=dataset, kind=kind, metric=metric,
                                      rms_seed_difference=np.sqrt(np.mean(difference ** 2)),
                                      max_absolute_seed_difference=abs(difference).max(),
                                      data_sd=((a[metric] + b[metric]) / 2).std()))
    stability = pd.DataFrame(seed_rows)
    stability.to_csv(out / "monte_carlo_stability.csv", index=False)
    return results, cc, stability


def load_regional(plan: pd.DataFrame, scans: pd.DataFrame):
    """Recover original sorted IDs; verify each saved vector's scalar identity."""
    arrays = np.load(track(REGIONAL / "regional_modes.npz"))
    pp = {}
    for dataset in ["aachen_itbs", "sudmex_tms"]:
        group = plan.loc[plan.dataset.eq(dataset)]
        ids = sorted(group.subject_id.unique())
        value = {field: arrays[f"{dataset}_{field}"].copy() for field in ["pre", "post", "fd", "active"]}
        value["ids"] = ids
        assert value["pre"].shape == value["post"].shape == (len(ids), 6, 454)
        for i, subject in enumerate(ids):
            rows = group.loc[group.subject_id.eq(subject)]
            block = scans.loc[scans.subject_id.eq(subject)].set_index("session")
            for field, sessions in [("pre", rows.baseline_session), ("post", rows.followup_session)]:
                expected = block.loc[list(sessions), [f"whole_brain__{m}" for m in MODES]].mean().to_numpy()
                np.testing.assert_allclose(value[field][i].mean(axis=1), expected, atol=1e-12, rtol=1e-10)
            np.testing.assert_allclose(value["fd"][i], rows.delta_fd.mean(), atol=1e-12)
            assert value["active"][i] == (rows.arm.iloc[0] == "active")
        pp[dataset] = value
    return pp


def regional_inference(pp, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    aa, sud = pp["aachen_itbs"], pp["sudmex_tms"]
    na, ns = len(aa["ids"]), len(sud["ids"])
    np.savez_compressed(out / "regional_modes.npz", **{
        f"{d}_{field}": values[field] for d, values in pp.items() for field in ["pre", "post", "fd", "active", "ids"]})
    aa_template = (aa["post"] - aa["pre"]).mean(0)
    sud_template = np.empty((6, 454))
    for m in range(6):
        for roi in range(454):
            x = np.column_stack([np.ones(ns), sud["active"], sud["pre"][:, m, roi], sud["fd"]])
            assert np.linalg.matrix_rank(x) == x.shape[1]
            sud_template[m, roi] = np.linalg.lstsq(x, sud["post"][:, m, roi], rcond=None)[0][1]
    results, templates, scores, loo, correlations = [], [], [], [], []
    for branch, features in [("all_regions", np.arange(454)), ("cortical_nodes", np.arange(400))]:
        at, st = aa_template[:, features].copy(), sud_template[:, features].copy()
        for matrix in [at, st]:
            matrix -= matrix.mean(1, keepdims=True)
            matrix /= np.abs(matrix).sum(1, keepdims=True)
            assert np.isfinite(matrix).all()
        for donor, template in [("aachen", at), ("sudmex", st)]:
            for k, mode in enumerate(MODES):
                templates.extend(dict(branch=branch, donor=donor, mode=mode, roi=int(roi) + 1, weight=float(weight))
                                 for roi, weight in zip(features, template[k]))
        apre = np.einsum("nmr,mr->nm", aa["pre"][:, :, features], st)
        apost = np.einsum("nmr,mr->nm", aa["post"][:, :, features], st)
        spre = np.einsum("nmr,mr->nm", sud["pre"][:, :, features], at)
        spost = np.einsum("nmr,mr->nm", sud["post"][:, :, features], at)
        changes = apost - apre
        ap, amax = exact_max_t(changes)
        for k, mode in enumerate(MODES):
            correlations.append(dict(branch=branch, mode=mode, pearson_r=stats.pearsonr(at[k], st[k]).statistic,
                                     spearman_rho=stats.spearmanr(at[k], st[k]).statistic,
                                     inference="Descriptive only; no parcel-wise p value."))
            results.append(dict(branch=branch, dataset="aachen_itbs", mode=mode, model="paired_external_template",
                                p=ap[k], p_maxT6=amax[k], **regional_original.paired_fit(changes[:, k])))
            for model, adjust in [("baseline_motion_adjusted", True), ("baseline_adjusted", False)]:
                results.append(dict(branch=branch, dataset="sudmex_tms", mode=mode, model=model,
                                    n_active=int(sud["active"].sum()), n_sham=int((~sud["active"]).sum()),
                                    **regional_original.hc3(spost[:, k], spre[:, k], sud["active"], sud["fd"], adjust)))
            for dataset, info, before, after in [("aachen_itbs", aa, apre, apost), ("sudmex_tms", sud, spre, spost)]:
                for i, subject in enumerate(info["ids"]):
                    scores.append(dict(branch=branch, dataset=dataset, subject_id=subject, mode=mode,
                                       pre_score=before[i, k], post_score=after[i, k], delta=after[i, k] - before[i, k]))
        if branch == "all_regions":
            # Omit one receiving participant, retaining the already external donor template.
            for omit in range(na):
                x = np.delete(changes, omit, axis=0)
                ps, pm = exact_max_t(x)
                for k, mode in enumerate(MODES):
                    loo.append(dict(dataset="aachen_itbs", omission=f"omit_{omit + 1:02d}",
                                    subject_omitted=aa["ids"][omit], mode=mode, estimate=x[:, k].mean(), p=ps[k], p_maxT6=pm[k]))
            for omit in range(ns):
                keep = np.arange(ns) != omit
                for k, mode in enumerate(MODES):
                    loo.append(dict(dataset="sudmex_tms", omission=f"omit_{omit + 1:02d}", subject_omitted=sud["ids"][omit], mode=mode,
                                    **regional_original.hc3(spost[keep, k], spre[keep, k], sud["active"][keep], sud["fd"][keep])))
    result = pd.DataFrame(results)
    result["p_holm6"] = result.groupby(["branch", "dataset", "model"]).p.transform(lambda p: multipletests(p, method="holm")[1])
    result["p_holm12"] = result.groupby(["dataset", "model"]).p.transform(lambda p: multipletests(p, method="holm")[1])
    result.to_csv(out / "external_template_tests.csv", index=False)
    conjunction = []
    for branch in ["all_regions", "cortical_nodes"]:
        for mode in MODES:
            a = result.loc[result.branch.eq(branch) & result["mode"].eq(mode) & result.dataset.eq("aachen_itbs")].iloc[0]
            b = result.loc[result.branch.eq(branch) & result["mode"].eq(mode) & result.dataset.eq("sudmex_tms") & result.model.eq("baseline_motion_adjusted")].iloc[0]
            positive = a.estimate > 0 and b.estimate > 0
            conjunction.append(dict(branch=branch, mode=mode, both_positive=bool(positive),
                                    aachen_estimate=a.estimate, sudmex_estimate=b.estimate,
                                    p_aachen=a.p, p_sudmex=b.p, p=max(a.p, b.p) if positive else 1.))
    cc = pd.DataFrame(conjunction)
    cc["p_holm6"] = cc.groupby("branch").p.transform(lambda p: multipletests(p, method="holm")[1])
    cc["p_holm12"] = multipletests(cc.p, method="holm")[1]
    cc.to_csv(out / "cross_dataset_conjunction.csv", index=False)
    outputs = {"pattern_correlations": pd.DataFrame(correlations), "external_template_weights": pd.DataFrame(templates),
               "participant_projected_scores": pd.DataFrame(scores), "leave_one_out": pd.DataFrame(loo)}
    outputs["leave_one_out"]["p_holm6"] = outputs["leave_one_out"].groupby(["dataset", "omission"]).p.transform(
        lambda p: multipletests(p, method="holm")[1])
    for name, data in outputs.items():
        data.to_csv(out / f"{name}.csv", index=False)
    return result, cc, outputs


def load_audit(path: Path, plan: pd.DataFrame):
    receipt = read_json(path / "audit_receipt.json")
    if (receipt.get("state") != "complete" or not receipt.get("all_193_windows_verified")
            or receipt.get("changed_source_files")):
        raise ValueError("FD audit not complete; refusing QC inference")
    for name, digest in receipt["output_sha256"].items():
        candidate = path / name
        if not candidate.exists():
            candidate = ROOT / name
        if not candidate.is_file() or hashlib.sha256(candidate.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Audit output hash mismatch: {name}")
    for dataset in ["aachen_itbs", "sudmex_tms"]:
        bounds = receipt["registered_thresholds"][dataset]
        assert bounds["mean_fd"] == .5 and bounds["peak_fd"] == 3. and bounds["min_timepoints"] == 100
    strict_bounds = receipt["strict_sensitivity"]
    assert strict_bounds["mean_fd_max"] == .3 and strict_bounds["fraction_fd_above_0_5_max"] == .2
    audit = read_csv(path / "scan_audit.csv")
    read_csv(path / "pair_audit.csv")
    subjects = read_json(path / "cohort_subjects.json")
    required = {"dataset", "subject_id", "session", "arm", "prefix", "retained_mean_fd",
                "registered_pass", "strict_pass", "bounds_verified", "legacy_displacement_pass"}
    if not required.issubset(audit.columns) or audit.duplicated(["dataset", "subject_id", "session"]).any():
        raise ValueError("Invalid audit schema/scan identity")
    for col in ["registered_pass", "strict_pass", "bounds_verified", "legacy_displacement_pass"]:
        audit[col] = bool_column(audit[col])
    if not audit.bounds_verified.all() or not np.isfinite(audit.retained_mean_fd).all():
        raise ValueError("FD bounds/values not verified for every scan")
    assert (~audit.strict_pass | audit.registered_pass).all()
    assert (~audit.registered_pass | audit.legacy_displacement_pass).all()
    lookup = audit.set_index(["dataset", "subject_id", "session"])
    cohort_counts = []
    for branch in ["registered", "strict"]:
        if branch not in subjects:
            raise ValueError(f"Missing branch in audited whitelist: {branch}")
        for dataset in ["aachen_itbs", "sudmex_tms"]:
            expected = []
            block = plan.loc[plan.dataset.eq(dataset)]
            for subject in sorted(block.subject_id.unique()):
                rows = block.loc[block.subject_id.eq(subject)]
                sessions = set(rows.baseline_session) | set(rows.followup_session)
                assert len(sessions) == (4 if dataset == "aachen_itbs" else 2)
                checks = [lookup.loc[(dataset, subject, session), branch + "_pass"] for session in sessions]
                if all(checks):
                    expected.append(subject)
            assert sorted(subjects[branch][dataset]) == expected, (branch, dataset, "audit cohort mismatch")
            selected = block.loc[block.subject_id.isin(expected)].drop_duplicates("subject_id")
            cohort_counts.append(dict(qc_policy=branch, dataset=dataset, participants=len(selected),
                                      active=int(selected.arm.eq("active").sum()), sham=int(selected.arm.eq("sham").sum()),
                                      iTBS=int(selected.arm.eq("iTBS").sum())))
    return audit, subjects, pd.DataFrame(cohort_counts), receipt


def filtered_inputs(scans, plan, regional, audit, subject_whitelist):
    ids = {subject for group in subject_whitelist.values() for subject in group}
    kept_plan = plan.loc[plan.subject_id.isin(ids)].copy()
    lookup = audit.set_index(["dataset", "subject_id", "session"])
    kept_scans = scans.loc[scans.subject_id.isin(ids)].copy()
    kept_scans["mean_fd_legacy"] = kept_scans.mean_fd
    kept_scans["mean_fd"] = [lookup.loc[(r.dataset, r.subject_id, r.session), "retained_mean_fd"]
                             for r in kept_scans.itertuples()]
    kept_plan["delta_fd"] = [
        lookup.loc[(r.dataset, r.subject_id, r.followup_session), "retained_mean_fd"] -
        lookup.loc[(r.dataset, r.subject_id, r.baseline_session), "retained_mean_fd"] for r in kept_plan.itertuples()]
    pp = {}
    for dataset, old in regional.items():
        keep = np.array([subject in subject_whitelist[dataset] for subject in old["ids"]])
        pp[dataset] = {field: old[field][keep].copy() for field in ["pre", "post", "fd", "active"]}
        pp[dataset]["ids"] = [subject for subject, yes in zip(old["ids"], keep) if yes]
        pp[dataset]["fd"] = np.array([kept_plan.loc[kept_plan.subject_id.eq(subject), "delta_fd"].mean()
                                        for subject in pp[dataset]["ids"]])
        assert len(pp[dataset]["ids"]) == len(subject_whitelist[dataset])
    return kept_scans, kept_plan, pp


def summarize(policy, family, results, cc):
    if "reference_seed" in results:
        results = results.loc[results.reference_seed.eq("combined")]
        cc = cc.loc[cc.reference_seed.eq("combined")]
    primary = results.loc[~results.model.eq("baseline_adjusted")]
    scope = "branch" if family == "regional_templates" else "scope"
    corrected_all = "p_holm12" if family == "regional_templates" else "p_holm42"
    return dict(qc_policy=policy, analysis=family, tests_per_model=int(primary.groupby("dataset").size().iloc[0]),
                aachen_n=int(primary.loc[primary.dataset.eq("aachen_itbs"), "n"].iloc[0]),
                sudmex_n=int(primary.loc[primary.dataset.eq("sudmex_tms"), "n"].iloc[0]),
                aachen_min_p=float(primary.loc[primary.dataset.eq("aachen_itbs"), "p"].min()),
                sudmex_min_p=float(primary.loc[primary.dataset.eq("sudmex_tms"), "p"].min()),
                single_cohort_holm6_passes=int((primary.p_holm6 < .05).sum()),
                single_cohort_all_scope_passes=int((primary[corrected_all] < .05).sum()),
                aachen_maxT6_passes=int((primary.loc[primary.dataset.eq("aachen_itbs"), "p_maxT6"] < .05).sum()),
                aachen_maxT42_passes=int((primary.loc[primary.dataset.eq("aachen_itbs"), "p_maxT42"] < .05).sum())
                if "p_maxT42" in primary else None,
                shared_uncorrected_passes=int((cc.p < .05).sum()),
                shared_holm6_passes=int((cc.p_holm6 < .05).sum()),
                shared_all_scope_passes=int((cc[corrected_all] < .05).sum()),
                minimum_shared_p=float(cc.p.min()), minimum_shared_holm6=float(cc.p_holm6.min()),
                primary_scope=scope)


def changes_from_original(policy, family, current, previous, conjunction=False):
    """Full test inventory, with estimates and corrected p values before/after QC."""
    keys = (["branch", "mode"] if family == "regional_templates" else ["metric"])
    if not conjunction:
        keys += ["dataset", "model"]
    if "reference_seed" in current:
        keys += ["reference_seed"]
    columns = [c for c in ["n", "n_active", "n_sham", "estimate", "ci_low", "ci_high", "p", "p_t",
                           "p_maxT6", "p_maxT42", "p_holm6", "p_holm12", "p_holm42", "same_direction",
                           "both_positive", "aachen_estimate", "sudmex_estimate", "p_aachen", "p_sudmex"]
               if c in current and c in previous]
    before = previous[keys + columns].rename(columns={c: "original_" + c for c in columns})
    after = current[keys + columns].rename(columns={c: "qc_" + c for c in columns})
    merged = before.merge(after, on=keys, how="outer", validate="one_to_one", indicator=True)
    assert merged._merge.eq("both").all()
    return merged.drop(columns="_merge").assign(qc_policy=policy, analysis=family)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--verify-only", action="store_true", help="Only reproduce old complete-sample tables; no QC inference.")
    args = parser.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(DEFAULT_OUT.resolve()):
        raise ValueError("Output must remain in the assigned modes directory")
    out.mkdir(parents=True, exist_ok=True)
    track(Path(__file__))
    for source in [Path(original.__file__), Path(regional_original.__file__),
                   ROOT / "setup/analyse_tms_familywise.py", ROOT / "setup/test_corrected_gwt_expectations.py",
                   ROOT / "setup/normalize_pid_processing_modes.py"]:
        track(source)
    for folder, manifest in [(RAW, "source_sha256.json"), (NORMALIZED, "run_source_sha256.json"), (REGIONAL, "source_sha256.json")]:
        track(folder / manifest)
        verification = read_json(folder / "verification.json")
        assert verification["passed"]
        track(folder / "analysis_specification.json")
    assert list(read_json(RAW / "scope_pair_counts.json")) == SCOPES
    plan = read_csv(CORE / "raw_paired_changes.csv")
    raw = read_csv(RAW / "scan_modes_and_atoms.csv")
    normalized = read_csv(NORMALIZED / "scan_modes_all_seeds.csv", dtype={"reference_seed": str})
    assert len(raw) == 140 and len(normalized) == 420 and set(normalized.reference_seed) == set(SEEDS)
    regional = load_regional(plan, raw)
    if not args.verify_only:
        audit, whitelists, cohort_counts, audit_receipt = load_audit(args.audit, plan)
        cohort_counts.to_csv(out / "cohort_counts.csv", index=False)
    spec = dict(created_utc=datetime.now(timezone.utc).isoformat(),
                purpose="Rerun all existing six-mode, seven-scope and external regional-template tests after retained-frame FD QC.",
                modes=original.MODES, scopes=SCOPES,
                audit=str(args.audit.relative_to(ROOT)),
                qc_registered="Existing intended rules: retained mean FD <=0.5 mm, peak FD <=3 mm, >=100 retained frames, plus original displacement screens.",
                qc_strict="Registered rules plus retained mean FD <=0.30 mm and fraction retained FD >0.5 mm <=0.20.",
                pairing="All four Aachen scans and both matched visits per person required; mean of two changes. SUDMEX paired t0/t1 active and sham.",
                motion_covariate="Post minus baseline mean FD on the exact retained-frame interval, replacing legacy mean FD only in QC analyses.",
                models="Aachen exact paired sign flips; SUDMEX post ~ active + baseline + deltaFD, HC3/t. Baseline-only SUDMEX sensitivity.",
                correction="Mode Holm6 per scope and Holm42 across scopes; Aachen maxT6/maxT42. Regional Holm6 per branch/Holm12 across all454 and cortical400; Aachen maxT6.",
                conjunction="max(two two-sided p), p=1 for discordant mode signs; regional requires both projected effects positive. Existing correction families retained.",
                normalized_reference="FIXED legacy baseline calibration and two seeds. No reference rebuilding with new QC. These are fixed-reference QC sensitivities.",
                regional_templates="Refit donor-cohort regional effects independently within each QC cohort; receiver uses only other-cohort template. LOO removes receiver only and retains external template.",
                limitations="Exploratory analyses. QC changes eligibility and precision; no claim of stability from nonsignificance or neuronal causation from predictive mode names.")
    save_json(out / "analysis_specification.json", spec)
    reference = out / "full_sample_reproduction"
    r, c, p = mode_inference(raw, plan, reference / "raw")
    compare_table(r, RAW / "mode_treatment_results.csv", ["dataset", "model", "metric"], "raw treatment")
    compare_table(c, RAW / "cross_dataset_conjunction.csv", ["metric"], "raw conjunction")
    compare_table(p, RAW / "participant_changes.csv", ["dataset", "subject_id", "comparison"], "raw paired changes")
    nr, nc, stability = normalized_inference(normalized, plan, reference / "normalized_fixed_reference")
    compare_table(nr, NORMALIZED / "mode_treatment_results_all_seeds.csv", ["reference_seed", "dataset", "model", "metric"], "normalized all-seed treatment")
    compare_table(nc, NORMALIZED / "cross_dataset_conjunction_all_seeds.csv", ["reference_seed", "metric"], "normalized all-seed conjunction")
    compare_table(stability, NORMALIZED / "monte_carlo_stability.csv", ["dataset", "kind", "metric"], "normalized seed stability")
    rr, rc, auxiliary = regional_inference(regional, reference / "regional_templates")
    compare_table(rr, REGIONAL / "external_template_tests.csv", ["branch", "dataset", "model", "mode"], "regional tests")
    compare_table(rc, REGIONAL / "cross_dataset_conjunction.csv", ["branch", "mode"], "regional conjunction")
    for name, keys in [("pattern_correlations", ["branch", "mode"]), ("external_template_weights", ["branch", "donor", "mode", "roi"]),
                       ("participant_projected_scores", ["branch", "dataset", "subject_id", "mode"]), ("leave_one_out", ["dataset", "omission", "mode"])]:
        compare_table(auxiliary[name], REGIONAL / f"{name}.csv", keys, "regional " + name)
    pd.DataFrame(COMPARISONS).to_csv(out / "full_sample_reproduction_checks.csv", index=False)
    print("Full-sample reproduction passed for all original result and diagnostic tables", flush=True)
    summaries, result_changes, conjunction_changes = [], [], []
    previous = {"raw_modes": (r, c), "normalized_modes_fixed_reference": (nr, nc), "regional_templates": (rr, rc)}
    if not args.verify_only:
        for policy in ["registered", "strict"]:
            scans, kept_plan, pp = filtered_inputs(raw, plan, regional, audit, whitelists[policy])
            norm, norm_plan, _ = filtered_inputs(normalized, plan, regional, audit, whitelists[policy])
            pd.testing.assert_frame_equal(kept_plan, norm_plan)
            destination = out / policy
            destination.mkdir(exist_ok=True)
            scans.to_csv(destination / "scan_modes_and_atoms.csv", index=False)
            norm.to_csv(destination / "scan_modes_all_seeds_fixed_reference.csv", index=False)
            kept_plan[["dataset", "subject_id", "arm", "comparison", "baseline_session", "followup_session", "delta_fd"]].to_csv(destination / "eligible_pair_plan.csv", index=False)
            r, c, _ = mode_inference(scans, kept_plan, destination / "raw")
            summaries.append(summarize(policy, "raw_modes", r, c))
            nr, nc, _ = normalized_inference(norm, kept_plan, destination / "normalized_fixed_reference")
            summaries.append(summarize(policy, "normalized_modes_fixed_reference", nr, nc))
            rr, rc, _ = regional_inference(pp, destination / "regional_templates")
            summaries.append(summarize(policy, "regional_templates", rr, rc))
            for family, current_results, current_cc in [("raw_modes", r, c),
                    ("normalized_modes_fixed_reference", nr, nc), ("regional_templates", rr, rc)]:
                before_results, before_cc = previous[family]
                result_changes.append(changes_from_original(policy, family, current_results, before_results))
                conjunction_changes.append(changes_from_original(policy, family, current_cc, before_cc, conjunction=True))
            print(f"Completed policy={policy}; Aachen n={len(pp['aachen_itbs']['ids'])}; SUDMEX n={len(pp['sudmex_tms']['ids'])}", flush=True)
    summary = pd.DataFrame(summaries)
    summary.to_csv(out / "headline_summary.csv", index=False)
    if result_changes:
        pd.concat(result_changes, ignore_index=True).to_csv(out / "full_sample_vs_qc_results.csv", index=False)
        pd.concat(conjunction_changes, ignore_index=True).to_csv(out / "full_sample_vs_qc_conjunctions.csv", index=False)
    save_json(out / "source_sha256.json", SOURCES)
    for relative, digest in SOURCES.items():
        if hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() != digest:
            raise AssertionError(f"Source changed during rerun: {relative}")
    save_json(out / "verification.json", dict(passed=True, state="reproduction_only" if args.verify_only else "complete",
               original_tables_reproduced=len(COMPARISONS), source_files_hashed=len(SOURCES),
               sources_unchanged=True, normalized_calibration_rebuilt=False,
               regional_donor_templates_refitted_per_qc_policy=not args.verify_only,
               policies=[] if args.verify_only else ["registered", "strict"],
               full_sample_max_absolute_difference=max(c["max_absolute_numeric_difference"] for c in COMPARISONS)))
    print(summary.to_string(index=False) if not summary.empty else "Reproduction-only verification complete", flush=True)
    print("COMPLETE", out, flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
