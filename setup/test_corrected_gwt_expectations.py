"""Test the four source-grounded GWT expectations from corrected TMS tables.

This separate report reads saved participant changes. It does not change the
core analyses or processing. Both raw and normalized branches use the same
four endpoints, participant pairing, and two-sided inference.
"""
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
METRICS = ("gateways_FPN", "gateways_DAN", "mean_redundancy_pc", "mean_synergy_pc")
SUBJECTS = [f"aachen_itbs_{i:02d}" for i in range(1, 17) if i != 2]
INVENTORY = ROOT / "output/corrected_tms_analysis_review/gwt_expectation_inventory.md"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def holm(p):
    """Keep unestimable endpoints in the declared correction burden."""
    p = np.asarray(p, float)
    finite = np.isfinite(p)
    values = np.where(finite, p, 1.)
    order = np.argsort(values)
    corrected = np.minimum(1., np.maximum.accumulate(values[order] * np.arange(len(p), 0, -1)))
    result = np.empty_like(values)
    result[order] = corrected
    result[~finite] = np.nan
    return result


def exact_signflip(values):
    """Enumerate signs and recompute SD on signed observations independently."""
    x = np.asarray(values, float)
    if x.ndim != 2 or not 3 <= len(x) <= 20 or not np.isfinite(x).all():
        raise ValueError("Exact sign flips require complete participant-by-endpoint data, 3 <= n <= 20")
    usable = x.std(axis=0, ddof=1) > 1e-13
    raw, adjusted = np.full(x.shape[1], np.nan), np.full(x.shape[1], np.nan)
    if not usable.any():
        return raw, adjusted
    data = x[:, usable]
    observed = np.abs(data.mean(0) / (data.std(axis=0, ddof=1) / np.sqrt(len(x))))
    signs = np.asarray(list(product((-1., 1.), repeat=len(x))))
    counts = np.zeros(data.shape[1], dtype=int)
    maxima = np.zeros(data.shape[1], dtype=int)
    for first in range(0, len(signs), 256):
        signed = signs[first:first + 256, :, None] * data[None, :, :]
        mean, sd = signed.mean(axis=1), signed.std(axis=1, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            null_t = np.abs(mean / (sd / np.sqrt(len(x))))
        null_t[(sd == 0) & (mean == 0)] = 0.
        if np.isnan(null_t).any():
            raise ValueError("Undefined sign-flip statistic")
        counts += (null_t >= observed - 1e-10).sum(0)
        maxima += (null_t.max(axis=1)[:, None] >= observed - 1e-10).sum(0)
    raw[usable], adjusted[usable] = counts / len(signs), maxima / len(signs)
    return raw, adjusted


def direct_hc3(frame, outcome, predictors, term):
    """OLS + explicit HC3 sandwich, with t(df=n-p) reference."""
    data = frame[[outcome, *predictors]].replace([np.inf, -np.inf], np.nan).dropna()
    keep = [c for c in predictors if c == "intercept" or data[c].nunique() > 1]
    if term not in keep:
        return {"n": len(data), "status": "constant_predictor", "p": np.nan}
    x, y = data[keep].to_numpy(float), data[outcome].to_numpy(float)
    if len(y) <= x.shape[1] + 2 or np.linalg.matrix_rank(x) != x.shape[1]:
        return {"n": len(data), "status": "rank_deficient_or_small_sample", "p": np.nan}
    pinv = np.linalg.pinv(x)
    beta = pinv @ y
    leverage = np.sum(x * pinv.T, axis=1)
    if np.any(1 - leverage <= 1e-10):
        return {"n": len(data), "status": "unit_leverage", "p": np.nan}
    residual = (y - x @ beta) / (1 - leverage)
    cov = (pinv * residual[None, :]) @ (pinv * residual[None, :]).T
    ix, df = keep.index(term), len(y) - x.shape[1]
    se = np.sqrt(cov[ix, ix])
    if not np.isfinite(se) or se <= 0:
        return {"n": len(data), "status": "invalid_standard_error", "p": np.nan}
    t = beta[ix] / se
    width = stats.t.ppf(.975, df) * se
    return {"n": len(data), "estimate": beta[ix], "ci_low": beta[ix] - width,
            "ci_high": beta[ix] + width, "se": se, "t": t, "df": df,
            "p": 2 * stats.t.sf(abs(t), df), "max_leverage": float(leverage.max()), "status": "ok"}


def self_check():
    """Small independent fixtures for the numerical helpers only."""
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(2026112604)
    x = rng.normal(size=(6, 4)) + [.1, -.2, .3, 0.]
    null = np.asarray([np.abs(stats.ttest_1samp(x * np.asarray(s)[:, None], 0, axis=0).statistic)
                       for s in product((-1., 1.), repeat=len(x))])
    obs = abs(stats.ttest_1samp(x, 0, axis=0).statistic)
    p, pc = exact_signflip(x)
    np.testing.assert_array_equal(p, (null >= obs - 1e-10).mean(0))
    np.testing.assert_array_equal(pc, (null.max(1)[:, None] >= obs - 1e-10).mean(0))
    np.testing.assert_allclose(holm([.01, .03, .04, .8]), multipletests([.01, .03, .04, .8], method="holm")[1])
    frame = pd.DataFrame({"y": rng.normal(size=40), "intercept": 1., "active": np.arange(40) % 2,
                          "baseline": rng.normal(size=40), "delta_fd": rng.normal(size=40)})
    columns = ["intercept", "active", "baseline", "delta_fd"]
    a = direct_hc3(frame, "y", columns, "active")
    b = sm.OLS(frame.y, frame[columns]).fit(cov_type="HC3", use_t=True)
    np.testing.assert_allclose([a["estimate"], a["se"], a["p"], a["ci_low"], a["ci_high"]],
                               [b.params.active, b.bse.active, b.pvalues.active, *b.conf_int().loc["active"]], atol=1e-12)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--branch", choices=["raw", "normalized"], default="raw")
    args = parser.parse_args()
    run = args.run.resolve()
    if not run.is_relative_to(ROOT.resolve()):
        raise ValueError("Reports must stay inside TMS_GWT")
    self_check()
    source = run / f"{args.branch}_paired_changes.csv"
    pairs = pd.read_csv(source)
    key = ["dataset", "subject_id", "comparison"]
    if pairs.duplicated(key).any():
        raise ValueError("Duplicate participant contrasts")
    aachen = pairs[pairs.dataset.eq("aachen_itbs")]
    if sorted(aachen.subject_id.unique()) != SUBJECTS or len(aachen) != 30:
        raise ValueError("Expected 15 paired Aachen participants")
    columns = [m + "_delta" for m in METRICS]
    pooled = []
    for subject in SUBJECTS:
        rows = aachen[aachen.subject_id.eq(subject)].set_index("comparison")
        if set(rows.index) != {"post90_minus_pre90", "post120_minus_pre120"}:
            raise ValueError("Both intensity visits are required")
        for contrast, pre, post in [("post90_minus_pre90", "pre90", "post90"),
                                    ("post120_minus_pre120", "pre120", "post120")]:
            row = rows.loc[contrast]
            if row.baseline_session != pre or row.followup_session != post:
                raise ValueError("Incorrect Aachen session labels")
        x = rows[columns].to_numpy(float)
        if not np.isfinite(x).all():
            raise ValueError("Incomplete endpoint pairing")
        pooled.append(x.mean(axis=0))
    pooled = np.asarray(pooled)
    pflip, pmax = exact_signflip(pooled)
    records = []
    for k, metric in enumerate(METRICS):
        x = pooled[:, k]
        se = x.std(ddof=1) / np.sqrt(len(x))
        if x.std(ddof=1) <= 1e-13:
            records.append({"dataset": "aachen_itbs", "metric": metric, "n": 15,
                            "model": "participant mean of two post-pre contrasts", "estimate": x.mean(),
                            "p": np.nan, "p_sign_flip": np.nan, "p_maxT4": np.nan,
                            "sign_vectors": 32768, "status": "constant_endpoint_not_tested"})
            continue
        t = x.mean() / se
        width = stats.t.ppf(.975, 14) * se
        records.append({"dataset": "aachen_itbs", "metric": metric, "n": 15,
                        "model": "participant mean of two post-pre contrasts", "estimate": x.mean(),
                        "ci_low": x.mean() - width, "ci_high": x.mean() + width,
                        "se": se, "t": t, "df": 14, "p": 2 * stats.t.sf(abs(t), 14),
                        "p_sign_flip": pflip[k], "p_maxT4": pmax[k], "sign_vectors": 32768, "status": "ok"})
    sud = pairs[pairs.dataset.eq("sudmex_tms")].copy()
    if len(sud) != 40 or sud.subject_id.duplicated().any() or sud.arm.value_counts().to_dict() != {"active": 22, "sham": 18}:
        raise ValueError("Expected 40 acute randomized SUDMEX pairs")
    if not (sud.baseline_session.eq("t0") & sud.followup_session.eq("t1") & sud.comparison.eq("t1_minus_t0")).all():
        raise ValueError("Only acute t0/t1 pairs are permitted")
    sud["active"], sud["intercept"] = sud.arm.eq("active").astype(float), 1.
    for metric in METRICS:
        predictors = ["intercept", "active", metric + "_baseline", "delta_fd"]
        required = [metric + "_followup", metric + "_delta", *predictors]
        if not np.isfinite(sud[required].to_numpy(float)).all():
            raise ValueError("Expected complete SUDMEX endpoint and covariate data")
        np.testing.assert_allclose(sud[metric + "_followup"] - sud[metric + "_baseline"],
                                   sud[metric + "_delta"], atol=1e-12, rtol=1e-10)
        result = direct_hc3(sud, metric + "_followup", predictors, "active")
        records.append({"dataset": "sudmex_tms", "metric": metric, "n_active": 22, "n_sham": 18,
                        "model": "followup ~ active + baseline outcome + delta FD; HC3/t", **result})
    result = pd.DataFrame(records)
    result["p_holm4"] = result.groupby("dataset").p.transform(holm)
    result["branch"], result["expected_direction"], result["test_sidedness"] = args.branch, "increase", "two-sided"
    result["family_size"] = 4
    result.to_csv(run / f"gwt_expectations_{args.branch}.csv", index=False)
    receipt = {"created_utc": datetime.now(timezone.utc).isoformat(), "metrics": METRICS,
               "branch": args.branch, "output_rows": len(result), "aachen_participants": 15,
               "sudmex_participants": 40, "sign_flip_vectors": 32768,
               "input_sha256": {str(p): sha256(p) for p in [source, INVENTORY, Path(__file__)]},
               "inference": "Holm4 of two-sided t p-values per dataset; Aachen also exact joint sign-flip maxT4. CIs are marginal. Sign-flip inference assumes independent participants and joint sign symmetry. HC3 uses a small-sample t approximation; motion-adjusted SUDMEX is a sensitivity estimand.",
               "scope": "Expectation family defined before corrected effect review; not a preregistration. No treatment source metrics or core analysis files changed."}
    (INVENTORY.parent / f"gwt_expectations_{args.branch}_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(result[["dataset", "metric", "estimate", "p", "p_holm4", "p_maxT4"]].to_string(index=False))


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
