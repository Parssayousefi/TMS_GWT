"""Independently verify the DAN gateway finding shown in corrected Figure 2.

Read scan-level outputs; preserve all analysis inputs and the existing report.
Reconstruct participant changes, t intervals, exhaustive sign flips and the
four-outcome family. Make a separately labelled figure for interpretation.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
METRICS = ["gateways_FPN", "gateways_DAN", "mean_redundancy_pc", "mean_synergy_pc"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    out = ROOT / "output/dan_gateway_audit" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True)
    inputs = [Path(__file__)]
    for branch in ["raw", "normalized"]:
        inputs += [RUN / f"{branch}_scan_metrics.csv", RUN / f"{branch}_paired_changes.csv",
                   RUN / f"{branch}_aachen_participant_changes.csv", RUN / f"gwt_expectations_{branch}.csv"]
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    results, subject_rows, visits, sud_results = [], [], [], []
    for branch in ["raw", "normalized"]:
        scans = pd.read_csv(RUN / f"{branch}_scan_metrics.csv")
        pairs = pd.read_csv(RUN / f"{branch}_paired_changes.csv")
        archived = pd.read_csv(RUN / f"gwt_expectations_{branch}.csv").set_index(["dataset", "metric"])
        aa = scans.loc[scans.dataset.eq("aachen_itbs")].set_index(["subject_id", "session"])
        ids = sorted(aa.index.get_level_values(0).unique())
        assert len(ids) == 15 and len(aa) == 60
        pre = np.stack([aa.loc[(s, v), METRICS].to_numpy(float) for s in ids for v in ["pre90", "pre120"]]).reshape(15, 2, 4)
        post = np.stack([aa.loc[(s, v), METRICS].to_numpy(float) for s in ids for v in ["post90", "post120"]]).reshape(15, 2, 4)
        change = (post - pre).mean(axis=1)
        saved = pd.read_csv(RUN / f"{branch}_aachen_participant_changes.csv").set_index("subject_id").loc[ids]
        np.testing.assert_allclose(change, saved[[m + "_delta" for m in METRICS]], atol=1e-12)
        ttest = stats.ttest_1samp(change, 0, axis=0)
        mean, se = change.mean(0), stats.sem(change, axis=0)
        ci = np.stack([mean - stats.t.ppf(.975, 14) * se, mean + stats.t.ppf(.975, 14) * se])
        # Independent exact sign enumeration via bits; use sufficient statistics
        # instead of the published helper's signed-array standard deviations.
        signs = 2 * (((np.arange(2**15)[:, None] >> np.arange(15)) & 1).astype(float)) - 1
        null_mean = signs @ change / 15
        null_se = np.sqrt(np.maximum((np.sum(change**2, axis=0) - 15*null_mean**2) / (14*15), 0))
        null_t = np.abs(np.divide(null_mean, null_se, out=np.zeros_like(null_mean), where=null_se > 0))
        nominal_flip = (null_t >= abs(ttest.statistic)[None, :] - 1e-10).mean(0)
        max_t = (null_t.max(1)[:, None] >= abs(ttest.statistic)[None, :] - 1e-10).mean(0)
        holm = multipletests(ttest.pvalue, method="holm")[1]
        for j, m in enumerate(METRICS):
            row = dict(branch=branch, dataset="aachen_itbs", metric=m, n=15,
                       pre_mean=pre[:, :, j].mean(), post_mean=post[:, :, j].mean(),
                       estimate=mean[j], ci_low=ci[0, j], ci_high=ci[1, j], p=ttest.pvalue[j],
                       p_sign_flip=nominal_flip[j], p_holm4=holm[j], p_maxT4=max_t[j],
                       n_decrease=int((change[:, j] < 0).sum()), n_zero=int((change[:, j] == 0).sum()),
                       n_increase=int((change[:, j] > 0).sum()))
            old = archived.loc[("aachen_itbs", m)]
            for key in ["estimate", "ci_low", "ci_high", "p", "p_sign_flip", "p_holm4", "p_maxT4"]:
                np.testing.assert_allclose(row[key], old[key], atol=1e-12)
            results.append(row)
        for i, s in enumerate(ids):
            subject_rows.append(dict(branch=branch, subject_id=s, pre90=pre[i,0,1], post90=post[i,0,1],
                                     pre120=pre[i,1,1], post120=post[i,1,1], mean_change=change[i,1]))
        for j, label in enumerate(["90%", "120%"]):
            x = post[:, j, 1] - pre[:, j, 1]
            visits.append(dict(branch=branch, visit=label, pre_mean=pre[:,j,1].mean(),
                               post_mean=post[:,j,1].mean(), change_mean=x.mean(),
                               n_decrease=int((x < 0).sum()), n_zero=int((x == 0).sum()), n_increase=int((x > 0).sum())))
        sud = pairs.loc[pairs.dataset.eq("sudmex_tms")].copy()
        assert len(sud) == 40 and sud.arm.value_counts().to_dict() == {"active": 22, "sham": 18}
        for m in METRICS:
            # Rebuild pre/post counts from scans, independently of the paired table.
            source = scans.loc[scans.dataset.eq("sudmex_tms")].set_index(["subject_id", "session"])
            baseline = np.asarray([source.loc[(s, "t0"), m] for s in sud.subject_id], float)
            followup = np.asarray([source.loc[(s, "t1"), m] for s in sud.subject_id], float)
            np.testing.assert_allclose(baseline, sud[m + "_baseline"])
            np.testing.assert_allclose(followup, sud[m + "_followup"])
            x = pd.DataFrame({"constant": 1., "active": sud.arm.eq("active").astype(float).to_numpy(),
                              "baseline": baseline, "delta_fd": sud.delta_fd.to_numpy()})
            fit = sm.OLS(followup, x).fit(cov_type="HC3", use_t=True)
            old = archived.loc[("sudmex_tms", m)]
            np.testing.assert_allclose([fit.params.active, fit.pvalues.active, *fit.conf_int().loc["active"]],
                                      old[["estimate", "p", "ci_low", "ci_high"]].to_numpy(float), atol=1e-12)
            sud_results.append(dict(branch=branch, metric=m, n=40, estimate=fit.params.active,
                                    ci_low=fit.conf_int().loc["active", 0], ci_high=fit.conf_int().loc["active", 1],
                                    p=fit.pvalues.active,
                                    active_change=np.mean((followup-baseline)[x.active.eq(1)]),
                                    sham_change=np.mean((followup-baseline)[x.active.eq(0)])))
    frame = pd.DataFrame(results)
    frame.to_csv(out / "aachen_four_outcome_audit.csv", index=False)
    pd.DataFrame(subject_rows).to_csv(out / "dan_participant_counts.csv", index=False)
    pd.DataFrame(visits).to_csv(out / "dan_visit_descriptives.csv", index=False)
    sud = pd.DataFrame(sud_results)
    sud["p_holm4"] = sud.groupby("branch").p.transform(lambda p: multipletests(p, method="holm")[1])
    sud.to_csv(out / "sudmex_four_outcome_audit.csv", index=False)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    titles = ["FPN gateways", "DAN gateways", "Redundancy participation", "Synergy participation"]
    for ax, m, title in zip(axes.flat, METRICS, titles):
        for y, branch, color in [(1, "raw", "#1c4f65"), (0, "normalized", "#8c5b9d")]:
            r = frame.loc[frame.branch.eq(branch) & frame.metric.eq(m)].iloc[0]
            ax.errorbar(r.estimate, y, xerr=[[r.estimate-r.ci_low], [r.ci_high-r.estimate]], fmt="o", color=color, capsize=4)
            ax.text(.03, .62 if y == 1 else .17,
                    f"t-test p={r.p:.3f}; family-wise maxT p={r.p_maxT4:.3f}",
                    transform=ax.transAxes, fontsize=9, color=color)
        ax.axvline(0, color=".55", lw=1, ls="--")
        ax.set_yticks([0, 1], ["Normalized", "Raw"]); ax.set_ylim(-.65, 1.55)
        ax.set_title(title)
        ax.set_xlabel("Post minus pre (regions)" if m.startswith("gateways") else "Post minus pre (coefficient)")
    fig.suptitle("Aachen: unadjusted 95% intervals; family-wise tests cover all four outcomes", fontsize=12)
    fig.savefig(out / "gwt_expectations_with_inference.png", dpi=160)
    plt.close(fig)
    assert all(sha(ROOT / p) == h for p, h in hashes.items())
    receipt = {"passed": True, "source_sha256": hashes, "n_aachen": 15, "aachen_scans": 60,
               "n_sudmex": 40, "sign_patterns": 32768,
               "verified": ["Counts and changes reconstructed from scan-level tables", "All eight Aachen rows match the published intervals and p-values", "SUDMEX HC3 results match via statsmodels", "Inputs unchanged"],
               "interpretation": "Normalized Aachen DAN reduction is nominal under a t test only; exact unadjusted sign-flip p=.056640625, Holm4 p=.152730, maxT4 p=.162109. SUDMEX does not reproduce a reduction."}
    (out / "verification.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    (out.parent / "LATEST.txt").write_text(out.name + "\n", encoding="utf-8")
    print("OUTPUT", out)
    print(frame.loc[frame.metric.eq("gateways_DAN")].to_string(index=False))
    print(sud.loc[sud.metric.eq("gateways_DAN")].to_string(index=False))
    print(pd.DataFrame(visits).to_string(index=False))
    print("All checks passed")


if __name__ == "__main__":
    main()
