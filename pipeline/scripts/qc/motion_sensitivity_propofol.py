"""Motion sensitivity analysis for the propofol cohort — reviewer-proofing.

THE PROBLEM A REVIEWER WILL RAISE
    In this cohort head motion is confounded with the condition of interest:
    deep sedation has ~2x the mean FD of the awake baseline. Any sedation effect
    could therefore be a motion effect. Removing ART scrubbing (the uniform
    pipeline) removes one line of defence, so the question sharpens: are the
    propofol results motion-driven?

WHAT THIS DOES, per variant (ART vs uniform) and per contrast:
    1. EFFECT          paired t on the per-subject change in each measure
    2. MOTION COUPLING Pearson r between per-subject dFD and dMeasure
                       (dX = drug - awake, so this is the within-subject question,
                       not the between-subject one)
    3. dFD-ADJUSTED    lm(dMeasure ~ dFD); the intercept test is the effect with
                       the motion slope partialled out. Reported with the
                       unadjusted effect so shrinkage is visible.
    4. MOTION-MATCHED  re-test on the subjects whose |dFD| is at or below the
                       median, i.e. the half of the sample where motion barely
                       differs between conditions
    5. CROSS-VARIANT   the two pipelines differ ONLY in motion handling, so
                       agreement between them is itself a motion control

Inputs: results/qc/{dataset}/motion_synergy_{label}/per_scan.csv
        (written by scripts/qc/cambridge_motion_synergy_baseline.py)

Usage:
    python scripts/qc/motion_sensitivity_propofol.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
QC = PROJECT_DIR / "results" / "qc"

VARIANTS = {
    "ART":     (QC / "cambridge_propofol" / "motion_synergy_post_art" / "per_scan.csv"),
    "UNIFORM": (QC / "cambridge_propofol_uniform" / "motion_synergy_uniform_rerun" / "per_scan.csv"),
}
CONTROL = "restawake"
DRUGS = ["restdeep", "restlight", "restrecovery"]
MEASURES = ["synergy_workspace_mean", "redundancy_workspace_mean",
            "phi_r_workspace_mean", "phi_workspace_mean"]
MOTION = "fd_mean"


def paired_frame(df: pd.DataFrame, drug: str, measure: str) -> pd.DataFrame:
    """One row per subject with dMeasure and dFD for this contrast."""
    d = df[df.task == drug][["subject", measure, MOTION]].rename(
        columns={measure: "m_drug", MOTION: "fd_drug"})
    c = df[df.task == CONTROL][["subject", measure, MOTION]].rename(
        columns={measure: "m_ctrl", MOTION: "fd_ctrl"})
    p = d.merge(c, on="subject").dropna()
    p["d_measure"] = p.m_drug - p.m_ctrl
    p["d_fd"] = p.fd_drug - p.fd_ctrl
    return p


def analyse(p: pd.DataFrame) -> dict:
    n = len(p)
    out = {"n": n}
    if n < 4:
        return out
    t, pv = stats.ttest_rel(p.m_drug, p.m_ctrl)
    out["t"], out["p"] = float(t), float(pv)
    out["mean_d"] = float(p.d_measure.mean())
    out["dz"] = float(p.d_measure.mean() / p.d_measure.std(ddof=1))

    r, rp = stats.pearsonr(p.d_fd, p.d_measure)
    out["r_motion"], out["p_motion"] = float(r), float(rp)

    # effect with the dFD slope partialled out: intercept of dMeasure ~ dFD,
    # evaluated at the sample-mean dFD is just the mean; the informative test is
    # the intercept at dFD = 0 (i.e. a hypothetical motion-free subject).
    X = np.column_stack([np.ones(n), p.d_fd.to_numpy()])
    y = p.d_measure.to_numpy()
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - 2
    mse = float(resid @ resid) / dof
    cov = mse * np.linalg.inv(X.T @ X)
    se_b0 = float(np.sqrt(cov[0, 0]))
    out["adj_intercept"] = float(beta[0])
    out["adj_t"] = float(beta[0] / se_b0) if se_b0 > 0 else np.nan
    out["adj_p"] = float(2 * stats.t.sf(abs(out["adj_t"]), dof)) if se_b0 > 0 else np.nan

    # motion-matched half: |dFD| <= median(|dFD|)
    keep = p[p.d_fd.abs() <= p.d_fd.abs().median()]
    out["n_matched"] = len(keep)
    if len(keep) >= 4:
        tm, pm = stats.ttest_rel(keep.m_drug, keep.m_ctrl)
        out["matched_t"], out["matched_p"] = float(tm), float(pm)
        out["matched_mean_d"] = float(keep.d_measure.mean())
        out["matched_mean_dfd"] = float(keep.d_fd.mean())
    out["mean_dfd"] = float(p.d_fd.mean())
    return out


def main() -> int:
    frames = {}
    for name, path in VARIANTS.items():
        if not path.exists():
            print(f"MISSING: {path}")
            return 1
        frames[name] = pd.read_csv(path)

    # Third arm: ART as the Python analyses actually apply it — scrubbing AND the
    # >=20%-outlier scan exclusion (n 11/12/16). This is the only configuration in
    # which the full ART motion protocol is in force, so it is the fair test of
    # "does aggressive motion control change the answer?".
    incl_path = QC / "cambridge_propofol" / "analysis_inclusion.csv"
    if incl_path.exists():
        incl = pd.read_csv(incl_path, comment="#")
        keep = set(incl.loc[incl.included.astype(str).str.lower().eq("true"), "scan_id"])
        art_gated = frames["ART"][frames["ART"].scan_id.isin(keep)].copy()
        frames["ART+excl"] = art_gated
        print(f"ART+excl arm: {len(art_gated)} of {len(frames['ART'])} scans retained "
              f"by the >=20% ART exclusion\n")
    else:
        print(f"note: {incl_path} absent — skipping the ART+excl arm "
              f"(run scripts/qc/export_analysis_inclusion.py)\n")

    rows = []
    for variant, df in frames.items():
        for drug in DRUGS:
            for measure in MEASURES:
                p = paired_frame(df, drug, measure)
                res = analyse(p)
                res.update(variant=variant, contrast=f"{drug}_vs_{CONTROL}",
                           measure=measure.replace("_workspace_mean", ""))
                rows.append(res)
    res = pd.DataFrame(rows)

    out_dir = PROJECT_DIR / "results" / "qc" / "propofol_motion_sensitivity"
    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "motion_sensitivity.csv", index=False)

    pd.set_option("display.width", 220)

    print("=" * 108)
    print("1. IS MOTION CONFOUNDED WITH CONDITION?  (per-subject dFD = drug - awake, mm)")
    print("=" * 108)
    mo = (res[res.measure == "synergy"]
          .pivot(index="contrast", columns="variant", values="mean_dfd"))
    print(mo.to_string(float_format=lambda v: f"{v:+.4f}"))

    print("\n" + "=" * 108)
    print("2. MOTION COUPLING: r(dFD, dMeasure) — does motion track the effect?")
    print("=" * 108)
    mc = res.pivot_table(index=["contrast", "measure"], columns="variant",
                         values=["r_motion", "p_motion"])
    print(mc.to_string(float_format=lambda v: f"{v:+.3f}"))

    print("\n" + "=" * 108)
    print("3. EFFECT, UNADJUSTED vs dFD-ADJUSTED vs MOTION-MATCHED HALF")
    print("=" * 108)
    for variant in VARIANTS:
        print(f"\n--- {variant} ---")
        sub = res[res.variant == variant][
            ["contrast", "measure", "n", "mean_d", "t", "p",
             "adj_t", "adj_p", "n_matched", "matched_t", "matched_p"]]
        print(sub.to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    print("\n" + "=" * 108)
    print("4. SIGNIFICANCE STABILITY (does any call change under motion control?)")
    print("=" * 108)
    flips = []
    for _, r in res.iterrows():
        base = r.get("p", np.nan) < 0.05
        for label, col in (("dFD-adjusted", "adj_p"), ("motion-matched", "matched_p")):
            v = r.get(col, np.nan)
            if pd.notna(v) and (v < 0.05) != base:
                flips.append(f"  {r.variant:8s} {r.contrast:26s} {r.measure:12s} "
                             f"raw p={r['p']:.4f} -> {label} p={v:.4f}")
    print("\n".join(flips) if flips else "  none — every significance call survives both motion controls")

    print(f"\nwritten: {out_dir / 'motion_sensitivity.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
