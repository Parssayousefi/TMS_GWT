# B5: analyze NuMIT results - per-arm paired drug-vs-control tests on
# raw / NMI / NuMIT values (workspace + brainwide means), the fMRI analog of
# Liardi et al. Fig 4. Plus QC summary. Runs on whatever results exist
# (checkpointed driver), so it can be run mid-stream.
# TEMPORARY - tmp-only investigation.
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

OUT = r"D:\projects\TMS_GWT\numit-phiid\numit\results"
MEASURES = ["synergy", "redundancy", "synergy_total", "phi", "phi_r"]

ARMS = [
    ("monash_psilo", "psilocybin", "baseline", 1),
    ("monash_psilo", "psilocybin", "baseline", 3),
    ("basel_lpm", "psilocybin", "placebo", 1),
    ("basel_lpm", "lsd", "placebo", 1),
    ("basel_lpm", "mescaline", "placebo", 1),
    ("basel_lam", "lsd", "placebo", 1),
    ("basel_lam", "mdma", "placebo", 1),
    ("basel_lam", "amph", "placebo", 1),
    ("london_lsd", "LSD", "PCB", 1),
    ("berlin_ketamine", "ketamine", "placebo", 1),
    ("cambridge_propofol_uniform", "restdeep", "restawake", 1),
]


def cond_files(ds, cond, lag):
    """scan npz files for one condition (session- or task-named)."""
    out = {}
    for f in glob.glob(os.path.join(OUT, f"dataset-{ds}_sub-*_lag{lag}.npz")):
        base = os.path.basename(f)
        if f"_ses-{cond}_" not in base and f"_task-{cond}_" not in base:
            continue
        out[re.search(r"sub-([^_]+)", base).group(1)] = f
    return out


rows = []
qc_rows = []
for ds, drug, ctrl, lag in ARMS:
    fd = cond_files(ds, drug, lag)
    fc = cond_files(ds, ctrl, lag)
    subs = sorted(set(fd) & set(fc))
    if len(subs) < 6:
        rows.append({"arm": f"{ds}/{drug}", "lag": lag, "n": len(subs),
                     "note": "insufficient pairs yet"})
        continue
    dd = {s: dict(np.load(fd[s])) for s in subs}
    cc = {s: dict(np.load(fc[s])) for s in subs}
    for s in subs:
        for tag, blob in (("drug", dd[s]), ("ctrl", cc[s])):
            qc_rows.append({"arm": f"{ds}/{drug}", "lag": lag,
                            "tdmi_err": float(blob["qc_tdmi_err_max"]),
                            "sat": float(blob["qc_syn_q_sat"])})
    d_tdmi = np.array([dd[s]["tdmi_raw_ws"] - cc[s]["tdmi_raw_ws"]
                       for s in subs], dtype=float)
    for m in MEASURES:
        row = {"arm": f"{ds}/{drug}", "lag": lag, "n": len(subs),
               "measure": m}
        # ANCOVA sensitivity check: Delta-measure ~ a + b*Delta-TDMI;
        # test the intercept (effect at zero TDMI change). Linear - kept as
        # a robustness column, NOT a substitute for the null models.
        dy = np.array([dd[s][f"{m}_raw_ws"] - cc[s][f"{m}_raw_ws"]
                       for s in subs], dtype=float)
        X = np.column_stack([np.ones(len(subs)), d_tdmi])
        beta, res_ss = np.linalg.lstsq(X, dy, rcond=None)[:2]
        dof = len(subs) - 2
        sig2 = float(res_ss[0]) / dof if len(res_ss) else np.nan
        cov = sig2 * np.linalg.inv(X.T @ X)
        t_int = beta[0] / np.sqrt(cov[0, 0])
        p_int = 2 * stats.t.sf(abs(t_int), dof)
        row["ancova_t"] = round(float(t_int), 2)
        row["ancova_p"] = float(f"{p_int:.2e}")
        row["ancova_slope"] = round(float(beta[1]), 3)
        for tag in ("raw", "nmi", "numit"):
            for reg in ("ws", "brain"):
                key = f"{m}_{tag}_{reg}"
                dv = np.array([dd[s][key] for s in subs], dtype=float)
                cv = np.array([cc[s][key] for s in subs], dtype=float)
                t, p = stats.ttest_rel(dv, cv)
                if tag == "numit":
                    eff = float((dv - cv).mean())        # quantile shift
                else:
                    eff = 100 * float((dv - cv).mean() / cv.mean())
                if reg == "ws":
                    row[f"{tag}_eff"] = round(eff, 3 if tag == "numit" else 1)
                    row[f"{tag}_t"] = round(float(t), 2)
                    row[f"{tag}_p"] = float(f"{p:.2e}")
                else:
                    row[f"{tag}_t_brain"] = round(float(t), 2)
        rows.append(row)
    # TDMI raw for reference
    dv = np.array([dd[s]["tdmi_raw_ws"] for s in subs], dtype=float)
    cv = np.array([cc[s]["tdmi_raw_ws"] for s in subs], dtype=float)
    t, p = stats.ttest_rel(dv, cv)
    rows.append({"arm": f"{ds}/{drug}", "lag": lag, "n": len(subs),
                 "measure": "tdmi", "raw_eff":
                 round(100 * float((dv - cv).mean() / cv.mean()), 1),
                 "raw_t": round(float(t), 2), "raw_p": float(f"{p:.2e}")})

df = pd.DataFrame(rows)
df.to_csv(r"D:\projects\TMS_GWT\numit-phiid\numit\numit_arm_table.csv", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
cols = ["arm", "lag", "n", "measure", "raw_eff", "raw_t", "raw_p",
        "nmi_eff", "nmi_t", "nmi_p", "ancova_t", "ancova_p",
        "numit_eff", "numit_t", "numit_p"]
cols = [c for c in cols if c in df.columns]
print(df[cols].to_string(index=False))

if qc_rows:
    q = pd.DataFrame(qc_rows)
    print("\n=== QC (per arm): max TDMI match err / mean saturated-quantile "
          "fraction (synergy) ===")
    print(q.groupby(["arm", "lag"]).agg(
        tdmi_err=("tdmi_err", "max"), sat=("sat", "mean"),
        n_scans=("sat", "size")).round(4).to_string())
