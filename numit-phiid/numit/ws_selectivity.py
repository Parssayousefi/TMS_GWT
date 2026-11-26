# Point-4 analysis: is the NuMIT-normalized synergy increase WORKSPACE-
# SELECTIVE? Per scan, mean normalized quantile over three pair compartments:
#   ws-ws (both ROIs in the HCP workspace), ws-out (one in), out-out (none).
# Paired drug-vs-control t per compartment per arm + direct test of
# (ws-ws minus out-out) contrast-of-contrasts.
# TEMPORARY - tmp-only investigation.
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

RES = r"D:\projects\TMS_GWT\numit-phiid\numit\results"
ws = pd.read_csv(r"D:\projects\TMS_GWT\pipeline\results\group_analysis\hcp\rest"
                 r"\hcp_rest_workspace_summary.csv")
WS = ws.workspace_member.astype(bool).values
N = 454
iu = np.triu_indices(N, 1)
in_i = WS[iu[0]]
in_j = WS[iu[1]]
COMP = {"ws-ws": in_i & in_j, "ws-out": in_i ^ in_j,
        "out-out": ~in_i & ~in_j}

ARMS = [
    ("basel_lpm", "psilocybin", "placebo", 1),
    ("basel_lpm", "lsd", "placebo", 1),
    ("basel_lpm", "mescaline", "placebo", 1),
    ("basel_lam", "lsd", "placebo", 1),
    ("london_lsd", "LSD", "PCB", 1),
    ("monash_psilo", "psilocybin", "baseline", 3),
    ("berlin_ketamine", "ketamine", "placebo", 1),
]


def cond_files(ds, cond, lag, tag):
    out = {}
    for f in glob.glob(os.path.join(
            RES, f"dataset-{ds}_sub-*_lag{lag}_q_{tag}.npy")):
        b = os.path.basename(f)
        if f"_ses-{cond}_" not in b and f"_task-{cond}_" not in b:
            continue
        out[re.search(r"sub-([^_]+)", b).group(1)] = f
    return out


def comp_means(f):
    q = np.load(f).astype(np.float64)[iu]
    return {c: float(np.nanmean(q[m])) for c, m in COMP.items()}


rows = []
for ds, drug, ctrl, lag in ARMS:
    for tag in ("sts", "rtr"):
        fd = cond_files(ds, drug, lag, tag)
        fc = cond_files(ds, ctrl, lag, tag)
        subs = sorted(set(fd) & set(fc))
        if len(subs) < 6:
            continue
        md = [comp_means(fd[s]) for s in subs]
        mc = [comp_means(fc[s]) for s in subs]
        row = {"arm": f"{ds}/{drug}", "lag": lag, "meas": tag, "n": len(subs)}
        for c in COMP:
            d = np.array([x[c] for x in md]) - np.array([x[c] for x in mc])
            t, p = stats.ttest_rel(np.array([x[c] for x in md]),
                                   np.array([x[c] for x in mc]))
            row[f"{c}_d"] = round(float(d.mean()), 4)
            row[f"{c}_t"] = round(float(t), 2)
            row[f"{c}_p"] = float(f"{p:.3f}")
        sel = (np.array([x["ws-ws"] for x in md])
               - np.array([x["ws-ws"] for x in mc])
               - (np.array([x["out-out"] for x in md])
                  - np.array([x["out-out"] for x in mc])))
        t, p = stats.ttest_1samp(sel, 0.0)
        row["wsSEL_d"] = round(float(sel.mean()), 4)
        row["wsSEL_t"] = round(float(t), 2)
        row["wsSEL_p"] = float(f"{p:.3f}")
        rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(r"D:\projects\TMS_GWT\numit-phiid\numit\ws_selectivity.csv", index=False)
pd.set_option("display.width", 240)
print(df.to_string(index=False))
print("\nwsSEL = (drug-ctrl in ws-ws) minus (drug-ctrl in out-out):")
print("positive & significant = normalized effect concentrated in workspace")
