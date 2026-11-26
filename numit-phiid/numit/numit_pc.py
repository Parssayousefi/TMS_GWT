# Topology under normalization: participation-coefficient analysis on the
# NuMIT quantile matrices (q_sts / q_rtr) - does Basel's synergy-PC
# redistribution survive normalization, does Monash gain one (lag 3)?
# Mirrors the notebook PC pipeline (8-network partition, NaN-aware, HCP
# workspace, within-workspace stats). Runs on whatever results exist.
# TEMPORARY - tmp-only investigation.
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats as sps

PIPE = r"D:\projects\TMS_GWT\pipeline"
RES = r"D:\projects\TMS_GWT\numit-phiid\numit\results"
rng = np.random.default_rng(42)

roi = pd.read_csv(os.path.join(PIPE, "atlases", "roi_to_network_mapping.csv"))
roi = roi[roi.ROI > 0]
modules = roi.Network_ID.values
ws = pd.read_csv(os.path.join(PIPE, "results", "group_analysis", "hcp", "rest",
                              "hcp_rest_workspace_summary.csv"))
ws_idx = np.where(ws.workspace_member.astype(bool).values)[0]
N = 454

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


def compute_pc(mat):
    m = mat.copy()
    np.fill_diagonal(m, 0.0)
    all_nf = ~np.isfinite(m).any(axis=1)
    k = np.nansum(m, axis=1)
    ksafe = np.where(k == 0, np.inf, k)
    ss = np.zeros(N)
    for s in np.unique(modules):
        ss += (np.nansum(m[:, modules == s], axis=1) / ksafe) ** 2
    pc = 1 - ss
    pc[np.isinf(ksafe)] = np.nan
    pc[all_nf] = np.nan
    return pc


def split_half(prof, n_splits=300):
    n = prof.shape[0]
    if n < 6:
        return np.nan
    if n % 2:
        prof = prof[:-1]
        n -= 1
    rs = []
    for _ in range(n_splits):
        idx = rng.permutation(n)
        a = np.nanmean(prof[idx[: n // 2]], axis=0)
        b = np.nanmean(prof[idx[n // 2:]], axis=0)
        ok = np.isfinite(a) & np.isfinite(b)
        rs.append(np.corrcoef(a[ok], b[ok])[0, 1])
    return float(np.mean(rs))


def cond_files(ds, cond, lag, tag):
    out = {}
    for f in glob.glob(os.path.join(
            RES, f"dataset-{ds}_sub-*_lag{lag}_q_{tag}.npy")):
        base = os.path.basename(f)
        if f"_ses-{cond}_" not in base and f"_task-{cond}_" not in base:
            continue
        out[re.search(r"sub-([^_]+)", base).group(1)] = f
    return out


rows = []
for ds, drug, ctrl, lag in ARMS:
    for tag, name in (("sts", "synergy_pc"), ("rtr", "redundancy_pc")):
        fd = cond_files(ds, drug, lag, tag)
        fc = cond_files(ds, ctrl, lag, tag)
        subs = sorted(set(fd) & set(fc))
        if len(subs) < 6:
            rows.append({"arm": f"{ds}/{drug}", "lag": lag, "measure": name,
                         "n": len(subs), "note": "insufficient"})
            continue
        prof = []
        for s in subs:
            prof.append((compute_pc(np.load(fd[s]).astype(np.float64))
                         - compute_pc(np.load(fc[s]).astype(np.float64)))[ws_idx])
        prof = np.array(prof)
        wsm = np.nanmean(prof, axis=1)
        t, p = sps.ttest_1samp(wsm, 0.0)
        rows.append({"arm": f"{ds}/{drug}", "lag": lag, "measure": name,
                     "n": len(subs), "wsmean_t": round(float(t), 2),
                     "wsmean_p": float(f"{p:.2e}"),
                     "SH": round(split_half(prof), 3)})

df = pd.DataFrame(rows)
df.to_csv(r"D:\projects\TMS_GWT\numit-phiid\numit\numit_pc_table.csv", index=False)
pd.set_option("display.width", 160)
print(df.to_string(index=False))
