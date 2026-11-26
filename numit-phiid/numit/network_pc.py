# Point-1 analysis: network-resolved synergy-PC drug shifts among workspace
# ROIs - does the FPN/network-specific gateway pattern survive normalization?
# Computes per-network mean Delta(synergy-PC) on RAW (production) matrices vs
# NuMIT quantile matrices for the pooled serotonergic arms + per-arm FPN test.
# TEMPORARY - tmp-only investigation.
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

PIPE = r"D:\projects\TMS_GWT\pipeline"
RAWD = os.path.join(PIPE, "results", "phiid_matrices")
RES = r"D:\projects\TMS_GWT\numit-phiid\numit\results"
roi = pd.read_csv(os.path.join(PIPE, "atlases", "roi_to_network_mapping.csv"))
roi = roi[roi.ROI > 0]
modules = roi.Network_ID.values
net_names = (roi.drop_duplicates("Network_ID")
             .sort_values("Network_ID").Network_Name.tolist())
ws = pd.read_csv(os.path.join(PIPE, "results", "group_analysis", "hcp", "rest",
                              "hcp_rest_workspace_summary.csv"))
ws_mask = ws.workspace_member.astype(bool).values
N = 454

ARMS = [
    ("basel_lpm", "psilocybin", "placebo"),
    ("basel_lpm", "lsd", "placebo"),
    ("basel_lpm", "mescaline", "placebo"),
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


def files(ds, cond, kind):
    """kind: 'raw' -> production synergy matrices; 'q' -> NuMIT quantiles."""
    if kind == "raw":
        pat = os.path.join(RAWD, f"dataset-{ds}_sub-*_ses-{cond}_task-*"
                                 "_run-*_synergy.npy")
    else:
        pat = os.path.join(RES, f"dataset-{ds}_sub-*_lag1_q_sts.npy")
    out = {}
    for f in glob.glob(pat):
        b = os.path.basename(f)
        if kind == "q" and f"_ses-{cond}_" not in b:
            continue
        if kind == "raw" and "_synergy_total" in b:
            continue
        out[re.search(r"sub-([^_]+)", b).group(1)] = f
    return out


# pooled per-subject Delta-PC profiles (one per subject per arm)
profiles = {"raw": [], "q": []}
for kind in ("raw", "q"):
    for ds, drug, ctrl in ARMS:
        fd, fc = files(ds, drug, kind), files(ds, ctrl, kind)
        subs = sorted(set(fd) & set(fc))
        for s in subs:
            dpc = (compute_pc(np.load(fd[s]).astype(np.float64))
                   - compute_pc(np.load(fc[s]).astype(np.float64)))
            profiles[kind].append(dpc)
    profiles[kind] = np.array(profiles[kind])
    print(f"{kind}: {profiles[kind].shape[0]} subject-arm profiles", flush=True)

rows = []
for nid, name in enumerate(net_names):
    sel = ws_mask & (modules == nid)
    if sel.sum() == 0:
        continue
    row = {"network": name, "n_ws_rois": int(sel.sum())}
    for kind in ("raw", "q"):
        v = np.nanmean(profiles[kind][:, sel], axis=1)
        t, p = stats.ttest_1samp(v[np.isfinite(v)], 0.0)
        row[f"{kind}_d"] = round(float(np.nanmean(v)), 4)
        row[f"{kind}_t"] = round(float(t), 2)
        row[f"{kind}_p"] = float(f"{p:.4f}")
    rows.append(row)
df = pd.DataFrame(rows)
df.to_csv(r"D:\projects\TMS_GWT\numit-phiid\numit\network_pc.csv", index=False)
pd.set_option("display.width", 160)
print("\nPooled Basel-LPM psychedelic arms (psilo+lsd+mescaline, 92 pairs):")
print("per-network mean Delta synergy-PC among WORKSPACE ROIs")
print(df.to_string(index=False))
# raw-vs-normalized network-profile agreement
a = df.raw_d.values
b = df.q_d.values
r = np.corrcoef(a, b)[0, 1]
print(f"\nnetwork-profile correlation raw vs normalized: r={r:.2f}")
