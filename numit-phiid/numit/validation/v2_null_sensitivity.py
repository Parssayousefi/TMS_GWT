# VALIDATION 2: null-ensemble sensitivity. Do the drug conclusions depend on
# how the nulls are anchored? Four configurations:
#   prod    - scan-level empirical-A pool, Wishart df=2 (production choice)
#   vanilla - A ~ N(0,1) iid entries (the paper's default null)
#   armpool - empirical-A pool drawn from OTHER scans of the same arm
#   df5     - scan-level pool, Wishart df=5 (different noise-shape dispersion)
# Arms: basel_lpm/psilocybin (strongest effect), berlin_ketamine (extreme
# regime), monash tau=1 (highest autocorrelation; expected-null case).
# 8000 subsampled workspace pairs/scan (seeded), N=100 nulls.
# Usage: v2_null_sensitivity.py SHARD NSHARDS   (results: one npz per scan)
# TEMPORARY - tmp-only investigation.
import glob
import hashlib
import os
import re
import sys
import time

import numpy as np
import pandas as pd

BASE = r"D:\projects\TMS_GWT\numit-phiid\numit"
sys.path.insert(0, BASE)
import phiid_gauss as pg  # noqa: E402
import numit_null as nn  # noqa: E402

TS = r"D:\projects\TMS_GWT\pipeline\results\timeseries"
OUT = os.path.join(BASE, "validation", "v2_results")
os.makedirs(OUT, exist_ok=True)

N_PAIR = 8000
N_NULL = 100

ws = pd.read_csv(r"D:\projects\TMS_GWT\pipeline\results\group_analysis\hcp\rest"
                 r"\hcp_rest_workspace_summary.csv")
WS_IDX = np.where(ws.workspace_member.astype(bool).values)[0]
iu = np.triu_indices(len(WS_IDX), 1)
WS_PAIRS = np.stack([WS_IDX[iu[0]], WS_IDX[iu[1]]], axis=1)

ARMS = {
    "basel_lpm": ("psilocybin", "placebo", 1),
    "berlin_ketamine": ("ketamine", "placebo", 1),
    "monash_psilo": ("psilocybin", "baseline", 1),
}


def arm_scan_files(ds):
    drug, ctrl, _ = ARMS[ds]
    out = []
    for cond in (drug, ctrl):
        out += sorted(glob.glob(os.path.join(
            TS, f"dataset-{ds}_sub-*_ses-{cond}_task-*_run-*"
                "_timeseries_deconv.npy")))
    return out


def scan_pairs_A_obs(path, lag):
    """Subsampled ws pairs: observed atoms + per-pair VAR params."""
    seed = int(hashlib.md5(os.path.basename(path).encode()).hexdigest()[:8],
               16)
    r = np.random.default_rng(seed)
    idx = r.choice(len(WS_PAIRS), size=N_PAIR, replace=False)
    prs = WS_PAIRS[idx]
    ts = np.load(path)
    with np.errstate(invalid="ignore", divide="ignore"):
        sigma = pg.lagged_sigma_from_timeseries(ts, lag, prs)
        obs = pg.phiid_atoms(sigma)
        A, _, _ = pg.var1_from_sigma(sigma)
    return obs, A


ARM_POOL_FILE = os.path.join(OUT, "arm_pools.npz")


def build_arm_pools():
    """Pool of empirical As per arm from 6 scans each (for 'armpool')."""
    if os.path.exists(ARM_POOL_FILE):
        return dict(np.load(ARM_POOL_FILE))
    pools = {}
    r = np.random.default_rng(3)
    for ds in ARMS:
        files = arm_scan_files(ds)
        picks = r.choice(len(files), size=min(6, len(files)), replace=False)
        As = []
        for k in picks:
            _, A = scan_pairs_A_obs(files[k], ARMS[ds][2])
            good = np.all(np.isfinite(A.reshape(len(A), 4)), axis=1)
            As.append(A[good])
        pools[ds] = np.concatenate(As)[:50000]
        print(f"arm pool {ds}: {len(pools[ds])} As", flush=True)
    tmpf = os.path.join(OUT, f"tmp_pools_{os.getpid()}.npz")
    np.savez(tmpf, **pools)
    os.replace(tmpf, ARM_POOL_FILE)
    return pools


def process(path, ds, lag, arm_pools, van_pool):
    scan = os.path.basename(path).replace("_timeseries_deconv.npy", "")
    out_f = os.path.join(OUT, f"{scan}.npz")
    if os.path.exists(out_f):
        return "skip"
    t0 = time.time()
    obs, A = scan_pairs_A_obs(path, lag)
    ok = (np.isfinite(obs["tdmi"]) & (obs["tdmi"] > 1e-6)
          & np.all(np.isfinite(A.reshape(len(A), 4)), axis=1))
    seed = int(hashlib.md5((scan + "null").encode()).hexdigest()[:8], 16)
    configs = {
        "prod": (A[ok], 2),
        "vanilla": (van_pool, 2),
        "armpool": (arm_pools[ds], 2),
        "df5": (A[ok], 5),
    }
    means = {}
    for cfg, (pool, df) in configs.items():
        r = np.random.default_rng(seed)
        with np.errstate(invalid="ignore", divide="ignore"):
            nulls = nn.null_atoms_for_targets(pool, obs["tdmi"][ok], N_NULL,
                                              r, wishart_df=df)
        for m in nn.MEASURES:
            q = nn.quantile_of(obs[m][ok], nulls[m])
            means[f"{m}_{cfg}"] = float(np.nanmean(q))
            if m == "synergy":
                means[f"sat_{cfg}"] = float(np.nanmean((q == 0) | (q == 1)))
    tmpf = os.path.join(OUT, f"tmp_{os.getpid()}.npz")
    np.savez(tmpf, **means)
    os.replace(tmpf, out_f)
    return f"{time.time() - t0:.0f}s"


if __name__ == "__main__":
    shard = int(sys.argv[1])
    nshards = int(sys.argv[2])
    arm_pools = build_arm_pools()
    van_pool = np.random.default_rng(9).normal(0, 1, (50000, 2, 2))
    jobs = []
    for ds in ARMS:
        for f in arm_scan_files(ds):
            jobs.append((f, ds, ARMS[ds][2]))
    jobs = sorted(jobs)
    mine = [j for k, j in enumerate(jobs) if k % nshards == shard]
    print(f"shard {shard}/{nshards}: {len(mine)}/{len(jobs)} jobs", flush=True)
    for k, (f, ds, lag) in enumerate(mine):
        try:
            r = process(f, ds, lag, arm_pools, van_pool)
        except Exception as e:
            r = f"ERROR {type(e).__name__}: {e}"
        print(f"[{shard}] {k + 1}/{len(mine)} {os.path.basename(f)[:55]}: {r}",
              flush=True)
    print(f"shard {shard} DONE", flush=True)
