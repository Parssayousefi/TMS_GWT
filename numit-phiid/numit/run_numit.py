# B4 driver: NuMIT-PhiID normalization over all paired arms, full-brain 454,
# N=100 nulls. Usage:  run_numit.py SHARD NSHARDS
# Jobs: all 7 paired arms at tau=1  +  monash additionally at tau=3 (the lag
# where the raw synergy effect emerges - 2026-08-25 finding).
# Per (scan, lag): saves
#   tmp/numit/results/{scan}_lag{L}.npz            (ws/brain means: raw, NMI,
#       NuMIT for 5 measures + tdmi + QC)
#   tmp/numit/results/{scan}_lag{L}_q_{sts|rtr|phi_r}.npy  (454x454 float32
#       quantile matrices)
# Checkpointed (skip if npz exists), atomic writes (tmp->rename), shard-safe.
# TEMPORARY - tmp-only investigation.
import glob
import hashlib
import os
import re
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phiid_gauss as pg  # noqa: E402
import numit_null as nn  # noqa: E402

PIPE = r"D:\projects\TMS_GWT\pipeline"
TS = os.path.join(PIPE, "results", "timeseries")
OUT = r"D:\projects\TMS_GWT\numit-phiid\numit\results"
os.makedirs(OUT, exist_ok=True)

N_NULL = 100
CHUNK = 8192
N = 454
MEASURES = nn.MEASURES  # synergy, redundancy, synergy_total, phi, phi_r

ws = pd.read_csv(os.path.join(PIPE, "results", "group_analysis", "hcp", "rest",
                              "hcp_rest_workspace_summary.csv"))
WS_MASK = ws.workspace_member.astype(bool).values

ARMS = [
    ("monash_psilo", "psilocybin", "baseline", "ses"),
    ("basel_lpm", "psilocybin", "placebo", "ses"),
    ("basel_lpm", "lsd", "placebo", "ses"),
    ("basel_lpm", "mescaline", "placebo", "ses"),
    ("basel_lam", "lsd", "placebo", "ses"),
    ("basel_lam", "mdma", "placebo", "ses"),
    ("basel_lam", "amph", "placebo", "ses"),
    ("london_lsd", "LSD", "PCB", "ses"),
    ("berlin_ketamine", "ketamine", "placebo", "ses"),
    ("cambridge_propofol_uniform", "restdeep", "restawake", "task"),
]


def arm_files(ds, cond, style):
    if style == "task":
        pat = f"dataset-{ds}_sub-*_task-{cond}_run-*_timeseries_deconv.npy"
    else:
        pat = f"dataset-{ds}_sub-*_ses-{cond}_task-*_run-*_timeseries_deconv.npy"
    return {re.search(r"sub-([^_]+)", os.path.basename(f)).group(1): f
            for f in sorted(glob.glob(os.path.join(TS, pat)))}


def build_jobs():
    jobs = []
    seen = set()
    for ds, drug, ctrl, style in ARMS:
        fd, fc = arm_files(ds, drug, style), arm_files(ds, ctrl, style)
        subs = sorted(set(fd) & set(fc))
        lags = (1, 3) if ds == "monash_psilo" else (1,)
        for lag in lags:
            for s in subs:
                for f in (fd[s], fc[s]):
                    key = (f, lag)
                    if key not in seen:
                        seen.add(key)
                        jobs.append(key)
    return sorted(jobs)


def process(path, lag):
    scan = os.path.basename(path).replace("_timeseries_deconv.npy", "")
    out_npz = os.path.join(OUT, f"{scan}_lag{lag}.npz")
    if os.path.exists(out_npz):
        return "skip"
    t0 = time.time()
    seed = int(hashlib.md5(f"{scan}_lag{lag}".encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)

    ts = np.load(path)
    valid = np.where(ts.std(0) > 0)[0]
    vi, vj = np.triu_indices(len(valid), 1)
    pairs = np.stack([valid[vi], valid[vj]], axis=1)
    P = len(pairs)

    with np.errstate(invalid="ignore", divide="ignore"):
        sigma = pg.lagged_sigma_from_timeseries(ts, lag, pairs)
        obs = pg.phiid_atoms(sigma)
        A_all, V_all, _ = pg.var1_from_sigma(sigma)

    tdmi = obs["tdmi"]
    ok = np.isfinite(tdmi) & (tdmi > 1e-6)
    ok &= np.all(np.isfinite(A_all.reshape(P, 4)), axis=1)
    pool = A_all[ok]

    quant = {m: np.full(P, np.nan, dtype=np.float32) for m in MEASURES}
    tdmi_err_max = 0.0
    idx_ok = np.where(ok)[0]
    for c0 in range(0, len(idx_ok), CHUNK):
        ci = idx_ok[c0:c0 + CHUNK]
        with np.errstate(invalid="ignore", divide="ignore"):
            nulls = nn.null_atoms_for_targets(A_all[ci], tdmi[ci], N_NULL, rng)
        tdmi_err_max = max(tdmi_err_max, float(np.nanmax(nulls["_tdmi_err"])))
        for m in MEASURES:
            quant[m][ci] = nn.quantile_of(obs[m][ci], nulls[m])

    # assemble outputs
    means = {}
    brain_mask_pairs = np.ones(P, dtype=bool)
    ws_pairs = WS_MASK[pairs[:, 0]] & WS_MASK[pairs[:, 1]]
    for m in MEASURES:
        raw = obs[m]
        nmi = raw / tdmi
        for tag, arr in (("raw", raw), ("nmi", nmi), ("numit", quant[m])):
            for reg, mask in (("brain", brain_mask_pairs), ("ws", ws_pairs)):
                means[f"{m}_{tag}_{reg}"] = float(np.nanmean(arr[mask]))
    for reg, mask in (("brain", brain_mask_pairs), ("ws", ws_pairs)):
        means[f"tdmi_raw_{reg}"] = float(np.nanmean(tdmi[mask]))
    qs = quant["synergy"]
    means["qc_tdmi_err_max"] = tdmi_err_max
    means["qc_frac_pairs_ok"] = float(ok.mean())
    means["qc_syn_q_sat"] = float(np.nanmean((qs == 0) | (qs == 1)))

    # quantile matrices for sts, rtr, phi_r
    for m, tag in (("synergy", "sts"), ("redundancy", "rtr"), ("phi_r", "phi_r")):
        Q = np.full((N, N), np.nan, dtype=np.float32)
        Q[pairs[:, 0], pairs[:, 1]] = quant[m]
        Q[pairs[:, 1], pairs[:, 0]] = quant[m]
        tmpf = os.path.join(OUT, f".tmp_{os.getpid()}_{tag}.npy")
        np.save(tmpf, Q)
        os.replace(tmpf, os.path.join(OUT, f"{scan}_lag{lag}_q_{tag}.npy"))

    tmpf = os.path.join(OUT, f".tmp_{os.getpid()}.npz")
    np.savez(tmpf, **means)
    os.replace(tmpf, out_npz)
    return f"{time.time() - t0:.0f}s"


if __name__ == "__main__":
    shard = int(sys.argv[1])
    nshards = int(sys.argv[2])
    jobs = build_jobs()
    mine = [j for k, j in enumerate(jobs) if k % nshards == shard]
    if len(sys.argv) > 3 and sys.argv[3] == "rev":
        mine = mine[::-1]  # sweeper: same partition, reverse order; skips
        # already-done jobs, meets the forward worker in the middle
    print(f"shard {shard}/{nshards}: {len(mine)} of {len(jobs)} jobs",
          flush=True)
    done = 0
    for path, lag in mine:
        try:
            r = process(path, lag)
        except Exception as e:  # keep the shard alive
            r = f"ERROR {type(e).__name__}: {e}"
        done += 1
        print(f"[{shard}] {done}/{len(mine)} "
              f"{os.path.basename(path)[:60]} lag{lag}: {r}", flush=True)
    print(f"shard {shard} DONE", flush=True)
