# VALIDATION 3: finite-sample calibration (parametric bootstrap).
# Production computes OBSERVED atoms from covariances estimated on T~200-500
# volumes, but compares them against ANALYTIC (noise-free) null atoms.
# Question: does that asymmetry distort the quantiles?
# Test: draw systems from the null ensemble, SIMULATE series of length T,
# estimate atoms exactly as production does (z-scored past/future blocks,
# sample covariance), quantile against analytic nulls matched to the
# ESTIMATED TDMI (mirroring production). If ~Uniform(0,1): calibrated.
# T = 200 (berlin-like), 350 (basel-like), 500 (monash-like).
# TEMPORARY - tmp-only investigation.
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = r"D:\projects\TMS_GWT\numit-phiid\numit"
sys.path.insert(0, BASE)
import phiid_gauss as pg  # noqa: E402
import numit_null as nn  # noqa: E402

rng = np.random.default_rng(23)
N_PROBE = 1200
N_NULL = 100

# empirical pool (same scan as production-style validation)
ws = pd.read_csv(r"D:\projects\TMS_GWT\pipeline\results\group_analysis\hcp\rest"
                 r"\hcp_rest_workspace_summary.csv")
ws_idx = np.where(ws.workspace_member.astype(bool).values)[0]
ts = np.load(r"D:\projects\TMS_GWT\pipeline\results\timeseries"
             r"\dataset-basel_lpm_sub-001_ses-placebo_task-rest_run-1"
             r"_timeseries_deconv.npy")
iu = np.triu_indices(len(ws_idx), 1)
pairs = np.stack([ws_idx[iu[0]], ws_idx[iu[1]]], axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    A_scan, _, _ = pg.var1_from_sigma(
        pg.lagged_sigma_from_timeseries(ts, 1, pairs))
POOL = A_scan[np.all(np.isfinite(A_scan.reshape(len(pairs), 4)), axis=1)]

# ---- draw probe systems from the null ensemble (production procedure) ----
tgt = rng.uniform(0.2, 1.8, N_PROBE)
with np.errstate(invalid="ignore", divide="ignore"):
    probe_atoms = nn.null_atoms_for_targets(POOL, tgt, 1, rng)
# regenerate the actual probe systems (A, V, f) by reproducing the draw is
# awkward; instead draw systems directly here the same way:
B = N_PROBE
A = POOL[rng.integers(0, len(POOL), size=B)]
d, P, Pinv, detP = nn._eig2(A)
rho = np.maximum(np.abs(d).max(axis=-1), 1e-12)
c11 = np.sqrt(rng.chisquare(2, B))
c22 = np.sqrt(rng.chisquare(1, B))
n21 = rng.normal(0, 1, B)
V = np.empty((B, 2, 2))
V[:, 0, 0] = c11 ** 2
V[:, 0, 1] = V[:, 1, 0] = c11 * n21
V[:, 1, 1] = n21 ** 2 + c22 ** 2

# rescale each A to hit its target TDMI via bisection on true system
Ahat = A / rho[:, None, None]
lo = np.full(B, 1e-4)
hi = np.full(B, 1 - 1e-6)
for _ in range(40):
    mid = 0.5 * (lo + hi)
    s4, S0 = pg.sigma4_from_var1(Ahat * mid[:, None, None], V, lag=1)
    td = pg.tdmi_from_sigma4(s4)
    up = td < tgt
    lo = np.where(up, mid, lo)
    hi = np.where(up, hi, mid)
f = 0.5 * (lo + hi)
Atrue = Ahat * f[:, None, None]

print(f"probe systems: {B}, target TDMI in [0.2, 1.8]")

results = []
for T in (200, 350, 500):
    t0 = time.time()
    est_sigma = np.empty((B, 4, 4))
    for k in range(B):
        # simulate T+100 steps, drop burn-in
        e = rng.multivariate_normal([0, 0], V[k], size=T + 100)
        x = np.zeros((T + 100, 2))
        for t in range(1, T + 100):
            x[t] = Atrue[k] @ x[t - 1] + e[t]
        x = x[100:]
        est_sigma[k] = pg.lagged_sigma_from_timeseries(x, 1,
                                                       np.array([[0, 1]]))[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        obs = pg.phiid_atoms(est_sigma)
        est_tdmi = obs["tdmi"]
        ok = np.isfinite(est_tdmi) & (est_tdmi > 1e-6)
        nulls = nn.null_atoms_for_targets(POOL, est_tdmi[ok], N_NULL,
                                          np.random.default_rng(T))
    row = {"T": T, "n": int(ok.sum())}
    for m in ("synergy", "redundancy", "phi_r"):
        q = nn.quantile_of(obs[m][ok], nulls[m])
        qo = q[np.isfinite(q)]
        ks = sps.kstest(qo, "uniform")
        row[f"{m}_mean"] = round(float(qo.mean()), 3)
        row[f"{m}_KS"] = round(float(ks.statistic), 3)
        row[f"{m}_KSp"] = float(f"{ks.pvalue:.3f}")
        row[f"{m}_tail"] = round(float(np.mean((qo < 0.05) | (qo > 0.95))), 3)
    results.append(row)
    print(f"T={T} done ({time.time() - t0:.0f}s)", flush=True)

df = pd.DataFrame(results)
pd.set_option("display.width", 200)
print()
print(df.to_string(index=False))
print("\nInterpretation: mean ~0.50 = unbiased; KS p > .01 = uniform;")
print("tail fraction ~0.10 expected for perfectly calibrated quantiles.")
print("Mild tail inflation = estimation noise widens observed atoms symmetrically;")
print("within-subject drug CONTRASTS difference this out (both conditions share T).")
