# VALIDATION 1: synthetic ground-truth recovery for PhiID-NuMIT
# (the PhiID analog of Liardi et al. Fig 1 + Fig 3).
#
# Probes (all analytic - no simulation noise):
#  1a NOISE-INVARIANCE: fix a system's dynamical STRUCTURE (shape of A, V),
#     sweep its signal strength (spectral radius f, the VAR analog of the
#     paper's noise gain g). Raw atoms must vary strongly (they track TDMI);
#     NuMIT quantiles should stay ~constant. Also: exact invariance of
#     everything under pure noise-covariance scaling V -> gamma*V (sanity).
#  1b STRUCTURE-TRACKING:
#     (rho-sweep) two self-coupled regions sharing correlated noise:
#         A = diag(0.55, 0.45), V = [[1, rho], [rho, 1]] - as rho rises the
#         pair becomes more redundant; NuMIT-rtr should rise monotonically.
#     (c-sweep) asymmetric cross-coupling A = [[0.5, c], [0.8c, 0.4]], V=I -
#         as c rises the joint dynamics carry more genuinely shared/joint
#         temporal structure; NuMIT synergy measures should rise.
#  Nulls: PRODUCTION procedure (empirical-A pool from a real basel_lpm
#  placebo scan, Wishart df=2, N=400); repeated with VANILLA N(0,1) pools
#  as a preview of Validation 2.
#
# TEMPORARY - tmp-only investigation.
import os
import sys

import numpy as np
import pandas as pd

BASE = r"D:\projects\TMS_GWT\numit-phiid\numit"
sys.path.insert(0, BASE)
import phiid_gauss as pg  # noqa: E402
import numit_null as nn  # noqa: E402

OUT = os.path.join(BASE, "validation")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(11)
N_NULL = 400

# ---- production-style pool from a real mid-regime scan ----
TS = r"D:\projects\TMS_GWT\pipeline\results\timeseries"
ws = pd.read_csv(r"D:\projects\TMS_GWT\pipeline\results\group_analysis\hcp\rest"
                 r"\hcp_rest_workspace_summary.csv")
ws_idx = np.where(ws.workspace_member.astype(bool).values)[0]
ts = np.load(os.path.join(
    TS, "dataset-basel_lpm_sub-001_ses-placebo_task-rest_run-1"
        "_timeseries_deconv.npy"))
iu = np.triu_indices(len(ws_idx), 1)
pairs = np.stack([ws_idx[iu[0]], ws_idx[iu[1]]], axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    sig_scan = pg.lagged_sigma_from_timeseries(ts, 1, pairs)
    A_scan, _, _ = pg.var1_from_sigma(sig_scan)
ok = np.all(np.isfinite(A_scan.reshape(len(pairs), 4)), axis=1)
POOL_PROD = A_scan[ok]
POOL_VAN = rng.normal(0, 1, (20000, 2, 2))
print(f"production pool: {len(POOL_PROD)} empirical As (basel_lpm placebo)")


def atoms_of(A, V):
    s4, _ = pg.sigma4_from_var1(A[None], V[None], lag=1)
    return pg.phiid_atoms(pg.normalize_sigma(s4))


def numit_of(A, V, pool, n=N_NULL, seed=0):
    ob = atoms_of(A, V)
    r = np.random.default_rng(seed)
    with np.errstate(invalid="ignore", divide="ignore"):
        nl = nn.null_atoms_for_targets(pool, ob["tdmi"], n, r)
    return ob, {m: float(nn.quantile_of(ob[m], nl[m])[0])
                for m in nn.MEASURES}


def rescale_radius(A, f):
    ev = np.max(np.abs(np.linalg.eigvals(A)))
    return A * (f / ev)


print("\n" + "=" * 74)
print("1a  NOISE-INVARIANCE: same structure, different signal strength")
print("    (raw atoms should vary; NuMIT quantiles should stay ~flat)")
print("=" * 74)
structures = {
    "designed asym": (np.array([[0.5, 0.3], [0.1, 0.4]]),
                      np.array([[1.0, 0.3], [0.3, 1.0]])),
    "empirical #1": (POOL_PROD[100], np.eye(2)),
    "empirical #2": (POOL_PROD[5000], np.eye(2)),
}
rows = []
for name, (A0, V0) in structures.items():
    for f in (0.35, 0.5, 0.65, 0.8, 0.9):
        A = rescale_radius(A0, f)
        ob, q = numit_of(A, V0, POOL_PROD, seed=int(f * 100))
        rows.append({"structure": name, "f": f,
                     "tdmi": round(float(ob["tdmi"][0]), 3),
                     "raw_sts": round(float(ob["synergy"][0]), 4),
                     "q_sts": round(q["synergy"], 3),
                     "raw_rtr": round(float(ob["redundancy"][0]), 4),
                     "q_rtr": round(q["redundancy"], 3)})
df = pd.DataFrame(rows)
print(df.to_string(index=False))
sums = []
for name in structures:
    d = df[df.structure == name]
    sums.append({
        "structure": name,
        "raw_sts spread (max/min)": round(d.raw_sts.max() / max(d.raw_sts.min(), 1e-9), 1),
        "q_sts SD": round(d.q_sts.std(), 3),
        "q_rtr SD": round(d.q_rtr.std(), 3)})
print("\nsummary (flat quantiles = SD near 0 while raw spreads x-fold):")
print(pd.DataFrame(sums).to_string(index=False))

# exact noise-scale invariance sanity
A0, V0 = structures["designed asym"]
a1 = atoms_of(A0, V0)
a2 = atoms_of(A0, 7.3 * V0)
d = max(abs(float(a1[m][0] - a2[m][0])) for m in nn.MEASURES + ["tdmi"])
print(f"\nexact V->7.3V invariance: max atom diff = {d:.2e} "
      f"({'OK' if d < 1e-12 else 'FAIL'})")

print("\n" + "=" * 74)
print("1b-rho  STRUCTURE-TRACKING: noise correlation -> redundancy")
print("=" * 74)
rows = []
A0 = np.diag([0.55, 0.45])
for rho in (0.0, 0.15, 0.3, 0.45, 0.6, 0.75):
    V = np.array([[1.0, rho], [rho, 1.0]])
    ob, q = numit_of(A0, V, POOL_PROD, seed=int(rho * 100) + 7)
    _, qv = numit_of(A0, V, POOL_VAN, seed=int(rho * 100) + 7)
    rows.append({"rho": rho, "tdmi": round(float(ob["tdmi"][0]), 3),
                 "raw_rtr": round(float(ob["redundancy"][0]), 4),
                 "q_rtr(prod)": round(q["redundancy"], 3),
                 "q_rtr(vanilla)": round(qv["redundancy"], 3),
                 "q_sts(prod)": round(q["synergy"], 3)})
df = pd.DataFrame(rows)
print(df.to_string(index=False))
from scipy.stats import spearmanr  # noqa: E402
for col in ("q_rtr(prod)", "q_rtr(vanilla)"):
    r_s = spearmanr(df.rho, df[col]).statistic
    print(f"  monotonicity rho vs {col}: Spearman={r_s:+.2f} "
          f"({'OK' if r_s > 0.9 else 'CHECK'})")

print("\n" + "=" * 74)
print("1b-c  STRUCTURE-TRACKING: cross-coupling -> synergy measures")
print("=" * 74)
rows = []
for c in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
    A = np.array([[0.5, c], [0.8 * c, 0.4]])
    V = np.eye(2)
    ob, q = numit_of(A, V, POOL_PROD, seed=int(c * 100) + 31)
    _, qv = numit_of(A, V, POOL_VAN, seed=int(c * 100) + 31)
    rows.append({"c": c, "tdmi": round(float(ob["tdmi"][0]), 3),
                 "raw_sts": round(float(ob["synergy"][0]), 4),
                 "q_sts(prod)": round(q["synergy"], 3),
                 "q_syntot(prod)": round(q["synergy_total"], 3),
                 "q_syntot(van)": round(qv["synergy_total"], 3),
                 "q_rtr(prod)": round(q["redundancy"], 3)})
df = pd.DataFrame(rows)
print(df.to_string(index=False))
for col in ("q_sts(prod)", "q_syntot(prod)", "q_syntot(van)"):
    r_s = spearmanr(df.c, df[col]).statistic
    print(f"  monotonicity c vs {col}: Spearman={r_s:+.2f} "
          f"({'OK' if r_s > 0.9 else 'CHECK'})")
print("\nNOTE: c=0 row shows the documented MMI quirk (independent symmetricish")
print("system -> sts absorbs self-storage); interpret sweeps, not the origin.")
