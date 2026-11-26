# V3 addendum: is the finite-sample synergy-quantile bias TDMI-dependent?
# If flat in TDMI, drug-induced TDMI shifts cannot leak into paired contrasts.
import sys
import numpy as np, pandas as pd
BASE = r"D:\projects\TMS_GWT\numit-phiid\numit"
sys.path.insert(0, BASE)
import phiid_gauss as pg, numit_null as nn
rng = np.random.default_rng(29)
ws = pd.read_csv(r"D:\projects\TMS_GWT\pipeline\results\group_analysis\hcp\rest\hcp_rest_workspace_summary.csv")
ws_idx = np.where(ws.workspace_member.astype(bool).values)[0]
ts = np.load(r"D:\projects\TMS_GWT\pipeline\results\timeseries\dataset-basel_lpm_sub-001_ses-placebo_task-rest_run-1_timeseries_deconv.npy")
iu = np.triu_indices(len(ws_idx), 1)
pairs = np.stack([ws_idx[iu[0]], ws_idx[iu[1]]], axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    A_scan, _, _ = pg.var1_from_sigma(pg.lagged_sigma_from_timeseries(ts, 1, pairs))
POOL = A_scan[np.all(np.isfinite(A_scan.reshape(len(pairs), 4)), axis=1)]
B, T = 3000, 350
tgt = rng.uniform(0.15, 1.8, B)
A = POOL[rng.integers(0, len(POOL), size=B)]
d, P, Pinv, detP = nn._eig2(A)
rho = np.maximum(np.abs(d).max(axis=-1), 1e-12)
c11 = np.sqrt(rng.chisquare(2, B)); c22 = np.sqrt(rng.chisquare(1, B)); n21 = rng.normal(0, 1, B)
V = np.empty((B, 2, 2)); V[:, 0, 0] = c11**2; V[:, 0, 1] = V[:, 1, 0] = c11*n21; V[:, 1, 1] = n21**2 + c22**2
Ahat = A / rho[:, None, None]
lo = np.full(B, 1e-4); hi = np.full(B, 1 - 1e-6)
for _ in range(40):
    mid = 0.5*(lo+hi)
    s4, _ = pg.sigma4_from_var1(Ahat*mid[:, None, None], V, lag=1)
    up = pg.tdmi_from_sigma4(s4) < tgt
    lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
Atrue = Ahat*(0.5*(lo+hi))[:, None, None]
est = np.empty((B, 4, 4))
for k in range(B):
    e = rng.multivariate_normal([0, 0], V[k], size=T+100)
    x = np.zeros((T+100, 2))
    for t in range(1, T+100): x[t] = Atrue[k] @ x[t-1] + e[t]
    est[k] = pg.lagged_sigma_from_timeseries(x[100:], 1, np.array([[0, 1]]))[0]
with np.errstate(invalid="ignore", divide="ignore"):
    obs = pg.phiid_atoms(est)
    ok = np.isfinite(obs["tdmi"]) & (obs["tdmi"] > 1e-6)
    nulls = nn.null_atoms_for_targets(POOL, obs["tdmi"][ok], 100, np.random.default_rng(5))
q = nn.quantile_of(obs["synergy"][ok], nulls["synergy"])
td = obs["tdmi"][ok]
bins = np.quantile(td, np.linspace(0, 1, 7))
lab = np.digitize(td, bins[1:-1])
rows = [{"tdmi_bin": f"{bins[i]:.2f}-{bins[i+1]:.2f}",
         "mean_q_sts": round(float(np.nanmean(q[lab == i])), 3),
         "n": int((lab == i).sum())} for i in range(6)]
print(pd.DataFrame(rows).to_string(index=False))
import scipy.stats as sps
r = sps.spearmanr(td, q[np.isfinite(q)] if np.isfinite(q).all() else q).statistic
print(f"\nSpearman(TDMI, q_sts) = {r:+.3f}  (near 0 = bias flat in TDMI)")
