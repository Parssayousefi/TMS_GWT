# V1 addendum: rotational (antisymmetric) coupling is the canonical
# synergy-GENERATING structure in 2-node Gaussian VARs (future depends on
# the joint phase of the pair). NuMIT synergy should track it.
import os, sys
import numpy as np, pandas as pd
BASE = r"D:\projects\TMS_GWT\numit-phiid\numit"
sys.path.insert(0, BASE); sys.path.insert(0, os.path.join(BASE, "validation"))
import phiid_gauss as pg, numit_null as nn
from scipy.stats import spearmanr
ws = pd.read_csv(r"D:\projects\TMS_GWT\pipeline\results\group_analysis\hcp\rest\hcp_rest_workspace_summary.csv")
ws_idx = np.where(ws.workspace_member.astype(bool).values)[0]
ts = np.load(r"D:\projects\TMS_GWT\pipeline\results\timeseries\dataset-basel_lpm_sub-001_ses-placebo_task-rest_run-1_timeseries_deconv.npy")
iu = np.triu_indices(len(ws_idx), 1)
pairs = np.stack([ws_idx[iu[0]], ws_idx[iu[1]]], axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    A_scan, _, _ = pg.var1_from_sigma(pg.lagged_sigma_from_timeseries(ts, 1, pairs))
POOL = A_scan[np.all(np.isfinite(A_scan.reshape(len(pairs), 4)), axis=1)]
rows = []
for c in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
    A = np.array([[0.4, c], [-c, 0.4]])
    s4, _ = pg.sigma4_from_var1(A[None], np.eye(2)[None], lag=1)
    ob = pg.phiid_atoms(pg.normalize_sigma(s4))
    r = np.random.default_rng(int(c * 100) + 3)
    with np.errstate(invalid="ignore", divide="ignore"):
        nl = nn.null_atoms_for_targets(POOL, ob["tdmi"], 400, r)
    q = {m: float(nn.quantile_of(ob[m], nl[m])[0]) for m in nn.MEASURES}
    rows.append({"c_rot": c, "tdmi": round(float(ob["tdmi"][0]), 3),
                 "raw_sts": round(float(ob["synergy"][0]), 4),
                 "q_sts": round(q["synergy"], 3),
                 "q_syntot": round(q["synergy_total"], 3),
                 "q_rtr": round(q["redundancy"], 3)})
df = pd.DataFrame(rows)
print(df.to_string(index=False))
for col in ("q_sts", "q_syntot"):
    r_s = spearmanr(df.c_rot, df[col]).statistic
    print(f"monotonicity c_rot vs {col}: Spearman={r_s:+.2f} "
          f"({'OK' if r_s > 0.9 else 'CHECK'})")
