# Unit tests for numit_null: closed-form Lyapunov/TDMI vs brute force,
# TDMI-match accuracy, timing, quantile sanity on one real scan.
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phiid_gauss as pg  # noqa: E402
import numit_null as nn  # noqa: E402

rng = np.random.default_rng(7)

print("TEST 1: closed-form stationary S0 + TDMI vs brute-force Lyapunov")
worst = 0.0
for _ in range(200):
    while True:
        A = rng.normal(0, 0.6, (2, 2))
        ev = np.linalg.eigvals(A)
        if np.max(np.abs(ev)) < 0.97 and np.abs(ev[0] - ev[1]) > 1e-6:
            break
    V = np.array([[1.0, 0.4], [0.4, 2.0]])
    # brute force: S0 = sum A^k V A'^k
    S0_bf = np.zeros((2, 2))
    Ak = np.eye(2)
    for _k in range(3000):
        S0_bf += Ak @ V @ Ak.T
        Ak = Ak @ A
    # closed form via the module internals: use unit-radius scaling with
    # f = rho(A) to reproduce exactly A
    d, P, Pinv, _ = nn._eig2(A[None])
    rho = np.abs(d).max()
    d_hat = d / rho
    W = Pinv @ V[None].astype(np.complex128) @ np.conj(np.swapaxes(Pinv, -1, -2))
    dd = d_hat[..., :, None] * np.conj(d_hat[..., None, :])
    f = np.array([rho])
    M = W / (1.0 - (f * f)[..., None, None] * dd)
    S0_cf = (P @ M @ np.conj(np.swapaxes(P, -1, -2))).real[0]
    worst = max(worst, np.max(np.abs(S0_cf - S0_bf)))
print(f"  max |S0_closed - S0_brute| over 200 systems: {worst:.2e}  "
      f"{'OK' if worst < 1e-8 else 'FAIL'}")

print()
print("TEST 2: TDMI Schur identity  0.5 log(detS0/detV) == sigma4 TDMI")
A = rng.normal(0, 0.4, (50, 2, 2))
ev = np.abs(np.linalg.eigvals(A)).max(axis=1)
A = A / np.maximum(ev, 1.0)[:, None, None] * 0.8
V = np.tile(np.eye(2), (50, 1, 1)) + 0.3
s4, S0 = pg.sigma4_from_var1(A, V, lag=1)
td1 = pg.tdmi_from_sigma4(s4)
detS0 = np.linalg.det(S0)
detV = np.linalg.det(V)
td2 = 0.5 * np.log(detS0 / detV)
err = np.max(np.abs(td1 - td2))
print(f"  max diff: {err:.2e}  {'OK' if err < 1e-10 else 'FAIL'}")

print()
print("TEST 3: TDMI matching + timing on a real scan chunk")
TS = r"D:\projects\TMS_GWT\pipeline\results\timeseries"
ts = np.load(os.path.join(
    TS, "dataset-monash_psilo_sub-PC001_ses-baseline_task-rest_run-1"
        "_timeseries_deconv.npy"))
valid = np.where(ts.std(0) > 0)[0]
pk = rng.choice(len(valid), size=(9000, 2))
pk = pk[pk[:, 0] != pk[:, 1]][:8192]
pairs = valid[pk]
sigma = pg.lagged_sigma_from_timeseries(ts, 1, pairs)
obs = pg.phiid_atoms(sigma)
A_all, V_all, _ = pg.var1_from_sigma(sigma)
tgt = obs["tdmi"]
ok = np.isfinite(tgt) & (tgt > 1e-6)
print(f"  chunk: {ok.sum()} valid pairs; tdmi range "
      f"[{tgt[ok].min():.3f}, {tgt[ok].max():.3f}]")

t0 = time.time()
nulls = nn.null_atoms_for_targets(A_all[ok], tgt[ok], n_null=100, rng=rng)
dt = time.time() - t0
print(f"  100 nulls x {ok.sum()} pairs: {dt:.1f} s "
      f"({dt / ok.sum() * 102831:.1f} s extrapolated per full-brain scan)")
me = nulls["_tdmi_err"]
print(f"  TDMI match |err|: median={np.median(me):.2e} "
      f"p99={np.percentile(me, 99):.2e} max={me.max():.2e}")

q_syn = nn.quantile_of(obs["synergy"][ok], nulls["synergy"])
q_red = nn.quantile_of(obs["redundancy"][ok], nulls["redundancy"])
for name, q in (("synergy", q_syn), ("redundancy", q_red)):
    qo = q[np.isfinite(q)]
    print(f"  NuMIT-{name}: mean={qo.mean():.3f} sd={qo.std():.3f} "
          f"frac@0={np.mean(qo == 0):.3f} frac@1={np.mean(qo == 1):.3f}")
print("  (degenerate saturation would be frac@0/frac@1 near 1)")

print()
print("TEST 4: CALIBRATION - systems drawn from the null ensemble must give")
print("  Uniform(0,1) quantiles (the rigorous end-to-end machinery test)")
pool = A_all[ok][:4000]
n_probe = 1500
probe_tdmi = rng.uniform(0.2, 2.0, n_probe)
probes = nn.null_atoms_for_targets(pool, probe_tdmi, n_null=1, rng=rng)
nulls2 = nn.null_atoms_for_targets(pool, probe_tdmi, n_null=100, rng=rng)
from scipy import stats as sps  # noqa: E402
for name in ("synergy", "redundancy", "phi_r"):
    obs_p = probes[name][:, 0].astype(np.float64)
    q = nn.quantile_of(obs_p, nulls2[name])
    qo = q[np.isfinite(q)]
    ks = sps.kstest(qo, "uniform")
    flag = "OK" if ks.pvalue > 0.01 else "FAIL"
    print(f"  {name:12s} n={len(qo)} mean={qo.mean():.3f} "
          f"KS p={ks.pvalue:.3f}  {flag}")
