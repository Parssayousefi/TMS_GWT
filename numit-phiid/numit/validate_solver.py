# B2 validation gate for the analytic Gaussian MMI-PhiID solver (phiid_gauss).
# Must pass before any NuMIT null runs. TEMPORARY - tmp-only investigation.
#
# 1. vs PRODUCTION: recompute tau=1 atoms from deconvolved timeseries for 4
#    scans; per-pair correlation vs results/phiid_matrices/*.npy
#    (accept r >= 0.999, median |rel diff| <= 1e-2).
# 2. vs phyid (tmp copy): agreement on synthetic VAR(1) data.
# 3. VAR wiring: atoms from the direct VAR(1) reconstruction of a data sigma
#    reproduce the data-sigma atoms exactly.
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phiid_gauss as pg  # noqa: E402

PIPE = r"D:\projects\TMS_GWT\pipeline"
TS = os.path.join(PIPE, "results", "timeseries")
PHIID = os.path.join(PIPE, "results", "phiid_matrices")
RNG = np.random.default_rng(2026)

SCANS = [
    "dataset-monash_psilo_sub-PC001_ses-baseline_task-rest_run-1",
    "dataset-monash_psilo_sub-PC001_ses-psilocybin_task-rest_run-1",
    "dataset-basel_lpm_sub-001_ses-psilocybin_task-rest_run-1",
    "dataset-berlin_ketamine_sub-01_ses-ketamine_task-rest_run-1",
]
MEASURES = ["synergy", "redundancy", "phi", "phi_r", "synergy_total"]

print("=" * 70)
print("CHECK 1: analytic solver vs production step04 matrices (tau=1)")
print("=" * 70)
all_ok = True
for scan in SCANS:
    ts_f = os.path.join(TS, scan + "_timeseries_deconv.npy")
    if not os.path.exists(ts_f):
        # berlin subject ids may differ - find any berlin ketamine scan
        import glob
        cands = sorted(glob.glob(os.path.join(
            TS, "dataset-berlin_ketamine_sub-*_ses-ketamine_*_timeseries_deconv.npy")))
        if "berlin" in scan and cands:
            ts_f = cands[0]
            scan = os.path.basename(ts_f).replace("_timeseries_deconv.npy", "")
        else:
            print(f"  SKIP (missing): {scan}")
            continue
    ts = np.load(ts_f)
    R = ts.shape[1]
    valid = np.where(ts.std(0) > 0)[0]
    # random subset of valid-ROI pairs for speed
    pk = RNG.choice(len(valid), size=(6000, 2))
    pk = pk[pk[:, 0] != pk[:, 1]][:3000]
    pairs = valid[pk]
    sigma = pg.lagged_sigma_from_timeseries(ts, 1, pairs)
    atoms = pg.phiid_atoms(sigma)
    print(f"\n  {scan}  (n_pairs={len(pairs)})")
    for meas in MEASURES:
        prod = np.load(os.path.join(PHIID, scan + f"_{meas}.npy"))
        pv = prod[pairs[:, 0], pairs[:, 1]]
        av = atoms[meas]
        ok = np.isfinite(pv) & np.isfinite(av)
        r = np.corrcoef(pv[ok], av[ok])[0, 1]
        scale = np.median(np.abs(pv[ok])) + 1e-12
        med_rel = np.median(np.abs(av[ok] - pv[ok])) / scale
        flag = "OK " if (r >= 0.999 and med_rel <= 1e-2) else "FAIL"
        if flag == "FAIL":
            all_ok = False
        print(f"    {meas:14s} r={r:.6f}  med|diff|/scale={med_rel:.2e}  "
              f"finite={ok.mean():.3f}  {flag}")

print()
print("=" * 70)
print("CHECK 2: analytic solver vs phyid on synthetic VAR(1) data")
print("=" * 70)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phyid_copy.calculate import calc_PhiID  # noqa: E402

for trial in range(3):
    # random stable 2-dim VAR(1)
    while True:
        A = RNG.normal(0, 0.5, (2, 2))
        if np.max(np.abs(np.linalg.eigvals(A))) < 0.95:
            break
    V = np.array([[1.0, 0.3], [0.3, 1.0]])
    T = 200_000
    x = np.zeros((T, 2))
    eps = RNG.multivariate_normal([0, 0], V, size=T)
    for t in range(1, T):
        x[t] = A @ x[t - 1] + eps[t]
    x = x[1000:]

    atoms_ph, _ = calc_PhiID(x[:, 0], x[:, 1], tau=1, kind="gaussian",
                             redundancy="MMI")
    sigma = pg.lagged_sigma_from_timeseries(x, 1, np.array([[0, 1]]))
    mine = pg.phiid_atoms(sigma)
    # also from the TRUE VAR parameters (stationary)
    s4_true, _ = pg.sigma4_from_var1(A[None], V[None], lag=1)
    truth = pg.phiid_atoms(pg.normalize_sigma(s4_true))
    print(f"\n  trial {trial}: spectral radius="
          f"{np.max(np.abs(np.linalg.eigvals(A))):.3f}")
    worst = 0.0
    for name in pg.ATOM_NAMES:
        ph = float(np.mean(atoms_ph[name][np.isfinite(atoms_ph[name])]))
        an = float(mine[name][0])
        tr = float(truth[name][0])
        worst = max(worst, abs(ph - an))
        if name in ("sts", "rtr", "xtx", "yty"):
            print(f"    {name}: phyid={ph:+.5f} analytic={an:+.5f} "
                  f"true-VAR={tr:+.5f}")
    flag = "OK" if worst < 5e-3 else "FAIL"
    if flag == "FAIL":
        all_ok = False
    print(f"    max |phyid - analytic| over 16 atoms: {worst:.2e}  {flag}")

print()
print("=" * 70)
print("CHECK 3: VAR(1) wiring - direct reconstruction reproduces data sigma")
print("=" * 70)
ts = np.load(os.path.join(
    TS, "dataset-monash_psilo_sub-PC001_ses-baseline_task-rest_run-1"
        "_timeseries_deconv.npy"))
pairs = np.array([[0, 100], [5, 250], [400, 453]])
sigma = pg.lagged_sigma_from_timeseries(ts, 1, pairs)
A, V, S0 = pg.var1_from_sigma(sigma)
C = S0 @ np.swapaxes(A, -1, -2)
top = np.concatenate([S0, C], axis=-1)
bot = np.concatenate([np.swapaxes(C, -1, -2),
                      A @ S0 @ np.swapaxes(A, -1, -2) + V], axis=-1)
sigma_rec = np.concatenate([top, bot], axis=-2)
err = np.max(np.abs(sigma_rec - sigma))
a1 = pg.phiid_atoms(sigma)
a2 = pg.phiid_atoms(sigma_rec)
aerr = max(np.max(np.abs(a1[n] - a2[n])) for n in pg.ATOM_NAMES)
flag = "OK" if (err < 1e-10 and aerr < 1e-10) else "FAIL"
if flag == "FAIL":
    all_ok = False
print(f"  max|sigma_rec - sigma|={err:.2e}  max atom diff={aerr:.2e}  {flag}")
# informational: stationary-vs-sample gap (nulls use the stationary route)
s4_st, _ = pg.sigma4_from_var1(A, V, lag=1)
gap = np.max(np.abs(pg.normalize_sigma(s4_st) - sigma))
print(f"  (info) stationary-VAR vs sample sigma gap: {gap:.3f}")

print()
print("=" * 70)
print("GATE:", "PASS" if all_ok else "FAIL")
print("=" * 70)
