# Vectorized analytic Gaussian MMI-PhiID from batched 4x4 lagged covariance
# matrices. Faithful port of the analytic (T->infinity) limit of Luppi/Mediano
# PhiIDFull.m + private/DoubleRedundancyMMI.m (tmp/matlab_deps copies, read
# 2026-08-25):
#   - variable order in the joint covariance: [x_past, y_past, x_fut, y_fut]
#   - single-target MMI redundancy = min of the marginal MIs
#   - double-redundancy rtr       = min of the four cross past-future MIs
#   - atoms = Minv @ reds with the constant 16x16 lattice matrix M
# PhiIDFull averages LOCAL (pointwise) atoms with plug-in Gaussian parameters;
# this module computes the closed-form expectation (its large-T limit).
# TEMPORARY (tmp-only investigation) - delete when the investigation closes.
import numpy as np

# 16x16 lattice matrix, copied verbatim from PhiIDFull.m (rows = reds, cols =
# atoms in order rtr,rtx,rty,rts,xtr,xtx,xty,xts,ytr,ytx,yty,yts,str,stx,sty,sts)
_M = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # rtr
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Rxyta
    [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Rxytb
    [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Rxytab
    [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Rabtx
    [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],  # Rabty
    [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],  # Rabtxy
    [1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Ixta
    [1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Ixtb
    [1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],  # Iyta
    [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],  # Iytb
    [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],  # Ixyta
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],  # Ixytb
    [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # Ixtab
    [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],  # Iytab
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],  # Ixytab
], dtype=np.float64)
_MINV = np.linalg.inv(_M)

ATOM_NAMES = ["rtr", "rtx", "rty", "rts", "xtr", "xtx", "xty", "xts",
              "ytr", "ytx", "yty", "yts", "str", "stx", "sty", "sts"]


def _logdet(sig, idx):
    """Batched log-determinant of the submatrix over index list idx.

    sig: (..., 4, 4), assumed SYMMETRIC (all our covariances are).
    Returns (...,) log|sub|; nan where det <= 0. Closed-form cofactor
    expansions for k<=4 - much faster than batched LAPACK slogdet.
    """
    k = len(idx)
    sub = sig[..., idx, :][..., :, idx]
    if k == 1:
        d = sub[..., 0, 0]
    elif k == 2:
        d = sub[..., 0, 0] * sub[..., 1, 1] - sub[..., 0, 1] * sub[..., 1, 0]
    elif k == 3:
        a = sub[..., 0, 0]
        b = sub[..., 1, 1]
        c = sub[..., 2, 2]
        dd = sub[..., 0, 1]
        e = sub[..., 0, 2]
        f = sub[..., 1, 2]
        d = a * b * c + 2 * dd * e * f - a * f * f - b * e * e - c * dd * dd
    else:
        a = sub[..., 0, 0]
        b = sub[..., 1, 1]
        c = sub[..., 2, 2]
        dg = sub[..., 3, 3]
        p = sub[..., 0, 1]
        q = sub[..., 0, 2]
        r = sub[..., 0, 3]
        s = sub[..., 1, 2]
        t = sub[..., 1, 3]
        u = sub[..., 2, 3]
        m0 = (b * c * dg + 2 * s * t * u - b * u * u - c * t * t - dg * s * s)
        m1 = p * (c * dg - u * u) - s * (q * dg - u * r) + t * (q * u - c * r)
        m2 = p * (s * dg - u * t) - b * (q * dg - u * r) + t * (q * t - s * r)
        m3 = p * (s * u - c * t) - b * (q * u - c * r) + s * (q * t - s * r)
        d = a * m0 - p * m1 + q * m2 - r * m3
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(d > 0, np.log(np.maximum(d, 1e-300)), np.nan)


def phiid_atoms(sigma):
    """All 16 MMI-PhiID atoms from batched 4x4 lagged covariances.

    sigma: (..., 4, 4) covariance of [x_past, y_past, x_fut, y_fut].
    Returns dict: 16 atoms + derived {tdmi, synergy, redundancy,
    synergy_total, phi, phi_r}, each shaped (...).
    """
    sigma = np.asarray(sigma, dtype=np.float64)
    ld = {}
    for idx in ([0], [1], [2], [3], [0, 1], [2, 3], [0, 2], [0, 3], [1, 2],
                [1, 3], [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3],
                [0, 1, 2, 3]):
        ld[tuple(idx)] = _logdet(sigma, list(idx))

    def mi(a, b):
        return 0.5 * (ld[tuple(a)] + ld[tuple(b)] - ld[tuple(sorted(a + b))])

    Ixta = mi([0], [2])
    Ixtb = mi([0], [3])
    Iyta = mi([1], [2])
    Iytb = mi([1], [3])
    Ixyta = mi([0, 1], [2])
    Ixytb = mi([0, 1], [3])
    Ixtab = mi([0], [2, 3])
    Iytab = mi([1], [2, 3])
    Ixytab = mi([0, 1], [2, 3])

    rtr = np.minimum(np.minimum(Ixta, Ixtb), np.minimum(Iyta, Iytb))
    reds = np.stack([
        rtr,
        np.minimum(Ixta, Iyta),        # Rxyta
        np.minimum(Ixtb, Iytb),        # Rxytb
        np.minimum(Ixtab, Iytab),      # Rxytab
        np.minimum(Ixta, Ixtb),        # Rabtx  (sources t1,t2 -> target x_past)
        np.minimum(Iyta, Iytb),        # Rabty
        np.minimum(Ixyta, Ixytb),      # Rabtxy
        Ixta, Ixtb, Iyta, Iytb,
        Ixyta, Ixytb, Ixtab, Iytab, Ixytab,
    ], axis=-1)                        # (..., 16)

    partials = reds @ _MINV.T          # (..., 16)
    out = {name: partials[..., i] for i, name in enumerate(ATOM_NAMES)}
    out["tdmi"] = Ixytab
    out["synergy"] = out["sts"]
    out["redundancy"] = out["rtr"]
    out["synergy_total"] = out["str"] + out["stx"] + out["sty"] + out["sts"]
    out["phi"] = Ixytab - Ixta - Iytb
    out["phi_r"] = out["phi"] + out["rtr"]
    return out


def lagged_sigma_from_timeseries(ts, lag, pairs):
    """Per-pair 4x4 lagged correlation matrices from one scan.

    Mirrors PhiIDFull's estimation: past block and future block are each
    z-scored over their own T-lag samples, covariance from the sample.
    ts: (T, R) array; pairs: (P, 2) int array of ROI indices.
    Returns sigma (P, 4, 4).
    """
    xp = ts[:-lag].astype(np.float64)
    xf = ts[lag:].astype(np.float64)
    n = xp.shape[0]
    xp = (xp - xp.mean(0)) / xp.std(0)
    xf = (xf - xf.mean(0)) / xf.std(0)
    c_pp = xp.T @ xp / n
    c_ff = xf.T @ xf / n
    c_pf = xp.T @ xf / n

    i, j = pairs[:, 0], pairs[:, 1]
    P = len(i)
    sig = np.empty((P, 4, 4))
    sig[:, 0, 0] = 1.0
    sig[:, 1, 1] = 1.0
    sig[:, 2, 2] = 1.0
    sig[:, 3, 3] = 1.0
    sig[:, 0, 1] = sig[:, 1, 0] = c_pp[i, j]
    sig[:, 2, 3] = sig[:, 3, 2] = c_ff[i, j]
    sig[:, 0, 2] = sig[:, 2, 0] = c_pf[i, i]   # x_p , x_f
    sig[:, 0, 3] = sig[:, 3, 0] = c_pf[i, j]   # x_p , y_f
    sig[:, 1, 2] = sig[:, 2, 1] = c_pf[j, i]   # y_p , x_f
    sig[:, 1, 3] = sig[:, 3, 1] = c_pf[j, j]   # y_p , y_f
    return sig


# ---------------- VAR(1) helpers (for NuMIT nulls) ----------------

def var1_from_sigma(sigma):
    """Per-pair 2-dim VAR(1) params implied by a 4x4 lagged covariance.

    sigma: (..., 4, 4) of [x_p, y_p, x_f, y_f]. Returns (A, V, S0):
    A (..., 2, 2) with future = A @ past; V residual covariance; S0
    contemporaneous covariance (past block).
    """
    S0 = sigma[..., :2, :2]
    C = sigma[..., :2, 2:]                 # cov(past, future): C[a, b] = cov(p_a, f_b)
    # future = A past  =>  cov(f, p) = A S0  =>  A = C^T S0^{-1}
    A = np.swapaxes(C, -1, -2) @ np.linalg.inv(S0)
    V = sigma[..., 2:, 2:] - A @ S0 @ np.swapaxes(A, -1, -2)
    return A, V, S0


def stationary_sigma0(A, V):
    """Solve the 2-dim discrete Lyapunov equation S0 = A S0 A^T + V, batched.

    vec(S0) = (I4 - kron(A, A))^{-1} vec(V), with vec in row-major order.
    A, V: (..., 2, 2). Returns S0 (..., 2, 2).
    """
    kron = np.einsum('...ij,...kl->...ikjl', A, A).reshape(A.shape[:-2] + (4, 4))
    eye = np.eye(4)
    rhs = V.reshape(V.shape[:-2] + (4,))
    s0 = np.linalg.solve(eye - kron, rhs[..., None])[..., 0]
    return s0.reshape(A.shape[:-2] + (2, 2))


def sigma4_from_var1(A, V, lag=1):
    """Joint 4x4 covariance [past, future] of a stationary 2-dim VAR(1).

    future is `lag` steps ahead: cov(p, f) = S0 (A^lag)^T.
    Returns (sigma4, S0). Not normalized to unit variance - callers that
    compare against correlation-based atoms should normalize with
    normalize_sigma().
    """
    S0 = stationary_sigma0(A, V)
    Al = A if lag == 1 else np.linalg.matrix_power(A, lag) if A.ndim == 2 \
        else _batched_matrix_power(A, lag)
    C = S0 @ np.swapaxes(Al, -1, -2)       # cov(past, future)
    top = np.concatenate([S0, C], axis=-1)
    bot = np.concatenate([np.swapaxes(C, -1, -2), S0], axis=-1)
    return np.concatenate([top, bot], axis=-2), S0


def _batched_matrix_power(A, k):
    out = A
    for _ in range(k - 1):
        out = out @ A
    return out


def normalize_sigma(sigma):
    """Rescale a covariance to a correlation matrix (unit diagonal), batched."""
    d = np.sqrt(np.maximum(np.diagonal(sigma, axis1=-2, axis2=-1), 1e-300))
    return sigma / (d[..., :, None] * d[..., None, :])


def tdmi_from_sigma4(sigma4):
    """TDMI of a 4x4 joint [past, future] covariance, batched."""
    return 0.5 * (_logdet(sigma4, [0, 1]) + _logdet(sigma4, [2, 3])
                  - _logdet(sigma4, [0, 1, 2, 3]))
