# B3: NuMIT null machinery for pairwise Gaussian MMI-PhiID (data-driven
# variant, Liardi et al. PLoS Comput Biol 2025, extended from PID to PhiID).
#
# Null model for one ROI pair with observed lag-tau TDMI:
#   - draw an empirical 2x2 VAR(1) coefficient matrix A from the SCAN's pool
#     of per-pair As (data-driven anchoring; vanilla N(0,1) As are nothing
#     like autocorrelation-dominated BOLD),
#   - draw residual covariance V ~ Wishart(I2, df=2)  (V's scale is
#     irrelevant: VAR TDMI is invariant to noise scale, per the paper),
#   - rescale A to spectral radius f, with f found by bisection so the null's
#     TDMI equals the observed TDMI,
#   - stationary S0 from the closed-form eigen solution of the 2-dim
#     Lyapunov equation (no linear solves needed):
#       A = P diag(d) P^-1,  W = P^-1 V P^-H,
#       S0(f) = P M(f) P^H with M_ij = W_ij / (1 - f^2 d_i conj(d_j))
#     (A scaled to unit spectral radius first, so A_f = f*A_hat).
#   - TDMI(f) = 0.5 log( det S0(f) / det V )   [Schur complement identity]
#   - joint sigma4 = [[S0, S0 A_f'], [A_f S0, S0]] -> normalize -> atoms.
#
# NuMIT value of an atom = quantile of the observed atom among N nulls.
# TEMPORARY - tmp-only investigation.
import numpy as np

import phiid_gauss as pg

MEASURES = ["synergy", "redundancy", "synergy_total", "phi", "phi_r"]


def _eig2(A):
    """Eigendecomposition of batched real 2x2 matrices (complex output).

    Returns d (..., 2) eigenvalues, P (..., 2, 2) eigenvector columns,
    Pinv (..., 2, 2). Near-defective matrices get a tiny jitter upstream.
    """
    a = A[..., 0, 0]
    b = A[..., 0, 1]
    c = A[..., 1, 0]
    e = A[..., 1, 1]
    tr = a + e
    det = a * e - b * c
    disc = np.sqrt((tr * tr - 4 * det).astype(np.complex128))
    d1 = (tr + disc) / 2
    d2 = (tr - disc) / 2
    # eigenvector for lambda: [b, lambda - a] unless b ~ 0, then [lambda-e, c]
    use_b = np.abs(b) > 1e-12
    v1 = np.where(use_b, b, d1 - e)
    w1 = np.where(use_b, d1 - a, c)
    v2 = np.where(use_b, b, d2 - e)
    w2 = np.where(use_b, d2 - a, c)
    # normalize columns
    n1 = np.sqrt(np.abs(v1) ** 2 + np.abs(w1) ** 2) + 1e-300
    n2 = np.sqrt(np.abs(v2) ** 2 + np.abs(w2) ** 2) + 1e-300
    P = np.empty(A.shape[:-2] + (2, 2), dtype=np.complex128)
    P[..., 0, 0] = v1 / n1
    P[..., 1, 0] = w1 / n1
    P[..., 0, 1] = v2 / n2
    P[..., 1, 1] = w2 / n2
    detP = P[..., 0, 0] * P[..., 1, 1] - P[..., 0, 1] * P[..., 1, 0]
    Pinv = np.empty_like(P)
    Pinv[..., 0, 0] = P[..., 1, 1] / detP
    Pinv[..., 0, 1] = -P[..., 0, 1] / detP
    Pinv[..., 1, 0] = -P[..., 1, 0] / detP
    Pinv[..., 1, 1] = P[..., 0, 0] / detP
    d = np.stack([d1, d2], axis=-1)
    return d, P, Pinv, np.abs(detP)


def null_atoms_for_targets(A_pool, target_tdmi, n_null, rng, chunk_atoms=True,
                           wishart_df=2):
    """Null-atom distributions for a batch of target TDMI values.

    A_pool: (Q, 2, 2) empirical VAR coefficient pool for this scan.
    target_tdmi: (P,) observed TDMI per pair (must be > 0).
    wishart_df: Bartlett degrees of freedom for V ~ Wishart(I2, df).
    Returns dict measure -> (P, n_null) float32 of null atom values.
    """
    P = target_tdmi.shape[0]
    B = P * n_null

    # --- sample A from pool + V ~ Wishart(I2, 2) via Bartlett ---
    A = A_pool[rng.integers(0, len(A_pool), size=B)]
    # jitter near-defective / near-zero matrices
    d, Pm, Pinv, detP = _eig2(A)
    bad = (detP < 1e-6) | (np.abs(d[..., 0] - d[..., 1]) < 1e-9)
    if bad.any():
        A = A.copy()
        A[bad] += rng.normal(0, 1e-4, (int(bad.sum()), 2, 2))
        d, Pm, Pinv, detP = _eig2(A)
    # scale to unit spectral radius
    rho = np.maximum(np.abs(d).max(axis=-1), 1e-12)
    d_hat = d / rho[..., None]

    c11 = np.sqrt(rng.chisquare(wishart_df, B))
    c22 = np.sqrt(rng.chisquare(wishart_df - 1, B))
    n21 = rng.normal(0, 1, B)
    V = np.empty((B, 2, 2))
    V[:, 0, 0] = c11 * c11
    V[:, 0, 1] = V[:, 1, 0] = c11 * n21
    V[:, 1, 1] = n21 * n21 + c22 * c22
    detV = V[:, 0, 0] * V[:, 1, 1] - V[:, 0, 1] ** 2
    # guard against degenerate V
    thin = detV < 1e-8
    if thin.any():
        V[thin, 0, 0] += 0.05
        V[thin, 1, 1] += 0.05
        detV = V[:, 0, 0] * V[:, 1, 1] - V[:, 0, 1] ** 2

    # W = Pinv V Pinv^H  (complex 2x2)
    W = Pinv @ V.astype(np.complex128) @ np.conj(np.swapaxes(Pinv, -1, -2))

    # d_i * conj(d_j) outer product, (B, 2, 2)
    dd = d_hat[..., :, None] * np.conj(d_hat[..., None, :])

    tgt = np.repeat(target_tdmi, n_null)

    # det(S0(f)) = |det P|^2 * det(M(f))  with  M_ij = W_ij/(1 - f^2 dd_ij),
    # so the bisection needs only elementwise complex arithmetic.
    W00 = W[..., 0, 0]
    W01 = W[..., 0, 1]
    W10 = W[..., 1, 0]
    W11 = W[..., 1, 1]
    dd00 = dd[..., 0, 0]
    dd01 = dd[..., 0, 1]
    dd10 = dd[..., 1, 0]
    dd11 = dd[..., 1, 1]
    detP2 = detP * detP  # |det P|^2 (detP from _eig2 is already abs)

    def tdmi_det_only(f2):
        m00 = W00 / (1.0 - f2 * dd00)
        m01 = W01 / (1.0 - f2 * dd01)
        m10 = W10 / (1.0 - f2 * dd10)
        m11 = W11 / (1.0 - f2 * dd11)
        detS0 = (detP2 * (m00 * m11 - m01 * m10)).real
        with np.errstate(invalid="ignore", divide="ignore"):
            return 0.5 * np.log(detS0 / detV)

    # bisection on f in (0, fmax)
    lo = np.full(B, 1e-6)
    hi = np.full(B, 1.0 - 1e-6)
    # nulls whose TDMI at fmax is still below target cannot reach it -> they
    # saturate at fmax (counted in QC); extremely rare since TDMI -> inf.
    for _ in range(32):
        mid = 0.5 * (lo + hi)
        go_up = tdmi_det_only(mid * mid) < tgt
        lo = np.where(go_up, mid, lo)
        hi = np.where(go_up, hi, mid)
    f = 0.5 * (lo + hi)
    t_f = tdmi_det_only(f * f)
    M = W / (1.0 - (f * f)[..., None, None] * dd)
    S0 = (Pm @ M @ np.conj(np.swapaxes(Pm, -1, -2))).real

    # assemble joint sigma4 of the stationary VAR(1) null
    Af = (Pm @ (f[..., None, None] * np.eye(2) * d_hat[..., None, :])
          @ Pinv).real       # P diag(f*d) Pinv
    C = S0 @ np.swapaxes(Af, -1, -2)
    sig4 = np.empty((B, 4, 4))
    sig4[:, :2, :2] = S0
    sig4[:, :2, 2:] = C
    sig4[:, 2:, :2] = np.swapaxes(C, -1, -2)
    sig4[:, 2:, 2:] = S0
    sig4 = pg.normalize_sigma(sig4)

    with np.errstate(invalid="ignore", divide="ignore"):
        atoms = pg.phiid_atoms(sig4)

    out = {m: atoms[m].reshape(P, n_null).astype(np.float32)
           for m in MEASURES}
    out["_tdmi_err"] = np.abs(t_f - tgt).reshape(P, n_null).astype(np.float32)
    return out


def quantile_of(obs, null):
    """Quantile of obs (P,) within null (P, N): mean(null < obs) + 0.5 ties."""
    less = (null < obs[:, None]).sum(axis=1)
    eq = (null == obs[:, None]).sum(axis=1)
    n_ok = np.isfinite(null).sum(axis=1)
    q = (less + 0.5 * eq) / np.maximum(n_ok, 1)
    q[~np.isfinite(obs)] = np.nan
    q[n_ok < 10] = np.nan
    return q.astype(np.float32)
