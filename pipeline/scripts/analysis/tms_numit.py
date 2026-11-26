"""TMS-only conditional-CDF evaluation for the exploratory NuMIT/PhiID branch.

The imported mathematical helpers are read-only. No cohort driver is imported.
References are shared across visits/arms; neither chunks nor scan order select
the empirical coefficient pool. Percentiles are not additive information atoms.
"""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'numit-phiid/numit'))
import numit_null as nn
import phiid_gauss as pg


def stable_phiid_atoms(sigma):
    """Independent LAPACK log determinants for the calibration observations.

    Explicit cofactor expansion can lose precision for nearly singular Wishart
    draws. This slower calculation retains those draws in the calibration.
    The redundancy definition and the fixed 16-atom lattice are unchanged.
    """
    sigma = np.asarray(sigma, dtype=float)
    def mi(a, b):
        def ld(indices):
            block = sigma[..., indices, :][..., :, indices]
            sign, value = np.linalg.slogdet(block)
            if np.any(sign <= 0): raise ValueError('Non-positive calibration covariance')
            return value
        return .5*(ld(a)+ld(b)-ld(sorted(a+b)))
    xa,xb,ya,yb = [mi([i],[j]) for i,j in [(0,2),(0,3),(1,2),(1,3)]]
    xya,xyb = mi([0,1],[2]),mi([0,1],[3])
    xab,yab = mi([0],[2,3]),mi([1],[2,3])
    total = mi([0,1],[2,3])
    redundancy = np.minimum.reduce([xa,xb,ya,yb])
    reds = np.stack([redundancy,np.minimum(xa,ya),np.minimum(xb,yb),np.minimum(xab,yab),
                     np.minimum(xa,xb),np.minimum(ya,yb),np.minimum(xya,xyb),
                     xa,xb,ya,yb,xya,xyb,xab,yab,total],axis=-1)
    atoms = reds@pg._MINV.T
    return {'tdmi':total,'synergy':atoms[...,15],'redundancy':atoms[...,0],
            'synergy_total':atoms[...,12:16].sum(axis=-1),
            'phi':total-xa-yb,'phi_r':total-xa-yb+redundancy}


def valid_reference(pool, target, count, seed):
    """Return exactly count valid draws, rejecting the same rows for all atoms."""
    rng = np.random.default_rng(seed)
    collected = {m: [] for m in nn.MEASURES}
    needed, rejected, generated, maximum = count, 0, 0, 0.
    for attempt in range(20):
        size = needed + 64
        result = nn.null_atoms_for_targets(pool, np.array([target]), size, rng)
        error = result['_tdmi_err'][0]
        good = np.isfinite(error) & (error <= 1e-6)
        for m in nn.MEASURES:
            good &= np.isfinite(result[m][0])
        indices = np.flatnonzero(good)[:needed]
        generated += size
        rejected += int((~good).sum())
        if len(indices):
            maximum = max(maximum, float(error[indices].max()))
            for m in nn.MEASURES:
                collected[m].append(result[m][0, indices])
            needed -= len(indices)
        if needed == 0:
            return {m: np.concatenate(a) for m, a in collected.items()}, {
                'generated': generated, 'rejected': rejected,
                'max_kept_tdmi_error': maximum, 'replenishment_batches': attempt}
    raise RuntimeError('Unable to obtain a valid conditional reference')


def empirical_cdf(sorted_values, observed):
    """Midrank CDF, including exact ties; reference contains finite values only."""
    return (np.searchsorted(sorted_values, observed, side='left') +
            np.searchsorted(sorted_values, observed, side='right')) / (2 * len(sorted_values))


def conditional_cdf(grid, sorted_reference, tdmi, observed):
    """Interpolate CDFs at the observed atom, never interpolate atom quantiles."""
    tdmi, observed = np.broadcast_arrays(np.asarray(tdmi), np.asarray(observed))
    shape = tdmi.shape
    tdmi, observed = tdmi.ravel(), observed.ravel()
    if not (np.isfinite(tdmi).all() and np.isfinite(observed).all()):
        raise ValueError('Nonfinite CDF query')
    if tdmi.min() < grid[0] or tdmi.max() > grid[-1]:
        raise ValueError(f'TDMI outside reference [{grid[0]}, {grid[-1]}]')
    lower = np.clip(np.searchsorted(grid, tdmi, side='right') - 1, 0, len(grid) - 2)
    fraction = (tdmi - grid[lower]) / (grid[lower + 1] - grid[lower])
    answer = np.empty(len(tdmi))
    for k in np.unique(lower):
        index = np.flatnonzero(lower == k)
        low = empirical_cdf(sorted_reference[k], observed[index])
        high = empirical_cdf(sorted_reference[k + 1], observed[index])
        answer[index] = low + fraction[index] * (high - low)
    assert ((answer >= 0) & (answer <= 1)).all()
    return answer.reshape(shape)
