"""Check saved TMS HRF outputs without rerunning the production pipeline."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def audit(source):
    prefix = source.name.removesuffix('_timeseries.npy')
    deconv = source.with_name(prefix + '_timeseries_deconv.npy')
    mat = source.with_name(prefix + '_hrf.mat')
    metadata = source.with_name(prefix + '_timeseries_deconv_metadata.yaml')
    x, y = np.load(source), np.load(deconv)
    if x.shape != y.shape or x.shape[1] != 454 or not np.isfinite(y).all():
        raise ValueError(f'Invalid deconvolution array: {source.name}')
    sigma = x.std(axis=0, ddof=1)
    valid = sigma >= 1e-10
    if np.any(y[:, ~valid] != 0) or np.any(y[:, valid].std(axis=0) == 0):
        raise ValueError(f'Invalid handling of constant ROIs: {source.name}')
    z = (x[:, valid] - x[:, valid].mean(axis=0)) / sigma[valid]
    if np.allclose(x, y) or np.allclose(z, y[:, valid]):
        raise ValueError(f'Output is only a copy or z-score: {source.name}')
    meta = yaml.safe_load(metadata.read_text(encoding='utf-8'))
    expected = {'tr': 2.0, 'basis_order': 3, 'thr': 1, 'hrf_len_s': 24, 'wiener_max_iter': 50}
    if any(meta.get(k) != v for k, v in expected.items()) or meta['deconv_shape'] != list(y.shape):
        raise ValueError(f'Incorrect saved parameters: {source.name}')
    with h5py.File(mat) as h:
        params = h['PARA'][:].T
        hrfs = h['hrfa'][:].T
        if params.shape != (3, 454) or hrfs.shape[1] != 454 or not np.isfinite(params).all() or not np.isfinite(hrfs).all():
            raise ValueError(f'Invalid HRF estimates: {source.name}')
        for field, expected_value in {'TR': 2., 'T': 3, 'order': 3, 'thr': 1, 'len': 24, 'AR_lag': 1}.items():
            if float(h['para'][field][0, 0]) != expected_value:
                raise ValueError(f'Incorrect MAT parameter {field}: {source.name}')
        if not np.array_equal(h['para']['lag'][:].ravel(), np.arange(6, 13)):
            raise ValueError('Incorrect onset lag grid')
        events = [h[r].size for r in h['event_bold'][:].ravel()]
    return {'file': source.name, 'n_timepoints': x.shape[0], 'n_rois': x.shape[1],
            'zero_variance_input_rois': int((~valid).sum()), 'output_std': float(y.std()),
            'rmse_vs_input_zscore': float(np.sqrt(np.mean((y[:, valid] - z)**2))),
            'ttp_min_s': float(params[1, valid].min()), 'ttp_max_s': float(params[1, valid].max()),
            'ttp_3_to_10_fraction': float(np.mean((params[1, valid] >= 3) & (params[1, valid] <= 10))),
            'ttp_outside_3_to_10_rois_one_based': ','.join(map(str, np.flatnonzero(valid & ((params[1] < 3) | (params[1] > 10))) + 1)),
            'mean_event_count': float(np.mean(events)),
            'hrf_peak_std_across_rois': float(params[0, valid].std()),
            'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'deconv_sha256': hashlib.sha256(deconv.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', choices=['aachen_itbs', 'sudmex_tms'], default='aachen_itbs')
    parser.add_argument('--out', type=Path, default=ROOT / 'logs/tms_pipeline_20261125/hrf_verification')
    args = parser.parse_args()
    sources = sorted((ROOT / 'pipeline/results/timeseries').glob(f'dataset-{args.dataset}_*_timeseries.npy'))
    if not sources:
        raise ValueError('No source scans found')
    rows = [audit(source) for source in sources]
    if len({r['deconv_sha256'] for r in rows}) != len(rows):
        raise ValueError('Duplicate outputs across scans')
    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out / (args.dataset + '_output_audit.csv'), index=False)
    inputs = args.out / 'inputs.json'
    checked_inputs = False
    if inputs.exists():
        before = {Path(r['file']).name: r['sha256'] for r in json.loads(inputs.read_text())}
        checked_inputs = all(r['file'] in before for r in rows)
        if checked_inputs and any(r['source_sha256'] != before[r['file']] for r in rows):
            raise ValueError('Input changed since preflight')
    summary = {'checked_utc': datetime.now(timezone.utc).isoformat(), 'dataset': args.dataset, 'scans_passed': len(rows),
               'all_sidecars_present': True, 'all_outputs_finite': True, 'all_outputs_differ_from_raw_and_zscore': True,
               'input_arrays_unchanged_verified': checked_inputs, 'timepoint_counts': sorted({r['n_timepoints'] for r in rows}),
               'scans_with_ttp_outside_3_to_10': sum(bool(r['ttp_outside_3_to_10_rois_one_based']) for r in rows),
               'minimum_physiological_ttp_fraction': min(r['ttp_3_to_10_fraction'] for r in rows),
               'ttp_range_s': [min(r['ttp_min_s'] for r in rows), max(r['ttp_max_s'] for r in rows)],
               'mean_events_range': [min(r['mean_event_count'] for r in rows), max(r['mean_event_count'] for r in rows)],
               'note': 'Parameter and numerical integrity checks; not proof of accurate latent neural activity.'}
    (args.out / (args.dataset + '_output_audit.json')).write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
