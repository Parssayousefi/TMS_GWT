"""Validated clinical-outcome input for TMS analyses (not the frozen ASC loader).

Use dataset/subject_id/comparison to join participant-level NEURAL CHANGES.
Default loading masks flagged measures, retaining all published values on disk.
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PATH = ROOT / 'data/master_clinical_effects.csv'
KEYS = ['dataset', 'subject_id', 'comparison']
METRICS = ('craving_vas', 'craving_ccq_general', 'craving_ccq_now', 'impulsivity_bis_total',
           'impulsivity_bis_cognitive', 'impulsivity_bis_motor', 'impulsivity_bis_nonplanning')


def subject_key(dataset, bids_subject_id):
    """Explicit mapping from SUDMEX BIDS identifiers to the GWT master convention."""
    value = str(bids_subject_id)
    if dataset != 'sudmex_tms' or not value.startswith('sub-') or not value[4:].isdigit():
        raise ValueError('Expected dataset sudmex_tms and a BIDS subject such as sub-001')
    return f'{dataset}_{int(value[4:]):03d}'


def validate_effects(frame):
    if frame.empty:
        raise ValueError('Clinical input has no participant rows')
    required = KEYS + ['subject_num', 'bids_subject_id', 'arm', 'baseline_session', 'followup_session',
                       'analysis_cohort', 'mri_pair_usable']
    required += [m + suffix for m in METRICS for suffix in ('_t0', '_t1', '_delta', '_qc')]
    missing = sorted(set(required) - set(frame))
    if missing:
        raise ValueError('Missing clinical columns: ' + ', '.join(missing))
    if frame[KEYS].isna().any().any() or frame.duplicated(KEYS).any():
        raise ValueError('Missing or duplicate clinical identity/contrast keys')
    if not frame.dataset.eq('sudmex_tms').all() or not frame.arm.isin(['active', 'sham']).all():
        raise ValueError('Unexpected dataset or randomized arm')
    if not frame.comparison.eq('t1_minus_t0').all() or not frame.baseline_session.eq('t0').all() or not frame.followup_session.eq('t1').all():
        raise ValueError('This loader expects the acute t1_minus_t0 contrast only')
    for _, row in frame.iterrows():
        if row.subject_id != subject_key(row.dataset, row.bids_subject_id) or row.subject_num != int(row.bids_subject_id[4:]):
            raise ValueError('Inconsistent subject identifiers')
    for col in ('analysis_cohort', 'mri_pair_usable'):
        if not frame[col].isin([0, 1, False, True]).all():
            raise ValueError('Invalid boolean column: ' + col)
    for metric in METRICS:
        cols = [metric + suffix for suffix in ('_t0', '_t1', '_delta')]
        values = frame[cols].apply(pd.to_numeric, errors='raise')
        if np.isinf(values.to_numpy(dtype=float)).any():
            raise ValueError('Infinite clinical score: ' + metric)
        expected = values.iloc[:, 1] - values.iloc[:, 0]
        if not np.allclose(values.iloc[:, 2], expected, equal_nan=True, atol=1e-10, rtol=1e-10):
            raise ValueError('Clinical change must equal t1 minus t0: ' + metric)
        if frame[metric + '_qc'].isna().any():
            raise ValueError('Missing source-quality flag: ' + metric)
        if (frame[metric + '_qc'].eq('ok') & values.isna().any(axis=1)).any():
            raise ValueError('A missing score is incorrectly marked ok: ' + metric)
    return frame


def load_clinical_effects(path=None, metric=None, require_mri_pair=True, include_flagged=False):
    """Read the master CSV or workbook; mask flagged deltas unless explicitly allowed.

    With metric set, return only complete observations for that outcome. Without
    metric, retain all eligible rows and use NaN for unavailable/flagged outcomes.
    No imputation, normalization, arm pooling or hypothesis tests are performed.
    """
    path = Path(path) if path is not None else DEFAULT_PATH
    if metric is not None and metric not in METRICS:
        raise ValueError('Unknown clinical metric: ' + str(metric))
    if path.suffix.lower() == '.xlsx':
        frame = pd.read_excel(path, sheet_name='Clinical changes', engine='openpyxl')
    elif path.suffix.lower() == '.csv':
        frame = pd.read_csv(path)
    else:
        raise ValueError('Clinical input must be CSV or XLSX')
    validate_effects(frame)
    frame = frame[frame.analysis_cohort == 1].copy()
    if require_mri_pair:
        frame = frame[frame.mri_pair_usable == 1].copy()
    for name in METRICS:
        flagged = frame[name + '_qc'] != 'ok'
        if flagged.any() and (metric is None or metric == name):
            action = 'retained by explicit request' if include_flagged else 'masked'
            warnings.warn(f'{name}: {int(flagged.sum())} source-flagged changes {action}', UserWarning, stacklevel=2)
        if not include_flagged:
            frame.loc[flagged, name + '_delta'] = np.nan
    if metric:
        frame = frame[frame[metric + '_delta'].notna()].copy()
    return frame.reset_index(drop=True)


def merge_clinical_neural(neural_changes, path=None, metric=None, include_flagged=False):
    """One-to-one join with participant-level neural changes, never individual scans.

    Neural input must have the same explicit three keys. The returned left join
    keeps neural observations with missing clinical data visible via _merge.
    """
    if not set(KEYS).issubset(neural_changes) or neural_changes[KEYS].isna().any().any() or neural_changes.duplicated(KEYS).any():
        raise ValueError('Neural changes require unique dataset/subject_id/comparison keys')
    clinical = load_clinical_effects(path, metric=metric, include_flagged=include_flagged)
    # If upstream supplied arms, check them before preserving only one arm column.
    if 'arm' in neural_changes:
        check = neural_changes[KEYS + ['arm']].merge(clinical[KEYS + ['arm']], on=KEYS, suffixes=('_neural', '_clinical'))
        if not check.arm_neural.eq(check.arm_clinical).all():
            raise ValueError('Neural/clinical randomized-arm mismatch')
        clinical = clinical.drop(columns='arm')
    overlap = (set(neural_changes) & set(clinical)) - set(KEYS)
    if overlap:
        raise ValueError('Ambiguous overlapping columns: ' + ', '.join(sorted(overlap)))
    return neural_changes.merge(clinical, on=KEYS, how='left', validate='one_to_one', indicator=True)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=DEFAULT_PATH)
    args = parser.parse_args()
    frame = load_clinical_effects(args.path)
    print(json.dumps({'input': str(args.path), 'mri_eligible_participants': len(frame),
                      'available_changes': {m: int(frame[m + '_delta'].notna().sum()) for m in METRICS}}, indent=2))
