"""Matched workspace contrasts and haemodynamic associations in both TMS cohorts.

Reads released corrected summaries and saved HRF parameters only. No production
processing, mask selection or reference fitting is performed. All declared tests
are retained, including the previously known Aachen spatial and HRF results.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import hashlib

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'setup'), str(ROOT / 'pipeline/scripts/analysis')]
import analyse_corrected_tms as core
import analyse_tms_workspace_scope as ws

RUN = ROOT / 'output/corrected_tms_analysis/20261126T112054_284547Z'
MEASURES = list(ws.MEASURES)
PARAMETERS = ['height', 'peak_seconds', 'width_seconds']
SCOPES = ['whole_brain', 'within_workspace', 'outside_workspace', 'difference']
INPUTS = {}
ERRORS = []


def track(path):
    path = Path(path)
    INPUTS[str(path.relative_to(ROOT))] = core.tms.sha256(path)
    return path


def load(name):
    return pd.read_csv(track(RUN / name))


def fit(frame, y, predictors, term):
    """Existing HC3 estimator, checked against an independent matrix formula."""
    result = ws.hc3(frame, y, predictors, term)
    assert result['status'] == 'ok', (term, result)
    design = frame[predictors].to_numpy(float)
    values = frame[y].to_numpy(float)
    inverse = np.linalg.pinv(design)
    beta = inverse @ values
    residual = values - design @ beta
    leverage = np.einsum('ij,ji->i', design, inverse)
    scaled = inverse * (residual / (1 - leverage))[None, :]
    covariance = scaled @ scaled.T
    k = predictors.index(term)
    se = np.sqrt(covariance[k, k])
    df = len(values) - design.shape[1]
    p = 2 * stats.t.sf(abs(beta[k] / se), df)
    ci = beta[k] + np.array([-1, 1]) * stats.t.ppf(.975, df) * se
    expected = np.array([beta[k], se, p, *ci])
    actual = np.array([result['estimate'], result['se'], result['p'], result['ci_low'], result['ci_high']])
    np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-10)
    ERRORS.append(float(np.max(abs(actual - expected))))
    return {**result, 'df': df, 'max_leverage': float(leverage.max())}


def zscore(x):
    assert np.isfinite(x).all() and x.std(ddof=1) > 1e-12
    return (x - x.mean()) / x.std(ddof=1)


def extract_hrf():
    scans = load('raw_scan_metrics.csv')
    plan = load('raw_paired_changes.csv')
    keys = set()
    for p in plan.itertuples():
        keys.update([(p.subject_id, p.baseline_session), (p.subject_id, p.followup_session)])
    mask = core.tms.load_reference(ROOT)['workspace']
    assert mask.sum() == 235 and len(mask) == 454
    track(ROOT / 'pipeline/results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy')
    records = []
    for scan in scans.itertuples():
        if (scan.subject_id, scan.session) not in keys:
            continue
        path = track(core.matrix_root(scan.dataset) / 'pipeline/results/timeseries' / (scan.prefix + '_hrf.mat'))
        with h5py.File(path) as f:
            parameters = f['PARA'][:]
            assert parameters.shape == (454, 3) and np.isfinite(parameters).all()
            assert float(f['para']['TR'][0, 0]) == 2.
        row = {k: getattr(scan, k) for k in ['dataset', 'subject_id', 'session', 'arm', 'mean_fd']}
        for scope, selected in [('whole_brain', np.ones(454, bool)), ('within_workspace', mask), ('outside_workspace', ~mask)]:
            for k, name in enumerate(PARAMETERS):
                row[scope + '__' + name] = parameters[selected, k].mean()
        for name in PARAMETERS:
            row['difference__' + name] = row['outside_workspace__' + name] - row['within_workspace__' + name]
        records.append(row)
    scans = pd.DataFrame(records)
    assert len(scans) == 140 and not scans.duplicated(['subject_id', 'session']).any()
    old = load('sensitivities/hrf_parameter_scans.csv').set_index(['subject_id', 'session']).sort_index()
    cols = [c for c in old if '__' in c]
    got = scans.query("dataset == 'aachen_itbs'").set_index(['subject_id', 'session']).sort_index()
    np.testing.assert_allclose(got[cols], old[cols], rtol=1e-12, atol=1e-12)
    index = scans.set_index(['subject_id', 'session'])
    pairs = []
    for p in plan.itertuples():
        a, b = index.loc[(p.subject_id, p.baseline_session)], index.loc[(p.subject_id, p.followup_session)]
        row = {k: getattr(p, k) for k in ['dataset', 'subject_id', 'arm', 'comparison']}
        for scope in SCOPES:
            for parameter in PARAMETERS:
                key = scope + '__' + parameter
                row[key + '_baseline'] = a[key]
                row[key + '_delta'] = b[key] - a[key]
        pairs.append(row)
    pairs = pd.DataFrame(pairs)
    aachen = pairs.query("dataset == 'aachen_itbs'").groupby('subject_id').mean(numeric_only=True).sort_index()
    sudmex = pairs.query("dataset == 'sudmex_tms'").set_index('subject_id').sort_index()
    assert len(aachen) == 15 and len(sudmex) == 40
    return scans, {'aachen_itbs': aachen, 'sudmex_tms': sudmex}


def information(branch):
    pairs = load(branch + '_paired_changes.csv')
    aa = pairs.query("dataset == 'aachen_itbs'").groupby('subject_id').mean(numeric_only=True).sort_index()
    sud = pairs.query("dataset == 'sudmex_tms'").set_index('subject_id').sort_index()
    for frame in [aa, sud]:
        for measure in MEASURES:
            for suffix in ['baseline', 'followup', 'delta']:
                frame['difference__' + measure + '_' + suffix] = frame['outside_workspace__' + measure + '_' + suffix] - frame['within_workspace__' + measure + '_' + suffix]
        frame['intercept'] = 1.
    sud['active'] = sud.arm.eq('active').astype(float)
    return {'aachen_itbs': aa, 'sudmex_tms': sud}


def spatial_tests(branch, data):
    rows = []
    aa, sud = data['aachen_itbs'], data['sudmex_tms']
    matrix = aa[['difference__' + m + '_delta' for m in MEASURES]]
    exact, maximum = ws.exact_max_t(matrix)
    for k, measure in enumerate(MEASURES):
        metric = 'difference__' + measure
        label = {'branch': branch, 'measure': measure, 'contrast': 'outside_minus_within_change'}
        rows.append({**label, 'dataset': 'aachen_itbs', 'model': 'paired_change', **core.tms.paired_summary(aa[metric + '_delta']), 'p_exact': exact[k], 'p_maxT5': maximum[k]})
        rows.append({**label, 'dataset': 'aachen_itbs', 'model': 'motion_adjusted_change', **fit(aa, metric + '_delta', ['intercept', 'delta_fd'], 'intercept')})
        rows.append({**label, 'dataset': 'sudmex_tms', 'model': 'unadjusted_active_minus_sham', **core.tms.welch_summary(sud.loc[sud.active.eq(1), metric + '_delta'], sud.loc[sud.active.eq(0), metric + '_delta'])})
        for name, predictors in [('baseline_adjusted', ['intercept', 'active', metric + '_baseline']), ('baseline_motion_adjusted', ['intercept', 'active', metric + '_baseline', 'delta_fd'])]:
            rows.append({**label, 'dataset': 'sudmex_tms', 'model': name, **fit(sud, metric + '_followup', predictors, 'active')})
        for arm in ['active', 'sham']:
            subset = sud.loc[sud.arm.eq(arm)]
            rows.append({**label, 'dataset': 'sudmex_tms', 'model': 'descriptive_' + arm + '_change', **core.tms.paired_summary(subset[metric + '_delta'])})
    return rows


def association_tests(branch, dataset, data, hrf):
    assert data.index.equals(hrf.index)
    labels = [(s, s, p) for s in SCOPES for p in PARAMETERS]
    # A shared whole-brain predictor allows an actual test of differential
    # within-versus-outside slopes, rather than comparing two p-values.
    labels += [('difference', 'whole_brain', p) for p in PARAMETERS]
    results, influence, observations = [], [], []
    for outcome_scope, predictor_scope, parameter in labels:
        label = {'branch': branch, 'dataset': dataset, 'outcome_scope': outcome_scope, 'predictor_scope': predictor_scope, 'parameter': parameter}
        frame = pd.DataFrame({'y': data[outcome_scope + '__phi_delta'], 'intercept': 1., 'delta_fd': data.delta_fd, 'x': hrf[predictor_scope + '__' + parameter + '_delta'], 'y_baseline': data[outcome_scope + '__phi_baseline'], 'x_baseline': hrf[predictor_scope + '__' + parameter + '_baseline']})
        frame['x_z'] = zscore(frame.x)
        predictors = ['intercept', 'delta_fd', 'x_z']
        if dataset == 'sudmex_tms':
            frame['active'] = data.active
            predictors.append('active')
        results.append({**label, 'model': 'motion_arm_adjusted', **fit(frame, 'y', predictors, 'x_z')})
        results.append({**label, 'model': 'baseline_motion_arm_adjusted', **fit(frame, 'y', predictors + ['y_baseline', 'x_baseline'], 'x_z')})
        rho, p = stats.spearmanr(frame.x.round(12), frame.y.round(12))
        results.append({**label, 'model': 'unadjusted_rank_sensitivity', 'n': len(frame), 'rho': rho, 'p': p})
        if dataset == 'sudmex_tms':
            frame['active_x_z'] = frame.active * frame.x_z
            results.append({**label, 'model': 'active_by_hrf_interaction', **fit(frame, 'y', predictors + ['active_x_z'], 'active_x_z')})
            for arm, active in [('active', 1), ('sham', 0)]:
                sub = frame.loc[frame.active.eq(active)].copy()
                # Keep a shared cohort-SD scale for comparing the two arms.
                results.append({**label, 'model': arm + '_motion_adjusted', **fit(sub, 'y', ['intercept', 'delta_fd', 'x_z'], 'x_z')})
        if branch == 'raw':
            for sid in frame.index:
                sub = frame.drop(sid).copy()
                sub['x_z'] = zscore(sub.x)
                influence.append({**label, 'omitted_subject': sid, **fit(sub, 'y', predictors, 'x_z')})
        for sid, row in frame.iterrows():
            observations.append({**label, 'subject_id': sid, 'phi_change': row.y, 'hrf_change': row.x, 'delta_fd': row.delta_fd, 'active': row.get('active', np.nan)})
    return results, influence, observations


def conjunction(table, keys, model):
    """Intersection-union evidence: both cohorts, same direction; no pooled mean."""
    selected = table.loc[table.model.eq(model)]
    merged = selected.query("dataset == 'aachen_itbs'").merge(selected.query("dataset == 'sudmex_tms'"), on=['branch', *keys], suffixes=('_aachen', '_sudmex'), validate='one_to_one')
    merged['same_direction'] = np.sign(merged.estimate_aachen) == np.sign(merged.estimate_sudmex)
    merged['p_conjunction'] = np.where(merged.same_direction, np.maximum(merged.p_aachen, merged.p_sudmex), 1.)
    merged['p_holm_family'] = merged.groupby('branch').p_conjunction.transform(ws.holm)
    return merged


def main():
    out = ROOT / 'output/cross_dataset_spatial_hrf' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True)
    print('OUTPUT:', out, flush=True)
    core.save_json(out / 'analysis_specification.json', {
        'question': 'Do changes differ inside versus outside the fixed workspace across both TMS datasets, and how are they related to estimated haemodynamics?',
        'status': 'Exploratory follow-up. Aachen spatial and nine HRF associations were known; this is not preregistration.',
        'spatial_family': 'Five direct outside-minus-within contrasts per dataset, branch and model. Aachen exact maxT5 plus Holm5; SUDMEX baseline-and-motion ANCOVA active-minus-sham Holm5. Raw analysis and normalized sensitivity retained.',
        'hrf_family': 'Fifteen associations per dataset/branch/model: Phi at four scopes against matching HRF height, peak and width (12), plus outside-minus-within Phi against whole-brain HRF parameters (3). Holm15.',
        'hrf_adjustment': 'Motion change in both cohorts; treatment arm also in SUDMEX. Baseline Phi and HRF adjustment is a declared sensitivity. Treatment-arm interaction and stratified estimates retained, each with the same 15-member family.',
        'replication': 'Same-sign intersection-union p=max(two cohort p-values), Holm across all fifteen matching endpoints; no pooling of treatment estimands or subjects.',
        'influence': 'All 15 primary raw associations repeated for each participant omission; family correction repeated for every omission.',
        'limits': 'No sham in Aachen; shared HRF/deconvolution data can cause model coupling. Association is not mediation or independent vascular evidence. A significant compartment and a nonsignificant compartment do not establish a spatial difference.'})
    track(Path(__file__))
    hrf_scans, hrf = extract_hrf()
    hrf_scans.to_csv(out / 'hrf_scans.csv', index=False)
    for dataset, frame in hrf.items():
        frame.to_csv(out / (dataset + '_hrf_participant_changes.csv'))
    spatial, associations, influence, observations = [], [], [], []
    for branch in ['raw', 'normalized']:
        data = information(branch)
        spatial += spatial_tests(branch, data)
        for dataset in ['aachen_itbs', 'sudmex_tms']:
            result, loo, obs = association_tests(branch, dataset, data[dataset], hrf[dataset])
            associations += result
            influence += loo
            observations += obs
            print('Complete:', branch, dataset, flush=True)
    spatial = pd.DataFrame(spatial)
    spatial['p_holm5'] = spatial.groupby(['branch', 'dataset', 'model']).p.transform(ws.holm)
    spatial.to_csv(out / 'spatial_results.csv', index=False)
    associations = pd.DataFrame(associations)
    associations['p_holm15'] = associations.groupby(['branch', 'dataset', 'model']).p.transform(ws.holm)
    associations.to_csv(out / 'hrf_association_results.csv', index=False)
    keys = ['outcome_scope', 'predictor_scope', 'parameter']
    joint = conjunction(associations, keys, 'motion_arm_adjusted')
    joint.to_csv(out / 'hrf_cross_dataset_conjunction.csv', index=False)
    pd.DataFrame(observations).to_csv(out / 'hrf_observations.csv', index=False)
    loo = pd.DataFrame(influence)
    loo['p_holm15'] = loo.groupby(['dataset', 'omitted_subject']).p.transform(ws.holm)
    loo.to_csv(out / 'hrf_leave_one_out.csv', index=False)
    loo.groupby(['dataset', *keys]).agg(estimate_min=('estimate', 'min'), estimate_max=('estimate', 'max'), worst_p_holm15=('p_holm15', 'max'), significant_omissions=('p_holm15', lambda x: int((x < .05).sum()))).to_csv(out / 'hrf_influence_summary.csv')
    old = load('aachen_spatial_interactions.csv').set_index(['branch', 'metric'])
    got = spatial.query("dataset == 'aachen_itbs' and model == 'paired_change'").set_index(['branch', 'measure']).loc[old.index]
    np.testing.assert_allclose(got[['estimate', 'ci_low', 'ci_high', 'p', 'p_maxT5']], old[['estimate', 'ci_low', 'ci_high', 'p', 'p_maxT_5']], rtol=1e-10, atol=1e-12)
    old = load('haemodynamic_followup/association_results.csv').rename(columns={'scope': 'outcome_scope'}).set_index(['outcome_scope', 'parameter'])
    got = associations.query("dataset == 'aachen_itbs' and branch == 'raw' and model == 'motion_arm_adjusted' and outcome_scope == predictor_scope").set_index(['outcome_scope', 'parameter']).loc[old.index]
    np.testing.assert_allclose(got[['estimate', 'ci_low', 'ci_high', 'p']], old[['estimate', 'ci_low', 'ci_high', 'p']], rtol=1e-9, atol=1e-12)
    assert all(core.tms.sha256(ROOT / p) == digest for p, digest in INPUTS.items())
    core.save_json(out / 'verification.json', {'status': 'passed', 'saved_Aachen_HRF_scans_match': True, 'prior_Aachen_spatial_and_HRF_estimates_reproduced': True, 'independent_HC3_checks': len(ERRORS), 'max_absolute_HC3_error': max(ERRORS), 'input_sha256': INPUTS})
    (out.parent / 'LATEST.txt').write_text(out.name + '\n')
    print('SPATIAL PRIMARY', flush=True)
    print(spatial.query("branch == 'raw' and ((dataset == 'aachen_itbs' and model == 'paired_change') or (dataset == 'sudmex_tms' and model == 'baseline_motion_adjusted'))")[[ 'dataset', 'measure', 'estimate', 'ci_low', 'ci_high', 'p', 'p_holm5', 'p_maxT5']].to_string(index=False), flush=True)
    print('CROSS-DATASET HRF', flush=True)
    print(joint[['branch', *keys, 'estimate_aachen', 'estimate_sudmex', 'p_conjunction', 'p_holm_family']].to_string(index=False), flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=2):
        main()
