"""Verify completed TMS notebook tables independently of the analysis helper.

Uses SciPy/statsmodels for t tests, intervals and FDR, and a direct NumPy HC3
calculation for the regression tables. Does not recalculate neural measures,
select hypotheses, edit source data, or import tms_publication.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.weightstats import CompareMeans, DescrStatsW


def same(actual, expected, name):
    np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-10,
                               equal_nan=True, err_msg=name)


def fdr(table):
    for family, block in table.groupby('family'):
        finite = block.p.notna()
        expected = multipletests(block.loc[finite, 'p'], method='fdr_bh')[1]
        same(block.loc[finite, 'q_family'], expected, family + ' FDR')


def check_welch(table, frame):
    for row in table.itertuples():
        assert row.status == 'ok', (row.metric, row.status)
        a = frame.loc[frame.arm == 'active', row.metric + '_delta'].dropna()
        b = frame.loc[frame.arm == 'sham', row.metric + '_delta'].dropna()
        assert (len(a), len(b)) == (row.n_active, row.n_sham)
        test = stats.ttest_ind(a, b, equal_var=False)
        comparison = CompareMeans(DescrStatsW(a), DescrStatsW(b))
        interval = comparison.tconfint_diff(usevar='unequal')
        same([row.estimate, row.mean_change_active, row.mean_change_sham,
              row.t, row.p, row.ci_low, row.ci_high],
             [a.mean() - b.mean(), a.mean(), b.mean(), test.statistic,
              test.pvalue, *interval], row.metric + ' Welch')
        ndf = len(a) + len(b) - 2
        pooled_sd = np.sqrt(((len(a)-1)*a.var(ddof=1) +
                             (len(b)-1)*b.var(ddof=1)) / ndf)
        g = (1-3/(4*ndf-1)) * (a.mean()-b.mean()) / pooled_sd
        same(row.effect_size, g, row.metric + ' Hedges g')
    fdr(table)


def check_aachen(table, pooled, pairs):
    original = pairs.loc[pairs.dataset == 'aachen_itbs']
    assert len(original) == 30 and len(pooled) == 15
    assert not pooled.subject_id.duplicated().any()
    assert original.groupby('subject_id').size().eq(2).all()
    for row in table.itertuples():
        metric = row.metric + '_delta'
        reconstructed = original.groupby('subject_id')[metric].mean()
        same(pooled.set_index('subject_id')[metric].sort_index(),
             reconstructed.sort_index(), row.metric + ' within-person pooling')
        values = pooled[metric].dropna()
        assert row.n == len(values) == 15 and row.status == 'ok'
        test = stats.ttest_1samp(values, 0)
        interval = DescrStatsW(values).tconfint_mean()
        same([row.estimate, row.t, row.p, row.ci_low, row.ci_high],
             [values.mean(), test.statistic, test.pvalue, *interval],
             row.metric + ' one-sample t')
    fdr(table)


def direct_hc3(frame, outcome, predictor, covariates, standardize):
    columns = list(dict.fromkeys([outcome, predictor, *covariates]))
    complete = frame[columns].replace([np.inf, -np.inf], np.nan).dropna()
    kept = [c for c in covariates if complete[c].nunique() >= 2]
    design = complete[[predictor, *kept]].astype(float).to_numpy(copy=True)
    if standardize:
        design[:, 0] = (design[:, 0] - design[:, 0].mean()) / design[:, 0].std(ddof=1)
    design = np.column_stack([np.ones(len(design)), design])
    y = complete[outcome].to_numpy(float)
    assert np.linalg.matrix_rank(design) == design.shape[1]
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    bread = np.linalg.inv(design.T @ design)
    hat = np.sum((design @ bread) * design, axis=1)
    adjusted_residual = (y - design @ beta) / (1-hat)
    meat = design.T @ (design * adjusted_residual[:, None]**2)
    covariance = bread @ meat @ bread
    se = np.sqrt(covariance[1, 1])
    df = len(y)-design.shape[1]
    half = stats.t.ppf(.975, df) * se
    p = 2 * stats.t.sf(abs(beta[1]/se), df)
    return len(y), [beta[1], beta[1]-half, beta[1]+half, p]


def check_regressions(table, frame, clinical=False):
    frame = frame.copy()
    frame['active'] = frame.arm.eq('active').astype(int)
    for row in table.itertuples():
        assert row.status == 'ok', row
        if clinical:
            outcome, predictor = row.clinical_metric + '_delta', row.neural_metric + '_delta'
            covariates = ['active', row.clinical_metric + '_t0',
                          row.neural_metric + '_baseline', 'delta_fd']
        else:
            outcome, predictor = row.metric + '_delta', 'active'
            covariates = [row.metric + '_baseline', 'delta_fd']
        n, estimates = direct_hc3(frame, outcome, predictor, covariates, clinical)
        assert n == row.n
        same([row.estimate, row.ci_low, row.ci_high, row.p], estimates,
             outcome + ' independent HC3')
    fdr(table)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--aachen-reference', type=Path)
    args = parser.parse_args()
    folder = args.run_directory.resolve()
    root = Path(__file__).resolve().parents[1]
    if root not in folder.parents:
        raise ValueError('Verification output must be inside TMS_GWT')
    status = json.loads((folder / 'run_status.json').read_text())
    assert all(status['cohort_readiness'].values())
    assert status['validated_scans'] == 193
    sources = {}

    def read(name):
        path = folder / (name + '.csv')
        sources[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        return pd.read_csv(path)

    pairs = read('participant_paired_changes')
    assert len(pairs) == 70 and pairs.pair_ready.all() and pairs.same_valid_rois.all()
    assert not pairs.duplicated(['dataset', 'subject_id', 'comparison']).any()
    sudmex = pairs.loc[pairs.dataset == 'sudmex_tms'].copy()
    assert len(sudmex) == 40 and not sudmex.subject_id.duplicated().any()
    assert sudmex.arm.value_counts().to_dict() == {'active': 22, 'sham': 18}
    assert sudmex.baseline_session.eq('t0').all() and sudmex.followup_session.eq('t1').all()
    coverage = read('scan_coverage')
    assert len(coverage) == 193 and coverage.n_valid_rois.eq(454).all()
    assert coverage.n_valid_workspace.eq(235).all()

    neural = read('sudmex_neural_results')
    assert len(neural) == 25
    check_welch(neural, sudmex)
    primary = read('primary_FPN_result')
    pd.testing.assert_frame_equal(primary.reset_index(drop=True),
                                  neural.loc[neural.metric == 'gateways_FPN'].reset_index(drop=True))
    aachen = read('aachen_neural_results')
    check_aachen(aachen, read('aachen_participant_pooled_changes'), pairs)
    if args.aachen_reference:
        pd.testing.assert_frame_equal(aachen, pd.read_csv(args.aachen_reference / 'aachen_neural_results.csv'))

    joined = read('clinical_neural_participant_changes')
    assert len(joined) == 40 and not joined.subject_id.duplicated().any()
    pd.testing.assert_series_equal(joined.set_index('subject_id').arm.sort_index(),
                                    sudmex.set_index('subject_id').arm.sort_index())
    clinical = read('clinical_treatment_changes')
    assert len(clinical) == 7
    check_welch(clinical, joined)
    expected_n = {'craving_vas': 40, 'craving_ccq_general': 38, 'craving_ccq_now': 40,
                  'impulsivity_bis_total': 39, 'impulsivity_bis_cognitive': 39,
                  'impulsivity_bis_motor': 39, 'impulsivity_bis_nonplanning': 39}
    for row in clinical.itertuples():
        assert row.n_active + row.n_sham == expected_n[row.metric]
    motion = read('motion_baseline_adjusted_neural_effects')
    assert len(motion) == 25
    check_regressions(motion, sudmex)
    associations = read('clinical_neural_associations')
    assert len(associations) == 35
    check_regressions(associations, joined, clinical=True)

    notebook_path = folder / 'TMS_Global_Workspace_Publication.ipynb'
    notebook = nbformat.read(notebook_path, as_version=4)
    cells = [c for c in notebook.cells if c.cell_type == 'code']
    assert len(cells) == 15 and all(c.execution_count is not None for c in cells)
    assert not any(o.output_type == 'error' for c in cells for o in c.outputs)
    if args.aachen_reference:
        previous = nbformat.read(args.aachen_reference / notebook_path.name, as_version=4)
        assert [c.source for c in cells] == [c.source for c in previous.cells if c.cell_type == 'code']
        prior_provenance = json.loads((args.aachen_reference / 'provenance.json').read_text())
        prior_hashes = {item['path'].replace('\\', '/'): item['sha256']
                        for item in prior_provenance['inputs']}
        for name in ['pipeline/scripts/analysis/tms_publication.py',
                     'pipeline/scripts/utils/clinical_data.py',
                     'outputs/tms_clinical_20261125/TMS_clinical_outcomes.xlsx']:
            assert hashlib.sha256((root / name).read_bytes()).hexdigest() == prior_hashes[name]
    receipt = {'verified_utc': datetime.now(timezone.utc).isoformat(), 'status': 'passed',
               'validated_scans': 193, 'sudmex_participants': 40, 'aachen_participants': 15,
               'sudmex_neural_tests': 25, 'aachen_neural_tests': 25,
               'clinical_treatment_tests': 7, 'motion_models': 25, 'clinical_association_models': 35,
               'notebook_executed_cells': len(cells), 'aachen_unchanged': bool(args.aachen_reference),
               'analysis_code_and_clinical_input_unchanged_since_aachen': bool(args.aachen_reference),
               'checks': 'SciPy t tests, statsmodels confidence intervals/BH, NumPy HC3 regression, joins/coverage',
               'verified_table_sha256': sources,
               'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (folder / 'independent_statistics_check.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if not key.endswith('sha256')}, indent=2))


if __name__ == '__main__':
    main()
