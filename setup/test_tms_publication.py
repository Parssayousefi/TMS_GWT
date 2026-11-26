"""Numerical parity and design checks; synthetic fixtures never enter real results."""
import ast
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pipeline/scripts/analysis'))
sys.path.insert(0, str(ROOT / 'pipeline/scripts/utils'))
import tms_publication as tms
from clinical_data import METRICS as CLINICAL


class PublicationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref = tms.load_reference(ROOT)
        notebook = json.loads((ROOT / 'pipeline/Global_Workspace_Complete_Pipeline_publication.ipynb').read_text(encoding='utf-8'))
        tree = ast.parse(''.join(notebook['cells'][2]['source']))
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in {'participation_coefficient', 'classify_roles'}]
        cls.original = {'np': np, 'stats': stats, 'modules': cls.ref['modules'], 'workspace': cls.ref['workspace'], 'N_ROI': 454}
        exec(compile(ast.Module(body=functions, type_ignores=[]), 'original_metric_functions_only', 'exec'), cls.original)
        rng = np.random.default_rng(531)
        def matrix():
            x = rng.gamma(2, .02, size=(454, 454))
            x = (x + x.T) / 2
            np.fill_diagonal(x, np.nan)
            return x
        syn, red = matrix(), matrix()
        cls.matrices = {'synergy': syn, 'redundancy': red, 'synergy_total': syn*2, 'phi': syn-red, 'phi_r': syn}

    def test_exact_metric_parity_with_original_notebook(self):
        syn = tms.participation_coefficient(self.matrices['synergy'], self.ref['modules'])
        red = tms.participation_coefficient(self.matrices['redundancy'], self.ref['modules'])
        np.testing.assert_array_equal(syn, self.original['participation_coefficient'](self.matrices['synergy']))
        np.testing.assert_array_equal(red, self.original['participation_coefficient'](self.matrices['redundancy']))
        for observed, expected in zip(tms.classify_roles(syn, red, self.ref['workspace']), self.original['classify_roles'](syn, red)):
            np.testing.assert_array_equal(observed, expected)

    def test_pc_known_value_and_zero_strength_missingness(self):
        matrix = np.ones((4, 4))
        pc = tms.participation_coefficient(matrix, np.array([0, 0, 1, 1]))
        np.testing.assert_allclose(pc, 4/9)
        matrix[0] = np.nan
        self.assertTrue(np.isnan(tms.participation_coefficient(matrix, np.array([0, 0, 1, 1]))[0]))

    def test_roles_and_ties_partition_fixed_workspace(self):
        values, _ = tms.metric_values(self.matrices, self.ref)
        self.assertEqual(values['gateways_total'] + values['broadcasters_total'] + values['ties_total'], 235)
        self.assertEqual(sum(values['gateways_' + n] for n in tms.NETWORKS), values['gateways_total'])
        pc = tms.participation_coefficient(self.matrices['synergy'], self.ref['modules'])
        gw, bc, score = tms.classify_roles(pc, pc, self.ref['workspace'])
        self.assertFalse(gw.any()); self.assertFalse(bc.any())
        self.assertTrue((score[self.ref['workspace']] == 0).all())

    def test_global_mean_unique_edges_and_missing_workspace(self):
        values, _ = tms.metric_values(self.matrices, self.ref)
        upper = self.matrices['synergy'][np.triu_indices(454, 1)]
        self.assertAlmostEqual(values['global_synergy'], upper.mean())
        matrices = {k: v.copy() for k, v in self.matrices.items()}
        roi = np.flatnonzero(self.ref['workspace'])[0]
        for matrix in matrices.values():
            matrix[roi, :] = np.nan; matrix[:, roi] = np.nan
        values, _ = tms.metric_values(matrices, self.ref)
        self.assertFalse(values['full_workspace'])
        self.assertTrue(np.isnan(values['gateways_total']))
        self.assertTrue(np.isnan(values['mean_synergy_pc']))

    def test_welch_and_paired_statistics_match_scipy(self):
        a, b = np.array([3, 5, 8, 2, 6.]), np.array([-1, 2, 0, 3, -2, 1.])
        got = tms.welch_summary(a, b)
        expected = stats.ttest_ind(a, b, equal_var=False)
        self.assertAlmostEqual(got['t'], expected.statistic)
        self.assertAlmostEqual(got['p'], expected.pvalue)
        self.assertAlmostEqual(got['estimate'], a.mean()-b.mean())
        paired = tms.paired_summary(a)
        expected = stats.ttest_1samp(a, 0)
        self.assertAlmostEqual(paired['p'], expected.pvalue)
        self.assertLess(paired['ci_low'], paired['estimate'])
        self.assertGreater(paired['ci_high'], paired['estimate'])
        self.assertTrue(np.isnan(tms.welch_summary([1]*4, [0]*4)['p']))

    def test_bh_known_values(self):
        np.testing.assert_allclose(tms.bh([.01, .04, .03, .2]), [.04, .053333333333, .053333333333, .2])

    def test_aachen_pooling_never_treats_visits_as_independent(self):
        rows = []
        for sid in ['aachen_itbs_01', 'aachen_itbs_03']:
            for comp, value in [('post90_minus_pre90', 2.), ('post120_minus_pre120', 6.)]:
                r = dict(dataset='aachen_itbs', subject_id=sid, comparison=comp, arm='iTBS', pair_ready=True, same_valid_rois=True, delta_fd=.1)
                for m in tms.METRICS:
                    r.update({m+'_baseline': 3., m+'_followup': 3.+value, m+'_delta': value})
                rows.append(r)
        pairs = pd.DataFrame(rows)
        pooled = tms.pool_aachen(pairs)
        self.assertEqual(len(pooled), 2)
        np.testing.assert_allclose(pooled.gateways_FPN_delta, 4)
        pairs.loc[0, 'gateways_FPN_delta'] = np.nan
        pooled = tms.pool_aachen(pairs)
        self.assertTrue(np.isnan(pooled.loc[pooled.subject_id == 'aachen_itbs_01', 'gateways_FPN_delta'].iloc[0]))

    def test_pending_cohort_has_no_inference(self):
        inv, _ = tms.inventory(ROOT)
        empty = pd.DataFrame(columns=['dataset', 'subject_id', 'session'])
        pairs = tms.paired_changes(inv, empty)
        self.assertEqual(len(pairs[pairs.dataset == 'sudmex_tms']), 40)
        self.assertEqual(len(tms.pool_aachen(pairs)), 15)
        self.assertFalse(pairs.pair_ready.any())
        self.assertTrue(tms.neural_tests(pairs[pairs.dataset == 'sudmex_tms'], True).empty)

    def test_clinical_interface_and_adjusted_models_with_synthetic_fixture(self):
        # Fixture data are in memory / a disposable directory, never a study output.
        rng = np.random.default_rng(221)
        neural, clinical = [], []
        for i in range(20):
            sid = f'sudmex_tms_{i+1:03d}'
            active = i % 2
            row = dict(dataset='sudmex_tms', subject_id=sid, comparison='t1_minus_t0', arm='active' if active else 'sham',
                       baseline_session='t0', followup_session='t1', pair_ready=True, same_valid_rois=True, delta_fd=rng.normal(0, .02))
            for metric in tms.METRICS:
                baseline, change = rng.normal(10, 2), rng.normal(active, 1)
                row.update({metric+'_baseline': baseline, metric+'_followup': baseline+change, metric+'_delta': change})
            clinical_row = {k: row[k] for k in ['dataset', 'subject_id', 'comparison', 'arm', 'baseline_session', 'followup_session']}
            clinical_row.update(subject_num=i+1, bids_subject_id=f'sub-{i+1:03d}', analysis_cohort=1, mri_pair_usable=1)
            for metric in CLINICAL:
                baseline = rng.normal(4, .5)
                change = rng.normal(-active + .2*row['gateways_FPN_delta'], .5)
                clinical_row.update({metric+'_t0': baseline, metric+'_t1': baseline+change, metric+'_delta': change, metric+'_qc': 'ok'})
            neural.append(row); clinical.append(clinical_row)
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp', prefix='tms_clinical_fixture_') as folder:
            source = Path(folder)/'synthetic_only.csv'
            pd.DataFrame(clinical).to_csv(source, index=False)
            joined, effects, models = tms.clinical_analyses(ROOT, pd.DataFrame(neural), source, ['gateways_FPN', 'mean_synergy_pc'])
        self.assertEqual(len(joined), 20)
        self.assertEqual(len(effects), 7)
        self.assertEqual(len(models), 14)
        self.assertTrue(models.status.eq('ok').all())
        self.assertTrue(models['n'].eq(20).all())
        self.assertTrue(np.isfinite(models.estimate).all())
        self.assertTrue(models.q_family.between(0, 1).all())
        adjusted = tms.motion_sensitivity(pd.DataFrame(neural))
        self.assertEqual(len(adjusted), 25)
        self.assertTrue(adjusted.status.eq('ok').all())

    def test_new_notebook_has_no_drug_analysis_code(self):
        notebook = json.loads((ROOT/'output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb').read_text(encoding='utf-8'))
        code = '\n'.join(''.join(c['source']) for c in notebook['cells'] if c['cell_type'] == 'code')
        for forbidden in ['basel_lam', 'basel_lpm', 'london_lsd', 'master_subjective_effects', 'publication-notebook-2026-09-26']:
            self.assertNotIn(forbidden, code)
        for i, cell in enumerate(notebook['cells']):
            if cell['cell_type'] == 'code':
                compile(''.join(cell['source']), f'cell{i}', 'exec')

    def test_unestimable_effect_does_not_get_a_direction(self):
        table = pd.DataFrame([{'metric': m, 'estimate': np.nan, 'q_family': np.nan}
                              for m in ['gateways_FPN', 'gateways_DAN', 'mean_redundancy_pc']])
        comparison = tms.comparison_to_paper(table, pd.DataFrame())
        self.assertTrue(comparison.SUDMEX.eq('Not estimable').all())
        table['estimate'] = 0.
        table['q_family'] = 1.
        comparison = tms.comparison_to_paper(table, pd.DataFrame())
        self.assertTrue(comparison.SUDMEX.str.startswith('zero estimated change').all())

    def test_regional_changes_preserve_pair_identity_and_coverage(self):
        regions, pairs = [], []
        for subject, same in [('sudmex_tms_001', True), ('sudmex_tms_002', False)]:
            pairs.append(dict(dataset='sudmex_tms', subject_id=subject, comparison='t1_minus_t0', arm='active',
                              baseline_session='t0', followup_session='t1', pair_ready=True, same_valid_rois=same))
            for session, offset in [('t0', 0), ('t1', 2)]:
                for roi in [3, 17, 100]:
                    regions.append(dict(dataset='sudmex_tms', subject_id=subject, session=session, roi=roi, network='FPN',
                                        synergy_pc=roi/100 + offset, redundancy_pc=roi/200 + offset, rank_score=roi + offset))
        result = tms.regional_changes(pd.DataFrame(regions), pd.DataFrame(pairs))
        self.assertEqual(result.subject_id.unique().tolist(), ['sudmex_tms_001'])
        self.assertEqual(result.roi.tolist(), [3, 17, 100])
        np.testing.assert_allclose(result[['synergy_pc', 'redundancy_pc', 'rank_score']], 2)

    def test_forest_plot_layout_with_explicit_synthetic_data(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        notebook = json.loads((ROOT/'output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb').read_text(encoding='utf-8'))
        setup = next(c for c in notebook['cells'] if c['cell_type'] == 'code')
        tree = ast.parse(''.join(setup['source']))
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'forest')
        namespace = {'np': np, 'plt': plt, 'tms': tms}
        exec(compile(ast.Module(body=[function], type_ignores=[]), 'forest_only', 'exec'), namespace)
        metrics = ['gateways_' + n for n in tms.NETWORKS]
        fake = pd.DataFrame([{'metric': m, 'estimate': i/4-1, 'ci_low': i/4-1.7, 'ci_high': i/4-.3} for i, m in enumerate(metrics)])
        figure = namespace['forest'](fake, metrics, 'SYNTHETIC FIXTURE — layout check only', 'Synthetic difference (regions; 95% CI)')
        output = ROOT / 'tmp/jupyter-notebook/SYNTHETIC_forest_layout.png'
        figure.savefig(output, dpi=120, bbox_inches='tight')
        self.assertTrue(output.is_file())
        plt.close(figure)


if __name__ == '__main__':
    unittest.main(verbosity=2)
