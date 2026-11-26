"""Integration and failure-mode checks for the TMS clinical workbook interface."""
from pathlib import Path
import sys
import tempfile
import unittest
import warnings

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pipeline/scripts/utils'))
from clinical_data import KEYS, METRICS, load_clinical_effects, merge_clinical_neural, subject_key, validate_effects


class ClinicalInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.csv = ROOT / 'data/master_clinical_effects.csv'
        cls.xlsx = ROOT / 'outputs/tms_clinical_20261125/TMS_clinical_outcomes.xlsx'
        cls.raw = pd.read_csv(cls.csv)

    def test_workbook_cached_formulas_match_csv_for_every_outcome(self):
        excel = pd.read_excel(self.xlsx, sheet_name='Clinical changes')
        validate_effects(excel)
        pd.testing.assert_frame_equal(excel[self.raw.columns], self.raw, check_dtype=False, atol=1e-10, rtol=1e-10)

    def test_zero_is_preserved_and_flags_are_measure_specific(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            frame = load_clinical_effects(self.csv).set_index('subject_id')
        self.assertEqual(len(frame), 40)
        self.assertEqual(frame.loc['sudmex_tms_002', 'craving_vas_t1'], 0)
        self.assertAlmostEqual(frame.loc['sudmex_tms_002', 'craving_vas_delta'], -0.3)
        self.assertTrue(np.isnan(frame.loc['sudmex_tms_001', 'craving_ccq_general_delta']))
        self.assertFalse(np.isnan(frame.loc['sudmex_tms_001', 'craving_vas_delta']))
        self.assertEqual(frame.impulsivity_bis_total_delta.notna().sum(), 39)

    def test_source_visits_and_missing_assessment_are_preserved(self):
        sessions = pd.read_csv(ROOT / 'data/master_clinical_sessions.csv')
        self.assertEqual(len(sessions), 154)
        self.assertEqual(sessions.subject_id.nunique(), 53)
        self.assertTrue(sessions.loc[sessions.source_visit == 'T1-4', 'session'].eq('t14').all())
        missing = sessions[sessions.impulsivity_bis_total.isna()]
        self.assertEqual(len(missing), 1)
        self.assertTrue(missing.impulsivity_bis_total_qc.eq('missing_assessment').all())
        excel = pd.read_excel(self.xlsx, sheet_name='Clinical sessions')
        pd.testing.assert_frame_equal(excel, sessions, check_dtype=False)

    def test_duplicate_identity_and_wrong_contrast_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            validate_effects(pd.concat([self.raw, self.raw.iloc[:1]], ignore_index=True))
        bad = self.raw.copy()
        bad.loc[0, 'followup_session'] = 't14'
        with self.assertRaisesRegex(ValueError, 'acute'):
            validate_effects(bad)

    def test_identifier_mismatch_and_corrupt_difference_are_rejected(self):
        self.assertEqual(subject_key('sudmex_tms', 'sub-001'), 'sudmex_tms_001')
        bad = self.raw.copy()
        bad.loc[0, 'subject_id'] = 'sudmex_tms_999'
        with self.assertRaisesRegex(ValueError, 'identifiers'):
            validate_effects(bad)
        bad = self.raw.copy()
        bad.loc[0, 'craving_vas_delta'] *= -1
        with self.assertRaisesRegex(ValueError, 't1 minus t0'):
            validate_effects(bad)

    def test_missing_score_cannot_become_zero_or_be_marked_ok(self):
        missing = self.raw.iloc[:1].copy()
        missing.loc[0, ['craving_vas_t1', 'craving_vas_delta']] = np.nan
        with self.assertRaisesRegex(ValueError, 'marked ok'):
            validate_effects(missing)
        missing.loc[0, 'craving_vas_qc'] = 't1:missing_assessment'
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp') as folder:
            file = Path(folder) / 'missing.csv'
            missing.to_csv(file, index=False)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                self.assertEqual(len(load_clinical_effects(file, metric='craving_vas')), 0)

    def test_join_uses_subject_and_contrast_not_row_order(self):
        # A key-only inventory, not synthetic neural measurements.
        inventory = self.raw.loc[self.raw.mri_pair_usable == 1, KEYS + ['arm']].iloc[::-1].copy()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            joined = merge_clinical_neural(inventory, self.xlsx, metric='craving_ccq_now')
        self.assertTrue(joined['_merge'].eq('both').all())
        self.assertEqual(joined.subject_id.tolist(), inventory.subject_id.tolist())
        first = joined.iloc[0]
        expected = self.raw.set_index('subject_id').loc[first.subject_id, 'craving_ccq_now_delta']
        self.assertEqual(first.craving_ccq_now_delta, expected)
        inventory.iloc[0, inventory.columns.get_loc('arm')] = 'sham' if inventory.iloc[0].arm == 'active' else 'active'
        with self.assertRaisesRegex(ValueError, 'arm mismatch'):
            merge_clinical_neural(inventory, self.csv, metric='craving_ccq_now')

    def test_flagged_records_require_explicit_opt_in(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            clean = load_clinical_effects(self.csv, metric='craving_ccq_general')
            reported = load_clinical_effects(self.csv, metric='craving_ccq_general', include_flagged=True)
        self.assertEqual(len(clean), 38)
        self.assertEqual(len(reported), 40)


if __name__ == '__main__':
    unittest.main(verbosity=2)
