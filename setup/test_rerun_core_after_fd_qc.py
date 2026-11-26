"""Bounded correctness checks for cohort filtering and FD reanalysis guards."""
import json
import unittest
from unittest import mock

import numpy as np
import pandas as pd

import rerun_core_after_fd_qc as core


class ProvenanceAndConjunctionTests(unittest.TestCase):
    def test_conjunction_requires_estimable_results_in_both_orders(self):
        self.assertTrue(np.isnan(core.conjunction_p(.001, np.nan)))
        self.assertTrue(np.isnan(core.conjunction_p(np.nan, .001)))
        self.assertTrue(np.isnan(core.conjunction_p(np.nan, np.nan)))
        self.assertEqual(core.conjunction_p(.01, .2), .2)
        self.assertEqual(core.conjunction_p(.2, .01), .2)
        adjusted = core.ws.holm([.02, np.nan, .4])
        np.testing.assert_allclose(adjusted, [.06, np.nan, .8], equal_nan=True)

    def test_unreleased_audit_stops_before_data_read(self):
        for receipt in ({}, {"state": "complete", "all_193_windows_verified": False,
                              "changed_source_files": []},
                        {"state": "complete", "all_193_windows_verified": True,
                         "changed_source_files": ["changed.tsv"]}):
            with self.subTest(receipt=receipt), mock.patch.object(core.Path, "read_text", return_value=json.dumps(receipt)):
                with self.assertRaisesRegex(ValueError, "has not released"):
                    core.load_verified_audit(core.ROOT / "output/test/scan_audit.csv")

    def test_modified_whitelist_rejected_even_with_complete_receipt(self):
        receipt = {"state": "complete", "all_193_windows_verified": True, "changed_source_files": [],
                   "output_sha256": {name: "original" for name in ("scan_audit.csv", "pair_audit.csv",
                       "cohort_subjects.json", "baseline_pool_audit.csv", "baseline_failures.csv")}}
        with mock.patch.object(core.Path, "read_text", return_value=json.dumps(receipt)), mock.patch.object(core, "sha", return_value="modified"):
            with self.assertRaisesRegex(ValueError, "output changed"):
                core.load_verified_audit(core.ROOT / "output/test/scan_audit.csv")

    def test_fresh_normalization_rejects_wrong_policy_and_partial_release(self):
        path = core.ROOT / "output/test/registered/normalized_scan_metrics.csv"
        with self.assertRaisesRegex(ValueError, "identify its policy"):
            core.verify_recalibrated_input(path, "strict", pd.DataFrame(), path)
        with mock.patch.object(core.Path, "read_text", return_value=json.dumps({"passed": True, "state": "prepared"})):
            with self.assertRaisesRegex(ValueError, "not completed"):
                core.verify_recalibrated_input(path, "registered", pd.DataFrame(), path)

    def test_fresh_normalization_rejects_changed_scalar_hash(self):
        path = core.ROOT / "output/test/registered/normalized_scan_metrics.csv"
        def contents(file, **kwargs):
            if file.name == "output_sha256.json":
                return json.dumps({str(core.Path("registered") / path.name): "expected"})
            return json.dumps({"passed": True, "state": "complete", "sources_unchanged": True,
                               "qc_policy": "registered"})
        with mock.patch.object(core.Path, "read_text", contents), mock.patch.object(core, "sha", return_value="changed"):
            with self.assertRaisesRegex(ValueError, "differs from its release hash"):
                core.verify_recalibrated_input(path, "registered", pd.DataFrame(), path)


class CohortAndAggregationTests(unittest.TestCase):
    def test_complete_intensities_and_duplicate_scan_guard(self):
        core.check_pair_filter()
        scans = pd.DataFrame([dict(dataset="aachen_itbs", subject_id="A", session="pre90")] * 2)
        plan = pd.DataFrame(columns=["dataset", "subject_id", "comparison"])
        with self.assertRaisesRegex(ValueError, "Duplicate scan"):
            core.retain_plan(scans, plan)

    def test_all_post_includes_t14_and_t4_equal_within_person_weight(self):
        rows = []
        for session, value, fd in [("t0", 10., .1), ("t1", 12., .2), ("t14", 20., .3),
                                   ("t2", 14., .4), ("t3", 16., .5), ("t4", 18., .6)]:
            rows.append({"dataset": "sudmex_tms", "subject_id": "sudmex_tms_001", "session": session,
                         "arm": "active", "valid_roi_signature": "fixed", "mean_fd": fd,
                         **{metric: value for metric in core.ws.METRICS}})
        scans = pd.DataFrame(rows)
        got = core.all_post(scans).iloc[0]
        self.assertEqual(got.n_post, 5)
        self.assertEqual(got.post_sessions, "t1|t14|t2|t3|t4")
        self.assertAlmostEqual(got.delta_fd, .3)
        self.assertEqual(got[core.ws.METRICS[0] + "_delta"], 6.)
        self.assertTrue(core.all_post(scans[scans.session.ne("t0")]).empty)
        # A passing later visit is still valid for this descriptive average
        # when an acute t1 has failed QC; it cannot enter the randomized pair.
        later = core.all_post(scans[scans.session.ne("t1")]).iloc[0]
        self.assertEqual(later.n_post, 4)
        self.assertEqual(later[core.ws.METRICS[0] + "_followup"], 17.)

    def test_all20_influence_inventory_and_known_omission_mean(self):
        rng = np.random.default_rng(82)
        pooled = pd.DataFrame(rng.normal(size=(6, 20)), index=list("ABCDEF"),
                              columns=[metric + "_delta" for metric in core.ws.INFO_METRICS])
        pooled.iloc[:, 0] = np.arange(1., 7.)
        result = core.information_leave_one_out(pooled, "raw")
        self.assertEqual(len(result), 120)
        self.assertTrue(result.groupby("omitted_subject").metric.nunique().eq(20).all())
        self.assertTrue(result.n.eq(5).all())
        self.assertTrue(result.diagnostic_only.all())
        selected = result[result.omitted_subject.eq("A") & result.metric.eq(core.ws.INFO_METRICS[0])].iloc[0]
        self.assertAlmostEqual(selected.estimate, 4.)
        self.assertTrue((result.p_maxT_all >= result.p_maxT_family).all())
        self.assertTrue((result.p_maxT_family >= result.p_sign_flip).all())


class GroupMotionPolicyTests(unittest.TestCase):
    def test_continuous_motion_predictors_rejected_before_fit(self):
        frame = pd.DataFrame({"y": range(10), "intercept": 1., "delta_fd": np.arange(10.)})
        for name in ("delta_fd", "mean_fd", "framewise_displacement", "motion", "rot_x", "trans_y"):
            with self.subTest(name=name), mock.patch.object(core.ws, "hc3") as fit:
                with self.assertRaisesRegex(ValueError, "excluded from group models"):
                    core.group_hc3(frame, "y", ["intercept", name], "intercept")
                fit.assert_not_called()

    def test_baseline_models_do_not_use_even_missing_motion_metadata(self):
        rng = np.random.default_rng(932)
        rows = []
        for subject in range(24):
            active = subject % 2
            row = {"dataset": "sudmex_tms", "subject_id": f"sudmex_tms_{subject:03d}",
                   "arm": "active" if active else "sham", "delta_fd": rng.normal()}
            for metric in core.ws.METRICS:
                baseline = rng.normal()
                followup = .6 * baseline + .4 * active + rng.normal()
                row.update({metric + "_baseline": baseline, metric + "_followup": followup,
                            metric + "_delta": followup-baseline})
            rows.append(row)
        pairs = pd.DataFrame(rows)
        original = core.sudmex_tests(pairs)
        perturbed = core.sudmex_tests(pairs.assign(delta_fd=np.nan))
        self.assertEqual(set(original.model), {"unadjusted_change_difference", "baseline_adjusted"})
        np.testing.assert_allclose(original[["estimate", "p", "p_holm_all40"]],
                                   perturbed[["estimate", "p", "p_holm_all40"]], equal_nan=True)
        metric = core.ws.METRICS[0]
        design = pairs.assign(intercept=1., active=pairs.arm.eq("active").astype(float))
        expected = core.gwt.direct_hc3(design, metric + "_followup", ["intercept", "active", metric + "_baseline"], "active")
        selected = original[original.metric.eq(metric) & original.model.eq("baseline_adjusted")].iloc[0]
        self.assertAlmostEqual(selected.estimate, expected["estimate"])
        self.assertAlmostEqual(selected.p, expected["p"])
        spatial = core.spatial_sudmex(pairs.assign(delta_fd=np.nan), "raw")
        self.assertEqual(set(spatial.model), {"baseline_adjusted"})
        self.assertFalse(spatial.predictors.str.contains("fd|motion", regex=True).any())

    def test_independent_numeric_self_check_has_no_motion_design(self):
        core.self_check_no_group_motion()
        for _, predictors, _ in core.GROUP_MODEL_DESIGNS:
            self.assertFalse(any("fd" in name or "motion" in name for name in predictors))


if __name__ == "__main__":
    unittest.main()
