"""Edge cases for the exact-input FD audit; no production data are modified."""
from __future__ import annotations

from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

import audit_existing_tms_fd as audit


LIMITS = dict(mean_fd=.5, peak_fd=3., min_timepoints=100,
              translation=3., rotation=3.15, first=5, last=0, tr=2.)


def confounds(n=105, fd=.1):
    frame = pd.DataFrame(0., index=range(n), columns=audit.MOTION)
    frame["framewise_displacement"] = fd
    return frame


class WindowEvidenceTests(unittest.TestCase):
    def test_truncation_keeps_original_indices_without_double_trimming(self):
        decision = dict(status="TRUNCATE", truncate_start=49, truncate_end=299)
        self.assertEqual(audit.retained_bounds(decision, 300, LIMITS), (49, 299))
        data = confounds(300)
        data.loc[:48, "framewise_displacement"] = 99.
        data.loc[299:, "framewise_displacement"] = 99.
        got = audit.motion_metrics(data, 49, 299, "TRUNCATE", LIMITS)
        self.assertEqual(got["n_retained"], 250)
        self.assertAlmostEqual(got["retained_mean_fd"], .1)
        self.assertTrue(got["registered_pass"])

    def test_missing_bounds_fail_even_when_length_would_be_known(self):
        for bad in (None, np.nan, 12.5):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                audit.retained_bounds(dict(status="TRUNCATE", truncate_start=bad,
                                           truncate_end=300), 300, LIMITS)

    def test_first_discard_only_clamps_start_not_adds_to_start(self):
        self.assertEqual(audit.retained_bounds(dict(status="TRUNCATE", truncate_start=2,
                                                   truncate_end=270), 300, LIMITS), (5, 270))

    def test_parallel_log_associates_source_with_adjacent_trim(self):
        log = """Starting subject A
Starting subject B
  func source: D:\\branch\\sub-A_bold.nii.gz
  trim:        first 5 + last 0 volumes
  func source: D:\\branch\\sub-B_bold.nii.gz
  trim:        QC truncation to volumes 49–299
"""
        got = audit.parse_prepare_log(log)
        self.assertIn("sub-A", str(got[0]["path"]))
        self.assertTrue(audit.verify_log_rule(got[0]["rule"], "INCLUDE", 5, 300, 300, LIMITS))
        self.assertIn("sub-B", str(got[1]["path"]))
        self.assertTrue(audit.verify_log_rule(got[1]["rule"], "TRUNCATE", 49, 299, 300, LIMITS))
        self.assertFalse(audit.verify_log_rule(got[1]["rule"], "TRUNCATE", 50, 300, 300, LIMITS))

    def test_missing_qc_evidence_fails_closed(self):
        row = SimpleNamespace(dataset="sudmex_tms", subject_id="sudmex_tms_001", session="t0",
                              arm="active", prefix="dataset-sudmex_tms_sub-001_ses-t0_task-rest_run-1")
        empty = pd.DataFrame(columns=["subject", "session", "task"])
        got = audit.audit_scan(row, {"qc": empty}, {}, {})
        self.assertFalse(got["bounds_verified"])
        self.assertFalse(got["registered_pass"])
        self.assertFalse(got["strict_pass"])
        self.assertIn("Expected one QC decision", got["evidence_error"])


class MotionGateTests(unittest.TestCase):
    def test_peak_failure_is_not_hidden_by_low_mean(self):
        frame = confounds()
        frame.loc[50, "framewise_displacement"] = 3.1
        got = audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)
        self.assertLess(got["retained_mean_fd"], .3)
        self.assertFalse(got["registered_pass"])
        self.assertEqual(got["registered_fail_reasons"], "registered_peak_fd_exceeded")
        self.assertFalse(got["strict_pass"])

    def test_axis_mm_and_radians_are_handled_separately(self):
        frame = confounds()
        frame.loc[50, "trans_z"] = 3.
        frame.loc[50, "rot_x"] = np.deg2rad(3.15)
        got = audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)
        self.assertTrue(got["registered_pass"])
        self.assertAlmostEqual(got["retained_max_rotation_deg"], 3.15)
        frame.loc[50, "rot_x"] = np.deg2rad(3.16)
        self.assertFalse(audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)["registered_pass"])
        frame.loc[50, "rot_x"] = 0.
        frame.loc[50, "trans_z"] = 3.01
        self.assertFalse(audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)["registered_pass"])

    def test_legacy_truncated_reference_is_preserved(self):
        frame = confounds()
        frame["trans_y"] = 4.
        frame.loc[70, "trans_y"] = 4.4
        absolute = audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)
        relative = audit.motion_metrics(frame, 5, 105, "TRUNCATE", LIMITS)
        self.assertFalse(absolute["registered_pass"])
        self.assertTrue(relative["registered_pass"])
        self.assertAlmostEqual(relative["retained_max_translation_mm"], .4)

    def test_strict_boundary_and_high_fd_threshold_are_inclusive(self):
        frame = confounds(fd=.5)
        got = audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)
        self.assertTrue(got["registered_pass"])
        self.assertEqual(got["retained_high_fd_count"], 0)
        self.assertFalse(got["strict_pass"])
        frame = confounds(fd=0.)
        frame.loc[5:24, "framewise_displacement"] = .6
        got = audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)
        self.assertEqual(got["retained_high_fd_count"], 20)
        self.assertTrue(got["strict_pass"])
        frame.loc[25, "framewise_displacement"] = .6
        self.assertFalse(audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)["strict_pass"])
        self.assertTrue(audit.motion_metrics(confounds(fd=.3), 5, 105, "INCLUDE", LIMITS)["strict_pass"])

    def test_missing_or_invalid_retained_values_cannot_pass(self):
        frame = confounds()
        frame.loc[:4, "framewise_displacement"] = np.nan
        self.assertTrue(audit.motion_metrics(frame, 5, 105, "INCLUDE", LIMITS)["registered_pass"])
        for column, invalid in (("rot_z", np.nan), ("framewise_displacement", -.01)):
            bad = frame.copy()
            bad.loc[30, column] = invalid
            with self.subTest(column=column), self.assertRaises(ValueError):
                audit.motion_metrics(bad, 5, 105, "INCLUDE", LIMITS)
        with self.assertRaises(ValueError):
            audit.motion_metrics(frame.drop(columns="rot_z"), 5, 105, "INCLUDE", LIMITS)

    def test_minimum_length_is_not_waived_for_low_motion(self):
        self.assertTrue(audit.motion_metrics(confounds(), 5, 105, "INCLUDE", LIMITS)["registered_pass"])
        got = audit.motion_metrics(confounds(), 5, 104, "INCLUDE", LIMITS)
        self.assertFalse(got["registered_pass"])
        self.assertIn("registered_min_timepoints_not_met", got["registered_fail_reasons"])


class PairingTests(unittest.TestCase):
    def test_aachen_requires_both_intensities_and_sud_requires_t0_t1(self):
        rows = []
        for subject in ("01", "02"):
            for session in ("pre90", "post90", "pre120", "post120"):
                rows.append(dict(dataset="aachen_itbs", subject_id="aachen_itbs_" + subject,
                                 session=session, prefix=subject + session,
                                 registered_pass=not(subject == "02" and session == "post120"),
                                 strict_pass=not(session == "post120")))
        for session in ("t0", "t1", "t2", "t3"):
            rows.append(dict(dataset="sudmex_tms", subject_id="sudmex_tms_001", session=session,
                             prefix=session, registered_pass=session != "t1", strict_pass=session != "t1"))
        pairs = pd.DataFrame([
            dict(dataset="aachen_itbs", subject_id="aachen_itbs_01", baseline_session="pre90",
                 followup_session="post90", comparison="90", arm="iTBS"),
            dict(dataset="aachen_itbs", subject_id="aachen_itbs_02", baseline_session="pre90",
                 followup_session="post90", comparison="90", arm="iTBS"),
            dict(dataset="sudmex_tms", subject_id="sudmex_tms_001", baseline_session="t0",
                 followup_session="t1", comparison="t1", arm="active")])
        pairs_out, cohorts = audit.build_pair_audit(pd.DataFrame(rows), pairs)
        self.assertEqual(cohorts["registered"]["aachen_itbs"], ["aachen_itbs_01"])
        self.assertEqual(cohorts["strict"]["aachen_itbs"], [])
        self.assertEqual(cohorts["registered"]["sudmex_tms"], [])
        incomplete = pairs_out.loc[pairs_out.branch.eq("registered") & pairs_out.subject_id.eq("aachen_itbs_02")].iloc[0]
        self.assertTrue(incomplete.eligible)
        self.assertFalse(incomplete.eligible_primary_cohort)


if __name__ == "__main__":
    unittest.main()
