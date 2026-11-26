"""Synthetic failures and boundary cases for the dataset's declared motion gate."""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pandas as pd

import qc_ds001927_derivatives as qc


POLICY = {
    "dataset": "ds001927", "tr_seconds": 1.25, "expected_volumes": 600,
    "sessions": ["ctr", "fro", "occ"], "tasks": ["restpre", "restpost"],
    "initial_volumes_to_remove": 8, "minimum_remaining_seconds": 480.0,
    "max_translation_mm": 3.0, "max_rotation_degrees": 3.15,
    "max_mean_fd_mm": .30, "high_fd_threshold_mm": .50,
    "max_high_fd_fraction": .20, "strict_sensitivity_max_mean_fd_mm": .20,
    "standardized_dvars_flag": 1.5, "auto_salvage": False,
    "delete_and_rejoin_scattered_frames": False,
}


def confounds(n=600, fd=.1):
    frame = pd.DataFrame({
        "trans_x": np.full(n, 10.), "trans_y": np.full(n, -7.),
        "trans_z": np.full(n, 6.), "rot_x": np.full(n, .2),
        "rot_y": np.full(n, -.1), "rot_z": np.full(n, .4),
        "framewise_displacement": np.full(n, fd), "std_dvars": np.ones(n),
    })
    if n:
        frame.loc[0, ["framewise_displacement", "std_dvars"]] = np.nan
    return frame


def leading_flags(frame, n):
    frame["non_steady_state_outlier00"] = (np.arange(len(frame)) < n).astype(int)


class MotionGateTests(unittest.TestCase):
    def test_mm_and_radians_relative_to_first_retained_not_origin_or_euclidean_norm(self):
        frame = confounds()
        frame.loc[:7, "trans_x"] = 10000
        frame.loc[599, "trans_x"] = 13
        frame.loc[599, "trans_y"] = -4
        frame.loc[599, "rot_z"] += np.deg2rad(3)
        metrics, flags = qc.assess_confounds(frame, POLICY)
        self.assertTrue(metrics["motion_gate_pass"])
        self.assertAlmostEqual(metrics["max_axis_translation_mm"], 3)
        self.assertAlmostEqual(metrics["max_axis_rotation_degrees"], 3)
        self.assertTrue(np.isnan(flags.loc[7, "max_axis_translation_relative_mm"]))
        self.assertEqual(flags.loc[8, "max_axis_translation_relative_mm"], 0)

    def test_exact_inclusive_boundaries_and_stricter_sensitivity(self):
        frame = confounds(fd=0)
        leading_flags(frame, 100)  # 500 retained; exactly 100/500 high-FD frames.
        frame.loc[100:199, "framewise_displacement"] = 1.5  # mean FD=.30.
        frame.loc[599, "trans_x"] += 3
        frame.loc[599, "rot_y"] += np.deg2rad(3.15)
        metrics, _ = qc.assess_confounds(frame, POLICY)
        self.assertTrue(metrics["motion_gate_pass"])
        self.assertAlmostEqual(metrics["fraction_fd_above_threshold_retained"], .2)
        self.assertAlmostEqual(metrics["mean_fd_mm"], .3)
        self.assertFalse(metrics["strict_mean_fd_sensitivity_pass"])
        self.assertFalse(metrics["strict_sensitivity_motion_gate_pass"])

    def test_each_excess_motion_threshold_can_fail_without_salvage(self):
        for column, amount, failed_check in [
            ("trans_z", 3.0001, "translation"),
            ("rot_x", np.deg2rad(3.1501), "rotation"),
        ]:
            with self.subTest(column=column):
                frame = confounds()
                frame.loc[599, column] += amount
                metrics, flags = qc.assess_confounds(frame, POLICY)
                self.assertFalse(metrics["motion_gate_pass"])
                self.assertFalse(metrics["motion_checks"][failed_check])
                self.assertEqual(metrics["n_retained_volumes"], 592)
                self.assertTrue(flags.loc[599, "retained_for_motion_assessment"])
                self.assertEqual(metrics["interior_frames_removed"], 0)

    def test_mean_fd_is_a_gate_independent_of_high_fd_fraction(self):
        metrics, _ = qc.assess_confounds(confounds(fd=.3001), POLICY)
        self.assertFalse(metrics["motion_gate_pass"])
        self.assertFalse(metrics["motion_checks"]["mean_fd"])
        self.assertTrue(metrics["motion_checks"]["high_fd_fraction"])

    def test_fraction_above_point_five_is_strict_and_uses_retained_denominator(self):
        frame = confounds(fd=0)
        leading_flags(frame, 100)
        frame.loc[100:200, "framewise_displacement"] = .5001
        frame.loc[201, "framewise_displacement"] = .5
        metrics, flags = qc.assess_confounds(frame, POLICY)
        self.assertEqual(metrics["fd_above_threshold_count_retained"], 101)
        self.assertAlmostEqual(metrics["fraction_fd_above_threshold_retained"], 101/500)
        self.assertFalse(metrics["motion_gate_pass"])
        self.assertTrue(metrics["motion_checks"]["mean_fd"])
        self.assertFalse(flags.loc[201, "fd_above_declared_threshold"])

    def test_leading_nonsteady_prefix_extends_trim_and_preserves_original_time(self):
        frame = confounds()
        # Standard one-hot fMRIPrep columns, contiguous only after their union.
        for index in range(12):
            frame[f"non_steady_state_outlier{index:02d}"] = (np.arange(600) == index).astype(int)
        metrics, flags = qc.assess_confounds(frame, POLICY)
        retained = flags.loc[flags.retained_for_motion_assessment]
        self.assertTrue(metrics["motion_gate_pass"])
        self.assertEqual(metrics["trim_first_volumes"], 12)
        self.assertEqual(retained.iloc[0].original_volume_index, 12)
        self.assertEqual(retained.iloc[0].original_volume_number_1based, 13)
        self.assertEqual(retained.iloc[0].time_seconds, 15)
        self.assertEqual(retained.iloc[-1].original_volume_index, 599)
        self.assertEqual(len(flags), 600)

    def test_nonprefix_flags_require_review_even_inside_default_trim(self):
        for flagged_index in (5, 100):
            with self.subTest(flagged_index=flagged_index):
                frame = confounds()
                frame["non_steady_state_outlier00"] = (np.arange(600) == flagged_index).astype(int)
                metrics, flags = qc.assess_confounds(frame, POLICY)
                self.assertFalse(metrics["motion_gate_pass"])
                self.assertTrue(metrics["review_reasons"])
                self.assertEqual(metrics["trim_first_volumes"], 8)
                self.assertEqual(metrics["n_retained_volumes"], 592)
                self.assertEqual(flags.original_volume_index.tolist(), list(range(600)))

    def test_duration_boundary_does_not_round_short_records_up(self):
        for prefix, expected in ((216, True), (217, False)):
            with self.subTest(prefix=prefix):
                frame = confounds()
                leading_flags(frame, prefix)
                metrics, _ = qc.assess_confounds(frame, POLICY)
                self.assertEqual(metrics["motion_gate_pass"], expected)
                self.assertAlmostEqual(metrics["retained_duration_seconds"], (600-prefix)*1.25)

    def test_missing_columns_and_retained_nonfinite_values_never_become_zero(self):
        frame = confounds().drop(columns="rot_y")
        metrics, _ = qc.assess_confounds(frame, POLICY)
        self.assertFalse(metrics["motion_gate_pass"])
        self.assertIn("missing_required_column:rot_y", metrics["input_errors"])
        self.assertIsNone(metrics["max_axis_rotation_degrees"])
        for missing_index, expected in ((7, True), (8, False), (599, False)):
            with self.subTest(missing_index=missing_index):
                frame = confounds()
                frame.loc[missing_index, "framewise_displacement"] = np.nan
                metrics, _ = qc.assess_confounds(frame, POLICY)
                self.assertEqual(metrics["motion_gate_pass"], expected)

    def test_negative_fd_and_invalid_nonsteady_flags_block_inclusion(self):
        frame = confounds()
        frame.loc[300, "framewise_displacement"] = -.1
        metrics, _ = qc.assess_confounds(frame, POLICY)
        self.assertFalse(metrics["motion_gate_pass"])
        self.assertIn("negative_fd_after_trim", metrics["input_errors"])
        frame = confounds()
        frame["non_steady_state_outlier00"] = 0.0
        frame.loc[0, "non_steady_state_outlier00"] = np.nan
        metrics, _ = qc.assess_confounds(frame, POLICY)
        self.assertFalse(metrics["motion_gate_pass"])
        self.assertTrue(any("invalid_nonsteady" in e for e in metrics["input_errors"]))

    def test_standardized_dvars_is_diagnostic_not_a_primary_exclusion(self):
        frame = confounds()
        frame.loc[8:, "std_dvars"] = 20
        frame.loc[8, "std_dvars"] = 1.5
        metrics, flags = qc.assess_confounds(frame, POLICY)
        self.assertTrue(metrics["motion_gate_pass"])
        self.assertEqual(metrics["std_dvars_above_threshold_count_retained"], 591)
        self.assertFalse(flags.loc[8, "std_dvars_above_declared_threshold"])
        metrics, _ = qc.assess_confounds(frame.drop(columns="std_dvars"), POLICY)
        self.assertTrue(metrics["motion_gate_pass"])
        self.assertIn("std_dvars_missing_diagnostic_only", metrics["warnings"])

    def test_partial_confounds_are_not_a_shorter_valid_run(self):
        metrics, flags = qc.assess_confounds(confounds(n=590), POLICY)
        self.assertFalse(metrics["motion_gate_pass"])
        self.assertIn("confounds_rows_590_expected_600", metrics["input_errors"])
        self.assertEqual(flags.iloc[-1].original_volume_index, 589)

    def test_policy_requires_explicit_valid_values(self):
        cases = []
        policy = copy.deepcopy(POLICY); policy.pop("max_mean_fd_mm"); cases.append(policy)
        for value in (None, True, ".3", float("nan"), -1):
            policy = copy.deepcopy(POLICY); policy["max_mean_fd_mm"] = value; cases.append(policy)
        policy = copy.deepcopy(POLICY); policy["auto_salvage"] = True; cases.append(policy)
        for index, policy in enumerate(cases):
            with self.subTest(index=index), self.assertRaises(ValueError):
                qc.validate_policy(policy)


class FileCompletionTests(unittest.TestCase):
    def setUp(self):
        test_root = qc.ROOT / "tmp/ds1927_qc_tests"
        test_root.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix="case_", dir=test_root)
        # Verify the exact recursive cleanup target stays inside the test area.
        Path(self.temp.name).resolve().relative_to(test_root.resolve())
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.derivatives = self.root / "derivatives"
        self.out = self.root / "qc"

    def image(self, path, shape=(3, 4, 5, 600), shift=0, tr=1.25):
        affine = np.diag([2., 2., 2., 1.]); affine[0, 3] = -96 + shift
        image = nib.Nifti1Image(np.zeros(shape, dtype=np.uint8), affine)
        image.header.set_xyzt_units("mm", "sec")
        image.header.set_zooms((2, 2, 2, tr) if len(shape) == 4 else (2, 2, 2))
        nib.save(image, path)

    def run_files(self, session="ctr", task="restpre", volumes=600):
        directory = self.derivatives / "sub-s001" / f"ses-{session}" / "func"
        directory.mkdir(parents=True, exist_ok=True)
        prefix = f"sub-s001_ses-{session}_task-{task}"
        spatial = f"{prefix}_space-{qc.STANDARD_SPACE}_res-2"
        bold = directory / f"{spatial}_desc-preproc_bold.nii.gz"
        mask = directory / f"{spatial}_desc-brain_mask.nii.gz"
        tsv = directory / f"{prefix}_desc-confounds_timeseries.tsv"
        self.image(bold, shape=(3, 4, 5, volumes))
        self.image(mask, shape=(3, 4, 5))
        confounds().to_csv(tsv, sep="\t", index=False, na_rep="n/a")
        return bold, mask, tsv

    def complete_fixture(self):
        return {(session, task): self.run_files(session, task)
                for session in POLICY["sessions"] for task in POLICY["tasks"]}

    def test_missing_pilot_outputs_write_all_six_pending_records(self):
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertEqual(len(result["runs"]), 6)
        self.assertEqual(result["files_complete_run_count"], 0)
        self.assertFalse(result["six_run_files_complete"])
        self.assertFalse(result["analysis_release"])
        self.assertTrue((self.out / "derivative_qc.json").is_file())
        self.assertEqual(len(pd.read_csv(self.out / "run_qc.csv")), 6)
        self.assertTrue(all(len(r["file_errors"]) == 3 for r in result["runs"]))

    def test_complete_headers_and_motion_do_not_release_anatomy_or_read_payloads(self):
        self.complete_fixture()
        policy_path = self.root / "policy.json"
        blob = json.dumps(POLICY).encode(); policy_path.write_bytes(blob)
        policy, digest = qc.load_policy(policy_path)
        with patch.object(nib.arrayproxy.ArrayProxy, "__array__", side_effect=AssertionError("Payload read")):
            result = qc.audit_subject(self.derivatives, "sub-s001", policy, self.out, policy_path, digest)
        self.assertTrue(result["six_run_files_complete"])
        self.assertTrue(result["six_run_motion_gate_pass"])
        self.assertEqual(result["motion_pass_run_count"], 6)
        self.assertFalse(result["analysis_release"])
        self.assertFalse(result["fmriprep_process_exit_verified"])
        self.assertEqual(result["policy_sha256"], hashlib.sha256(blob).hexdigest())
        for record in result["runs"]:
            flags = pd.read_csv(record["frame_flags"])
            self.assertEqual(len(flags), 600)
            self.assertEqual(flags.original_volume_index.iloc[-1], 599)
            self.assertEqual(flags.retained_for_motion_assessment.sum(), 592)

    def test_partial_cohort_does_not_silently_skip_five_runs(self):
        self.run_files()
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertEqual(result["files_complete_run_count"], 1)
        self.assertEqual(result["motion_pass_run_count"], 1)
        self.assertFalse(result["six_run_motion_gate_pass"])
        self.assertEqual(len(result["runs"]), 6)

    def test_duplicate_run_image_blocks_ambiguous_selection(self):
        runs = self.complete_fixture()
        bold = runs[("ctr", "restpre")][0]
        duplicate = bold.with_name(bold.name.replace("_space-", "_run-2_space-"))
        duplicate.write_bytes(bold.read_bytes())
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertFalse(result["six_run_files_complete"])
        self.assertIn("bold_count_2_expected_1", result["runs"][0]["file_errors"])

    def test_unexpected_seventh_run_blocks_exact_design_gate(self):
        self.complete_fixture(); self.run_files(task="restother")
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertEqual(result["files_complete_run_count"], 6)
        self.assertTrue(result["unexpected_run_files"])
        self.assertFalse(result["six_run_files_complete"])

    def test_mask_grid_mismatch_is_not_motion_release(self):
        runs = self.complete_fixture()
        self.image(runs[("ctr", "restpre")][1], shape=(3, 4, 5), shift=2)
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertIn("mask_does_not_match_bold_grid", result["runs"][0]["file_errors"])
        self.assertFalse(result["six_run_motion_gate_pass"])

    def test_common_standard_grid_is_checked_across_all_six_runs(self):
        runs = self.complete_fixture()
        bold, mask, _ = runs[("occ", "restpost")]
        self.image(bold, shift=2); self.image(mask, shape=(3, 4, 5), shift=2)
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertIn("bold_grid_differs_from_other_runs", result["runs"][-1]["file_errors"])
        self.assertFalse(result["six_run_files_complete"])

    def test_wrong_frame_count_and_tr_fail_header_gate(self):
        runs = self.complete_fixture()
        self.image(runs[("ctr", "restpre")][0], shape=(3, 4, 5, 599), tr=2)
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        errors = result["runs"][0]["file_errors"]
        self.assertIn("bold_volumes_599_expected_600", errors)
        self.assertIn("bold_wrong_tr_seconds:2.0", errors)
        self.assertIn("confounds_length_does_not_match_bold", errors)
        self.assertFalse(result["six_run_files_complete"])

    def test_nonfinite_affine_cannot_pass_or_break_the_failure_receipt(self):
        runs = self.complete_fixture()
        bold = runs[("ctr", "restpre")][0]
        blob = bytearray(gzip.decompress(bold.read_bytes()))
        # NIfTI-1 srow_x starts at byte 280; deliberately corrupt only a
        # synthetic fixture's selected sform, not any real imaging file.
        struct.pack_into("<f", blob, 280, float("nan"))
        bold.write_bytes(gzip.compress(blob))
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertFalse(result["six_run_files_complete"])
        self.assertTrue(any("invalid_affine" in error or "unreadable_header" in error
                            for error in result["runs"][0]["file_errors"]))
        text = (self.out / "derivative_qc.json").read_text()
        self.assertNotIn("NaN", text)
        self.assertFalse(json.loads(text)["analysis_release"])

    def test_truncated_confounds_and_duplicate_headers_are_detected(self):
        runs = self.complete_fixture()
        tsv = runs[("ctr", "restpre")][2]
        tsv.write_text("trans_x\ttrans_x\n0\t0\n", encoding="utf-8")
        result = qc.audit_subject(self.derivatives, "s001", POLICY, self.out)
        self.assertFalse(result["six_run_files_complete"])
        self.assertTrue(any("duplicate columns" in e for e in result["runs"][0]["file_errors"]))

    def test_unsafe_subject_and_output_in_readonly_derivatives_rejected(self):
        with self.assertRaises(ValueError):
            qc.audit_subject(self.derivatives, "../s001", POLICY, self.out)
        with self.assertRaises(ValueError):
            qc.audit_subject(self.derivatives, "s001", POLICY, self.derivatives / "qc")


if __name__ == "__main__":
    unittest.main(verbosity=2)
