"""Regression checks for QC exit handling and safe queue resumption."""
import copy
import unittest

import run_ds001927_cohort as runner


def qc_record(passed=True):
    return {"subject": "sub-s017", "six_run_files_complete": True,
            "six_run_motion_gate_pass": passed, "analysis_release": False,
            "state": "motion_pass_anatomy_pending" if passed else "motion_exclusion_or_review"}


def paused_state():
    return {"status": "paused_by_user", "selected_subjects": runner.SUBJECTS,
            "completed_subjects": runner.SUBJECTS[:16], "active_subject": None,
            "subject_pending_completion_reconciliation": None,
            "image_id": runner.IMAGE_ID, "analysis_release": False}


class QCExitTests(unittest.TestCase):
    def test_all_motion_pass(self):
        runner.validate_completed_qc(qc_record(), 0, "s017")

    def test_valid_motion_exclusion_does_not_stop_preprocessing(self):
        record = qc_record(False)
        original = copy.deepcopy(record)
        runner.validate_completed_qc(record, 3, "s017")
        self.assertEqual(record, original)
        self.assertFalse(record["six_run_motion_gate_pass"])
        self.assertFalse(record["analysis_release"])

    def test_incomplete_files_and_crashes_still_stop(self):
        for code in (0, 1, 2, 3, -1):
            record = qc_record(False)
            record["six_run_files_complete"] = False
            with self.subTest(code=code), self.assertRaises(RuntimeError):
                runner.validate_completed_qc(record, code, "s017")
        with self.assertRaises(RuntimeError):
            runner.validate_completed_qc(qc_record(False), 1, "s017")

    def test_inconsistent_or_wrong_receipts_stop(self):
        for record, code in [(qc_record(), 3), (qc_record(False), 0),
                             ({**qc_record(), "subject": "sub-s018"}, 0),
                             ({**qc_record(), "analysis_release": True}, 0),
                             ({**qc_record(), "state": "incomplete_or_invalid_files"}, 0)]:
            with self.subTest(record=record, code=code), self.assertRaises(RuntimeError):
                runner.validate_completed_qc(record, code, "s017")


class ResumeTests(unittest.TestCase):
    def test_resume_selects_only_five_unstarted_people(self):
        self.assertEqual(runner.validate_resume_state(paused_state(), runner.SUBJECTS),
                         ["s018", "s019", "s020", "s021", "s023"])

    def test_active_or_unreconciled_state_is_rejected(self):
        for update in [{"status": "running"}, {"active_subject": "s017"},
                       {"subject_pending_completion_reconciliation": "s017"}]:
            with self.subTest(update=update), self.assertRaises(RuntimeError):
                runner.validate_resume_state({**paused_state(), **update}, runner.SUBJECTS)

    def test_completed_subjects_must_be_unique_ordered_prefix(self):
        for completed in [["s003"], ["s002", "s002"], ["s002", "s004"]]:
            with self.subTest(completed=completed), self.assertRaises(RuntimeError):
                runner.validate_resume_state({**paused_state(), "completed_subjects": completed}, runner.SUBJECTS)

    def test_different_selection_or_image_is_rejected(self):
        with self.assertRaises(RuntimeError):
            runner.validate_resume_state(paused_state(), runner.SUBJECTS[1:])
        with self.assertRaises(RuntimeError):
            runner.validate_resume_state({**paused_state(), "image_id": "different"}, runner.SUBJECTS)


if __name__ == "__main__":
    unittest.main(verbosity=2)
