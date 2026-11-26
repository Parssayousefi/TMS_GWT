"""Regression checks for the observed Windows progress-file sharing violation."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from run_tms_processing import ROOT, write_json


class StatusWriteTests(unittest.TestCase):
    def test_transient_lock_retries_atomic_replace(self):
        temporary_root = (ROOT / 'tmp').resolve()
        with tempfile.TemporaryDirectory(dir=temporary_root) as folder:
            self.assertIn(temporary_root, Path(folder).resolve().parents)
            target = Path(folder) / 'status.json'
            target.write_text('{"status": "old"}')
            original_replace = Path.replace
            attempts = []

            def replace(source, destination):
                attempts.append(source)
                if len(attempts) <= 2:
                    self.assertEqual(json.loads(target.read_text()), {'status': 'old'})
                    raise PermissionError(5, 'simulated Windows sharing violation')
                return original_replace(source, destination)

            with patch.object(Path, 'replace', autospec=True, side_effect=replace), patch('run_tms_processing.time.sleep') as sleep:
                write_json(target, {'status': 'new', 'scans': 193})
            self.assertEqual(len(attempts), 3)
            self.assertEqual(sleep.call_count, 2)
            self.assertEqual(json.loads(target.read_text()), {'status': 'new', 'scans': 193})
            self.assertFalse(target.with_suffix('.json.tmp').exists())

    def test_persistent_lock_is_bounded_and_preserves_previous_state(self):
        temporary_root = (ROOT / 'tmp').resolve()
        with tempfile.TemporaryDirectory(dir=temporary_root) as folder:
            self.assertIn(temporary_root, Path(folder).resolve().parents)
            target = Path(folder) / 'status.json'
            target.write_text('{"status": "old"}')
            with patch.object(Path, 'replace', side_effect=PermissionError(5, 'still locked')) as replace, patch('run_tms_processing.time.sleep') as sleep:
                with self.assertRaises(PermissionError):
                    write_json(target, {'status': 'new'})
            self.assertEqual(replace.call_count, 20)
            self.assertEqual(sleep.call_count, 19)
            self.assertEqual(json.loads(target.read_text()), {'status': 'old'})


if __name__ == '__main__':
    unittest.main()
