"""Publish an aggregate receipt after the pinned dataset download is complete."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "output/ds001927_download_20261126"


def main():
    status = json.loads((SOURCE / "status.json").read_text())
    assert status["state"] == "complete" and not status["failed"]
    manifest = json.loads((SOURCE / "download_manifest.json").read_text())
    verified = json.loads((SOURCE / "verified_files.json").read_text())
    headers = json.loads((SOURCE / "nifti_header_checks.json").read_text())
    assert status["commit"] == manifest["commit"] == "3239a60164f7a230e02b6a42eb36a13ff4f97f57"
    assert len(verified) == status["completed_files"] == status["total_files"] == 433
    assert len(headers) == status["nifti_headers_verified"] == 207
    assert all(v["verified"] for v in verified)
    assert {v["path"]: v["expected_digest"] for v in verified} == {v["path"]: v["expected_digest"] for v in manifest["files"]}
    assert sum(v["bytes"] for v in verified) == status["verified_bytes"] == 22968601962
    for v in verified:
        assert (ROOT / "data/ds001927" / v["path"]).stat().st_size == v["bytes"]
    assert not list((ROOT / "data/ds001927").rglob("*.part"))
    published = ROOT / "PROVENANCE/additional_datasets_20261126/download_complete.json"
    receipt = {k: status[k] for k in ["state", "dataset", "version", "commit", "total_files", "verified_bytes",
               "completed_utc", "nifti_headers_verified", "rest_runs", "t1w_images", "participants",
               "sessions_per_participant", "content_checks", "quality_limit"]}
    receipt.update(local_data_directory="data/ds001927", imaging_committed_to_git=False,
                   preprocessing_started=False, pid_analysis_started=False,
                   source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [SOURCE / n for n in ["status.json", "download_manifest.json", "verified_files.json", "nifti_header_checks.json"]] + [Path(__file__), ROOT / "setup/download_ds001927.py"]})
    published.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    replacements = {
        "README.md": [("file headers were verified. Its full pinned release is now downloading to\n`data/ds001927/`, with per-file checksum verification; it has not been analysed.",
                       "file headers were verified. The full pinned release is downloaded in\n`data/ds001927/`: 433 verified files, 22.97 GB, including 138 resting runs and\n69 structural images. Preprocessing and PID analysis have not started. See the\n[download receipt](PROVENANCE/additional_datasets_20261126/download_complete.json).")],
        "docs/TMS_CURRENT_ANALYSIS_STATUS_20261126.md": [("full version 2.1.0 release is now downloading to `data/ds001927/` under the\nuser's instruction, with pinned content checksums. It has not been analysed.",
                       "full version 2.1.0 release is downloaded in `data/ds001927/`: 433 files,\n22.97 GB, all content hashes verified, with 138 resting runs and 69 structural\nimages. All 207 NIfTI headers were checked. It has not been preprocessed or\nanalysed; motion, registration and signal coverage still need QC.")],
        "docs/TMS_ADDITIONAL_DATASETS_20261126.md": [("the full ds001927 release. It is now downloading to `data/ds001927/` using\n`setup/download_ds001927.py`; the manifest is pinned to version 2.1.0 and\nverifies Git metadata hashes and annex content digests. Interrupted transfers\nare resumable. Current transfer status is recorded in\n`output/ds001927_download_20261126/status.json`.",
                       "the full ds001927 release. **The download is complete** in `data/ds001927/`:\n433 files and 22,968,601,962 bytes including metadata and the deposited study\ncode. Git metadata hashes and annex content digests all matched the pinned\nversion 2.1.0 manifest; SHA-256 hashes were also recorded for every local file.\nAll 207 NIfTI headers were checked, including 600 volumes at TR 1.25 s for each\nof the 138 resting runs. The study code was downloaded but not executed.\nThe public aggregate receipt is\n`PROVENANCE/additional_datasets_20261126/download_complete.json`.")]
    }
    for relative, edits in replacements.items():
        path = ROOT / relative
        text = path.read_text(encoding="utf-8")
        for old, new in edits:
            assert old in text or new in text, (relative, "Unrecognized documentation state")
            text = text.replace(old, new)
        path.write_bytes(text.encode("utf-8"))
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
