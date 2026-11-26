#!/usr/bin/env python3
"""
Selective download of ds003171 (Cambridge Propofol dataset) from OpenNeuro.

Dataset: "Modeling an auditory stimulated brain under altered states of consciousness
          using the generalized ising model"
URL: https://openneuro.org/datasets/ds003171/versions/2.0.1

This is the propofol anaesthesia dataset used in Luppi et al. (2024) eLife.
17 healthy adults scanned at 4 sedation levels: Awake, Light, Deep, Recovery.

Downloads:
  - Resting-state BOLD (task-rest*) func files
  - Audio task BOLD (task-audio*) func files
  - T1w anatomical images
  - Associated JSON sidecars and events files
  - Skips DWI (not needed for fMRI pipeline)
"""

import subprocess
import sys
import time
import logging
from pathlib import Path

# -- Configuration -----------------------------------------------------------

S3_BUCKET = "s3://openneuro.org/ds003171"
LOCAL_ROOT = Path("D:/projects/TMS_GWT/data/cambridge_propofol")

# All 17 subjects in the dataset
SUBJECTS = [
    "sub-02CB", "sub-04HD", "sub-04SG", "sub-08BC", "sub-08VR",
    "sub-10JR", "sub-13CA", "sub-16RA", "sub-17EK", "sub-17NA",
    "sub-19AK", "sub-19SA", "sub-22CY", "sub-22TK", "sub-2525JK",
    "sub-27HF", "sub-30AQ",
]

# Top-level files to download
TOP_LEVEL_FILES = [
    "dataset_description.json",
    "README",
    "CHANGES",
    "participants.tsv",
]

# -- Logging -----------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# -- Functions ---------------------------------------------------------------


def run_s3_sync(s3_path: str, local_path: str, include_filters: list[str],
                max_retries: int = 3) -> bool:
    """Run aws s3 sync with include/exclude filters and retry logic."""
    cmd = [
        "aws", "s3", "sync",
        s3_path, local_path,
        "--no-sign-request",
        "--exclude", "*",
    ]
    for pattern in include_filters:
        cmd.extend(["--include", pattern])

    for attempt in range(1, max_retries + 1):
        log.info(f"  Syncing (attempt {attempt}/{max_retries}): {s3_path}")
        log.info(f"  Filters: {include_filters}")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200,
            )
            if result.returncode == 0:
                if result.stdout.strip():
                    for line in result.stdout.strip().split("\n")[-5:]:
                        log.info(f"    {line}")
                return True
            else:
                log.warning(f"  aws s3 sync failed (code {result.returncode}): "
                            f"{result.stderr.strip()}")
        except subprocess.TimeoutExpired:
            log.warning(f"  Timeout on attempt {attempt}")
        except Exception as e:
            log.warning(f"  Error on attempt {attempt}: {e}")

        if attempt < max_retries:
            wait = 30 * attempt
            log.info(f"  Retrying in {wait}s...")
            time.sleep(wait)

    log.error(f"  FAILED after {max_retries} attempts: {s3_path}")
    return False


def download_top_level_files():
    """Download top-level dataset files."""
    log.info("Downloading top-level dataset files...")
    LOCAL_ROOT.mkdir(parents=True, exist_ok=True)

    for fname in TOP_LEVEL_FILES:
        local_file = LOCAL_ROOT / fname
        if local_file.exists():
            log.info(f"  Already exists: {fname}")
            continue
        s3_uri = f"{S3_BUCKET}/{fname}"
        cmd = ["aws", "s3", "cp", s3_uri, str(local_file), "--no-sign-request"]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            if result.returncode == 0:
                log.info(f"  Downloaded: {fname}")
            else:
                log.warning(f"  Failed to download {fname}: {result.stderr.strip()}")
        except Exception as e:
            log.warning(f"  Error downloading {fname}: {e}")


def download_subject(subject: str) -> bool:
    """Download anat + func for a single subject (skips DWI).

    Downloads:
      - anat/*T1w* (structural)
      - func/*rest* (resting-state BOLD at all sedation levels)
      - func/*audio* (auditory task BOLD at all sedation levels)
      - scans.tsv
    """
    s3_path = f"{S3_BUCKET}/{subject}/"
    local_path = str(LOCAL_ROOT / subject)

    log.info(f"{'=' * 60}")
    log.info(f"Downloading {subject}")
    log.info(f"{'=' * 60}")

    include_filters = [
        # T1w anatomical (nii.gz + json)
        "anat/*T1w*",
        # Resting-state BOLD at all sedation levels
        "func/*rest*",
        # Audio task BOLD at all sedation levels
        "func/*audio*",
        # Scans file
        "*.tsv",
    ]

    return run_s3_sync(s3_path, local_path, include_filters)


def download_all():
    """Download all subjects."""
    download_top_level_files()

    total = len(SUBJECTS)
    completed = 0
    failed = []

    log.info(f"\nTotal subjects to download: {total}")
    log.info(f"Subjects: {SUBJECTS}\n")

    for subject in SUBJECTS:
        ok = download_subject(subject)
        completed += 1
        if not ok:
            failed.append(subject)
        log.info(f"Progress: {completed}/{total} subjects done\n")

    # Summary
    log.info(f"\n{'=' * 60}")
    log.info("DOWNLOAD COMPLETE")
    log.info(f"{'=' * 60}")
    log.info(f"Subjects attempted: {total}")
    log.info(f"Subjects failed: {len(failed)}")
    if failed:
        log.warning(f"Failed subjects: {failed}")
    log.info(f"\nData saved to: {LOCAL_ROOT}")

    return len(failed) == 0


if __name__ == "__main__":
    success = download_all()
    sys.exit(0 if success else 1)
