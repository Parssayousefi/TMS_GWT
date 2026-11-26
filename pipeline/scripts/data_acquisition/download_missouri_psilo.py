#!/usr/bin/env python3
"""
Selective download of ds006072 (Missouri Psilocybin dataset) from OpenNeuro.

Downloads only:
  - Resting-state BOLD (BOLDREST) func files (all echoes)
  - T1w anatomical images
  - Fieldmaps (fmap)
  - Associated JSON sidecars

For baseline + drug sessions only, for 7 main subjects (P1-P7).
Skips task fMRI, DWI, NON_BIDS, BehavioralAssessments, replication subjects.
"""

import subprocess
import sys
import time
import logging
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

S3_BUCKET = "s3://openneuro.org/ds006072"
LOCAL_ROOT = Path("D:/projects/TMS_GWT/data/missouri_psilo")

# Sessions to download per subject: baseline sessions + drug sessions
# Drug sessions annotated with drug type for reference
SESSION_MAP = {
    "sub-P1": {
        "baseline": [1, 2, 3, 4, 5],
        "drug": [6, 11],  # ses-6=MTP, ses-11=PSIL
    },
    "sub-P2": {
        "baseline": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "drug": [10, 16],  # ses-10=PSIL, ses-16=MTP
    },
    "sub-P3": {
        "baseline": [1, 2, 3, 4, 5, 6, 7],
        "drug": [8, 12],  # ses-8=MTP, ses-12=PSIL
    },
    "sub-P4": {
        "baseline": [1, 2, 3, 4, 5],
        "drug": [6, 10],  # ses-6=MTP, ses-10=PSIL
    },
    "sub-P5": {
        "baseline": [1, 2, 3, 4, 5],
        "drug": [6, 10],  # ses-6=PSIL, ses-10=MTP
    },
    "sub-P6": {
        "baseline": [1, 2, 3, 4],
        "drug": [5, 9],  # ses-5=MTP, ses-9=PSIL
    },
    "sub-P7": {
        "baseline": [4, 5, 6],  # ses-1,2,3 are test/NORDIC, not baseline
        "drug": [7, 11],  # ses-7=PSIL, ses-11=MTP
    },
}

# Top-level files to download
TOP_LEVEL_FILES = [
    "dataset_description.json",
    "README",
    "CHANGES",
    "LICENSE",
    "session_data.csv",
    ".bidsignore",
    "PPFM_session_notes jss 20250402.xlsx",
    "2025-02-26_11h26m26_niix2bids.log",
]

# ── Logging ──────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Functions ────────────────────────────────────────────────────────────────


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
                timeout=7200,  # 2 hour timeout per sync
            )
            if result.returncode == 0:
                if result.stdout.strip():
                    for line in result.stdout.strip().split("\n")[-5:]:
                        log.info(f"    {line}")
                return True
            else:
                log.warning(f"  aws s3 sync failed (code {result.returncode}): {result.stderr.strip()}")
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
    """Download top-level dataset files (README, dataset_description.json, etc.)."""
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


def download_session(subject: str, ses_num: int, ses_type: str) -> bool:
    """Download a single subject-session with selective filters.

    Downloads:
      - func/*BOLDREST* (resting-state BOLD, all echoes + JSON)
      - anat/*T1w* (structural, + JSON)
      - fmap/* (fieldmaps, all files)
    """
    ses_label = f"ses-{ses_num}"
    s3_path = f"{S3_BUCKET}/{subject}/{ses_label}/"
    local_path = str(LOCAL_ROOT / subject / ses_label)

    log.info(f"{'='*60}")
    log.info(f"Downloading {subject}/{ses_label} ({ses_type})")
    log.info(f"{'='*60}")

    # Include filters for selective download
    include_filters = [
        # Resting-state BOLD (all echoes, .nii and .json)
        "func/*BOLDREST*",
        # T1w anatomical
        "anat/*T1w*",
        # All fieldmaps
        "fmap/*",
    ]

    return run_s3_sync(s3_path, local_path, include_filters)


def download_all():
    """Download all subjects and sessions defined in SESSION_MAP."""
    # Download top-level files first
    download_top_level_files()

    total_sessions = sum(
        len(v["baseline"]) + len(v["drug"])
        for v in SESSION_MAP.values()
    )
    completed = 0
    failed = []

    log.info(f"\nTotal sessions to download: {total_sessions}")
    log.info(f"Subjects: {list(SESSION_MAP.keys())}\n")

    for subject, sessions in SESSION_MAP.items():
        log.info(f"\n{'#'*60}")
        log.info(f"# Subject: {subject}")
        log.info(f"# Baseline sessions: {sessions['baseline']}")
        log.info(f"# Drug sessions: {sessions['drug']}")
        log.info(f"{'#'*60}\n")

        for ses_num in sessions["baseline"]:
            ok = download_session(subject, ses_num, "baseline")
            completed += 1
            if not ok:
                failed.append(f"{subject}/ses-{ses_num}")
            log.info(f"Progress: {completed}/{total_sessions} sessions done\n")

        for ses_num in sessions["drug"]:
            ok = download_session(subject, ses_num, "drug")
            completed += 1
            if not ok:
                failed.append(f"{subject}/ses-{ses_num}")
            log.info(f"Progress: {completed}/{total_sessions} sessions done\n")

    # Summary
    log.info(f"\n{'='*60}")
    log.info("DOWNLOAD COMPLETE")
    log.info(f"{'='*60}")
    log.info(f"Sessions attempted: {total_sessions}")
    log.info(f"Sessions failed: {len(failed)}")
    if failed:
        log.warning(f"Failed sessions: {failed}")

    return len(failed) == 0


if __name__ == "__main__":
    success = download_all()
    sys.exit(0 if success else 1)
