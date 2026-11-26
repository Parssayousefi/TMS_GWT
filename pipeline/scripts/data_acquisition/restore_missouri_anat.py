#!/usr/bin/env python3
"""
Restore anatomical data for Missouri ds006072 from _backup/ to active BIDS.

The rename_missouri_sessions.py script moved unused sessions (including ses-1
which contained anatomical scans) to _backup/. This script copies T1w/T2w
files from _backup/sub-PX/ses-1/anat/ into the active sub-PX/ses-baseline/anat/
directory, renaming ses-1 -> ses-baseline in filenames.

fMRIPrep requires T1w data to run. Without this step, fMRIPrep will fail.

Usage:
    python restore_missouri_anat.py              # Dry run (preview)
    python restore_missouri_anat.py --execute    # Apply changes
"""

import argparse
import logging
import shutil
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

DATA_ROOT = Path("D:/projects/TMS_GWT/data/missouri_psilo")
BACKUP_ROOT = DATA_ROOT / "_backup"

# Per-subject anat restore configuration.
# All subjects had anat data in ses-1 (the initial enrollment/screening session).
# P5 already has anat in ses-baseline (was present in the renamed session).
# P7 anat manually downloaded from OpenNeuro ses-0 and placed in ses-baseline.
RESTORE_MAP = {
    "P1": {"source_ses": "ses-1", "target_ses": "ses-baseline"},
    "P2": {"source_ses": "ses-1", "target_ses": "ses-baseline"},
    "P3": {"source_ses": "ses-1", "target_ses": "ses-baseline"},
    "P4": {"source_ses": "ses-1", "target_ses": "ses-baseline"},
    # P4 also has anat in ses-3, ses-4 but we use ses-1 for consistency
    "P6": {"source_ses": "ses-1", "target_ses": "ses-baseline"},
}

SKIP_SUBJECTS = {"P5", "P7"}  # Already have anat in ses-baseline/anat/

ALL_SUBJECTS = ["P1", "P2", "P3", "P4", "P5", "P6", "P7"]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Core logic ───────────────────────────────────────────────────────────────


def restore_anat_files(subject: str, source_ses: str, target_ses: str,
                       execute: bool) -> int:
    """Copy anat files from backup to active directory, renaming session entity.

    Returns count of files copied.
    """
    source_dir = BACKUP_ROOT / f"sub-{subject}" / source_ses / "anat"
    target_dir = DATA_ROOT / f"sub-{subject}" / target_ses / "anat"

    if not source_dir.exists():
        log.error(f"  Source not found: {source_dir.relative_to(DATA_ROOT)}")
        return 0

    if not source_dir.is_dir():
        log.error(f"  Source is not a directory: {source_dir.relative_to(DATA_ROOT)}")
        return 0

    # Create target directory
    if execute:
        target_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"  TARGET: sub-{subject}/{target_ses}/anat/")

    count = 0
    for fpath in sorted(source_dir.iterdir()):
        if not fpath.is_file():
            continue

        old_name = fpath.name
        # Replace _ses-1_ (or whichever source session) with _ses-baseline_
        new_name = old_name.replace(f"_{source_ses}_", f"_{target_ses}_")
        target_path = target_dir / new_name

        if target_path.exists():
            log.warning(f"  SKIP (exists): {new_name}")
            continue

        size_mb = fpath.stat().st_size / (1024 * 1024)
        log.info(f"  COPY: {old_name} -> {new_name} ({size_mb:.1f} MB)")

        if execute:
            shutil.copy2(str(fpath), str(target_path))
        count += 1

    return count


def verify_anat(subject: str, target_ses: str) -> bool:
    """Check that target anat dir has at least one T1w file."""
    target_dir = DATA_ROOT / f"sub-{subject}" / target_ses / "anat"
    if not target_dir.exists():
        return False
    t1w_files = list(target_dir.glob("*_T1w.nii*"))
    return len(t1w_files) > 0


def main():
    parser = argparse.ArgumentParser(
        description="Restore anatomical data from _backup/ to active BIDS directory."
    )
    parser.add_argument(
        "--execute", action="store_true",
        help="Actually copy files. Default is dry-run."
    )
    args = parser.parse_args()

    mode = "EXECUTE" if args.execute else "DRY RUN"
    log.info(f"{'='*60}")
    log.info(f"Missouri Anat Restore [{mode}]")
    log.info(f"Data root: {DATA_ROOT}")
    log.info(f"Backup root: {BACKUP_ROOT}")
    log.info(f"{'='*60}")

    if not DATA_ROOT.exists():
        log.error(f"Data root not found: {DATA_ROOT}")
        return

    if not BACKUP_ROOT.exists():
        log.error(f"Backup root not found: {BACKUP_ROOT}")
        return

    restored = 0
    skipped = 0
    warned = 0
    total_files = 0

    for subject in ALL_SUBJECTS:
        log.info(f"\n--- sub-{subject} ---")

        if subject in SKIP_SUBJECTS:
            has_anat = verify_anat(subject, "ses-baseline")
            if has_anat:
                log.info(f"  SKIP: Already has anat in ses-baseline/anat/")
                skipped += 1
            else:
                log.warning(f"  Expected anat in ses-baseline/anat/ but not found!")
                warned += 1
            continue

        if subject not in RESTORE_MAP:
            log.warning(f"  No restore mapping defined")
            warned += 1
            continue

        mapping = RESTORE_MAP[subject]
        source_ses = mapping["source_ses"]
        target_ses = mapping["target_ses"]

        log.info(f"  SOURCE: _backup/sub-{subject}/{source_ses}/anat/")
        n_files = restore_anat_files(subject, source_ses, target_ses, args.execute)
        total_files += n_files

        if n_files > 0:
            restored += 1
            # Verify after copy (only in execute mode)
            if args.execute:
                if verify_anat(subject, target_ses):
                    log.info(f"  VERIFIED: T1w file present in target")
                else:
                    log.error(f"  VERIFY FAILED: No T1w file found after copy!")
        else:
            log.warning(f"  No files copied for sub-{subject}")
            warned += 1

    # Summary
    log.info(f"\n{'='*60}")
    log.info(f"SUMMARY [{mode}]")
    log.info(f"  Subjects restored:  {restored}")
    log.info(f"  Subjects skipped:   {skipped}")
    log.info(f"  Subjects warned:    {warned}")
    log.info(f"  Total files copied: {total_files}")
    log.info(f"{'='*60}")

    if not args.execute:
        log.info("\nThis was a DRY RUN. Re-run with --execute to apply changes.")


if __name__ == "__main__":
    main()
