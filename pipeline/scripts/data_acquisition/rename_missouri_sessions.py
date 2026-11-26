#!/usr/bin/env python3
"""
Rename Missouri ds006072 BIDS sessions from numeric IDs to condition names.

Renames ses-X directories and all files within to use ses-{condition} names
(ses-baseline, ses-psilocybin, ses-methylphenidate), so the dataset matches
the pipeline's Pattern A convention (session name = condition).

Unused baseline sessions are moved to _backup/. Run with --execute to apply;
default is dry-run (preview only).

Usage:
    python rename_missouri_sessions.py              # Dry run
    python rename_missouri_sessions.py --execute    # Apply changes
"""

import argparse
import logging
import re
import shutil
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

DATA_ROOT = Path("D:/projects/TMS_GWT/data/missouri_psilo")
BACKUP_ROOT = DATA_ROOT / "_backup"

# Per-subject rename map: old numeric session → new condition name
# Baseline = last pre-drug session (temporally closest to drug administration)
# Triple-verified from README ORDER table, session_data.csv, and Excel Key sheet
RENAME_MAP = {
    "P1": {
        "ses-5":  "ses-baseline",
        "ses-11": "ses-psilocybin",
        "ses-6":  "ses-methylphenidate",
    },
    "P2": {
        "ses-9":  "ses-baseline",
        "ses-10": "ses-psilocybin",
        "ses-16": "ses-methylphenidate",
    },
    "P3": {
        "ses-7":  "ses-baseline",
        "ses-12": "ses-psilocybin",
        "ses-8":  "ses-methylphenidate",
    },
    "P4": {
        "ses-5":  "ses-baseline",
        "ses-10": "ses-psilocybin",
        "ses-6":  "ses-methylphenidate",
    },
    "P5": {
        "ses-5":  "ses-baseline",
        "ses-6":  "ses-psilocybin",
        "ses-10": "ses-methylphenidate",
    },
    "P6": {
        "ses-4":  "ses-baseline",
        "ses-9":  "ses-psilocybin",
        "ses-5":  "ses-methylphenidate",
    },
    "P7": {
        "ses-6":  "ses-baseline",
        "ses-7":  "ses-psilocybin",
        "ses-11": "ses-methylphenidate",
    },
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Core logic ───────────────────────────────────────────────────────────────


def rename_files_in_dir(session_dir: Path, old_ses: str, new_ses: str,
                        execute: bool) -> int:
    """Rename all files in a session directory, updating the _ses-X_ entity.

    Returns count of files renamed.
    """
    count = 0
    # Walk all subdirectories (func/, anat/, fmap/)
    for fpath in sorted(session_dir.rglob("*")):
        if not fpath.is_file():
            continue
        old_name = fpath.name
        # Replace _ses-X_ with _ses-{condition}_ in filename
        new_name = old_name.replace(f"_{old_ses}_", f"_{new_ses}_")
        if new_name == old_name:
            continue
        new_path = fpath.parent / new_name
        log.info(f"  FILE: {fpath.relative_to(DATA_ROOT)} -> {new_name}")
        if execute:
            fpath.rename(new_path)
        count += 1
    return count


def rename_session(subject: str, old_ses: str, new_ses: str,
                   execute: bool) -> bool:
    """Rename a session directory and its contents.

    Returns True if successful (or would be in dry-run).
    """
    sub_dir = DATA_ROOT / f"sub-{subject}"
    old_dir = sub_dir / old_ses
    new_dir = sub_dir / new_ses

    if not old_dir.exists():
        log.warning(f"  SKIP: {old_dir.relative_to(DATA_ROOT)} does not exist")
        return False

    if new_dir.exists():
        log.warning(f"  SKIP: {new_dir.relative_to(DATA_ROOT)} already exists "
                     "(already renamed?)")
        return False

    log.info(f"  DIR:  sub-{subject}/{old_ses} -> sub-{subject}/{new_ses}")

    # Rename files first (before moving directory)
    n_files = rename_files_in_dir(old_dir, old_ses, new_ses, execute)

    # Rename directory
    if execute:
        old_dir.rename(new_dir)

    log.info(f"  ({n_files} files renamed)")
    return True


def backup_session(subject: str, ses_name: str, execute: bool) -> bool:
    """Move an unused session to _backup/."""
    src = DATA_ROOT / f"sub-{subject}" / ses_name
    if not src.exists():
        return False

    dst_parent = BACKUP_ROOT / f"sub-{subject}"
    dst = dst_parent / ses_name

    log.info(f"  BACKUP: sub-{subject}/{ses_name} -> _backup/sub-{subject}/{ses_name}")
    if execute:
        dst_parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
    return True


def process_subject(subject: str, renames: dict, execute: bool) -> dict:
    """Process one subject: rename selected sessions, backup the rest.

    Returns stats dict.
    """
    sub_dir = DATA_ROOT / f"sub-{subject}"
    if not sub_dir.exists():
        log.warning(f"Subject directory not found: sub-{subject}")
        return {"renamed": 0, "backed_up": 0, "skipped": 0}

    stats = {"renamed": 0, "backed_up": 0, "skipped": 0}

    # Step 1: Rename selected sessions
    for old_ses, new_ses in renames.items():
        ok = rename_session(subject, old_ses, new_ses, execute)
        if ok:
            stats["renamed"] += 1
        else:
            stats["skipped"] += 1

    # Step 2: Backup remaining numeric sessions (not renamed, not already named)
    selected_old = set(renames.keys())
    # After renaming, check what numeric sessions remain
    ses_pattern = re.compile(r"^ses-\d+$")
    for entry in sorted(sub_dir.iterdir()):
        if not entry.is_dir():
            continue
        if not ses_pattern.match(entry.name):
            continue  # Skip already-named sessions (ses-baseline, etc.)
        if entry.name in selected_old:
            continue  # Should have been renamed already (only if rename failed)
        ok = backup_session(subject, entry.name, execute)
        if ok:
            stats["backed_up"] += 1

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Rename Missouri BIDS sessions from numeric to condition names."
    )
    parser.add_argument(
        "--execute", action="store_true",
        help="Actually rename files/directories. Default is dry-run."
    )
    args = parser.parse_args()

    mode = "EXECUTE" if args.execute else "DRY RUN"
    log.info(f"{'='*60}")
    log.info(f"Missouri Session Renaming [{mode}]")
    log.info(f"Data root: {DATA_ROOT}")
    log.info(f"{'='*60}")

    if not DATA_ROOT.exists():
        log.error(f"Data root not found: {DATA_ROOT}")
        return

    total = {"renamed": 0, "backed_up": 0, "skipped": 0}

    for subject, renames in sorted(RENAME_MAP.items()):
        log.info(f"\n--- sub-{subject} ---")
        for old_ses, new_ses in renames.items():
            log.info(f"  Plan: {old_ses} -> {new_ses}")

        stats = process_subject(subject, renames, args.execute)
        for k in total:
            total[k] += stats[k]

    # Summary
    log.info(f"\n{'='*60}")
    log.info(f"SUMMARY [{mode}]")
    log.info(f"  Sessions renamed:   {total['renamed']}")
    log.info(f"  Sessions backed up: {total['backed_up']}")
    log.info(f"  Sessions skipped:   {total['skipped']}")
    log.info(f"{'='*60}")

    if not args.execute:
        log.info("\nThis was a DRY RUN. Re-run with --execute to apply changes.")


if __name__ == "__main__":
    main()
