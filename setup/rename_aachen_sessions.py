#!/usr/bin/env python3
"""
Rename the Aachen iTBS (OpenNeuro ds001832) BIDS sessions from numeric IDs to condition names.

Session key from the dataset's derivatives/README.txt:
    ses-0 = anatomy scan
    ses-1 = pre 90% TMS rest      ses-2 = post 90% TMS rest
    ses-3 = pre 120% TMS rest     ses-4 = post 120% TMS rest

Renames the session directories and the _ses-X_ entity in every filename, so the
dataset follows the pipeline's Pattern A convention (session name = condition), like
Monash and Missouri. Same logic as pipeline/scripts/data_acquisition/rename_monash_sessions.py
(copied from GWT); the mapping is uniform across subjects.

Usage:
    python setup/rename_aachen_sessions.py              # Dry run
    python setup/rename_aachen_sessions.py --execute    # Apply changes
"""

import argparse
import logging
import time
from pathlib import Path

DATA_ROOT = Path("D:/projects/TMS_GWT/data/aachen_itbs")

RENAME_MAP = {
    "ses-0": "ses-anat",
    "ses-1": "ses-pre90",
    "ses-2": "ses-post90",
    "ses-3": "ses-pre120",
    "ses-4": "ses-post120",
}

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)


def rename_files_in_dir(session_dir: Path, old_ses: str, new_ses: str, execute: bool) -> int:
    count = 0
    for fpath in sorted(session_dir.rglob("*")):
        if not fpath.is_file():
            continue
        new_name = fpath.name.replace(f"_{old_ses}_", f"_{new_ses}_")
        if new_name == fpath.name:
            continue
        log.info(f"  FILE: {fpath.relative_to(DATA_ROOT)} -> {new_name}")
        if execute:
            fpath.rename(fpath.parent / new_name)
        count += 1
    return count


def rename_session(sub_dir: Path, old_ses: str, new_ses: str, execute: bool) -> bool:
    old_dir, new_dir = sub_dir / old_ses, sub_dir / new_ses
    if not old_dir.exists():
        return False
    if new_dir.exists():
        log.warning(f"  SKIP: {new_dir.relative_to(DATA_ROOT)} already exists (already renamed?)")
        return False
    log.info(f"  DIR:  {sub_dir.name}/{old_ses} -> {sub_dir.name}/{new_ses}")
    n_files = rename_files_in_dir(old_dir, old_ses, new_ses, execute)
    if execute:
        for attempt in range(5):
            try:
                old_dir.rename(new_dir)
                break
            except PermissionError:
                if attempt == 4:
                    log.error(f"  FAILED to rename {old_dir.relative_to(DATA_ROOT)} (locked?); "
                              "files inside ARE renamed; re-run to retry the directory.")
                    return False
                time.sleep(2)
    log.info(f"  ({n_files} files renamed)")
    return True


def main():
    ap = argparse.ArgumentParser(description="Rename ds001832 sessions to condition names.")
    ap.add_argument("--execute", action="store_true", help="Apply changes (default: dry run)")
    args = ap.parse_args()
    mode = "EXECUTE" if args.execute else "DRY RUN"
    log.info(f"Aachen iTBS session renaming [{mode}]  data root: {DATA_ROOT}")
    if not DATA_ROOT.exists():
        log.error(f"Data root not found: {DATA_ROOT}")
        return
    total = {"renamed": 0, "skipped": 0, "missing": 0}
    for sub_dir in sorted(DATA_ROOT.glob("sub-*")):
        if not sub_dir.is_dir():
            continue
        log.info(f"--- {sub_dir.name} ---")
        for old_ses, new_ses in RENAME_MAP.items():
            if not (sub_dir / old_ses).exists() and not (sub_dir / new_ses).exists():
                log.info(f"  (no {old_ses})")
                total["missing"] += 1
                continue
            total["renamed" if rename_session(sub_dir, old_ses, new_ses, args.execute) else "skipped"] += 1
    log.info(f"SUMMARY [{mode}]: renamed {total['renamed']}, skipped {total['skipped']}, "
             f"absent {total['missing']}")
    if not args.execute:
        log.info("This was a DRY RUN. Re-run with --execute to apply changes.")


if __name__ == "__main__":
    main()
