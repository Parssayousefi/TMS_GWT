#!/usr/bin/env python3
"""
Rename Monash PsiConnect (ds006110) BIDS sessions from numeric IDs to condition names.

Renames ses-01 -> ses-baseline and ses-02 -> ses-psilocybin for every sub-PC*
subject (directories and the _ses-XX_ entity in all filenames), so the dataset
matches the pipeline's Pattern A convention (session name = condition), like
Berlin Ketamine and Missouri ds006072.

Unlike rename_missouri_sessions.py there is no per-subject map and no backup
step: the mapping is uniform and every session is kept. Subjects without a
ses-02 (PC016, PC024, PC232) simply have nothing to rename there.

scans.tsv contents are left stale (filenames inside them still reference the
old session labels), matching the Missouri precedent; fMRIPrep does not need
them.

Usage:
    python rename_monash_sessions.py              # Dry run
    python rename_monash_sessions.py --execute    # Apply changes
"""

import argparse
import logging
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

DATA_ROOT = Path("D:/projects/TMS_GWT/data/monash_psilo")

RENAME_MAP = {
    "ses-01": "ses-baseline",
    "ses-02": "ses-psilocybin",
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
    """Rename all files in a session directory, updating the _ses-XX_ entity.

    Returns count of files renamed.
    """
    count = 0
    for fpath in sorted(session_dir.rglob("*")):
        if not fpath.is_file():
            continue
        old_name = fpath.name
        new_name = old_name.replace(f"_{old_ses}_", f"_{new_ses}_")
        if new_name == old_name:
            continue
        new_path = fpath.parent / new_name
        log.info(f"  FILE: {fpath.relative_to(DATA_ROOT)} -> {new_name}")
        if execute:
            fpath.rename(new_path)
        count += 1
    return count


def rename_session(sub_dir: Path, old_ses: str, new_ses: str,
                   execute: bool) -> bool:
    """Rename a session directory and its contents.

    Returns True if renamed (or would be in dry-run), False if skipped.
    """
    old_dir = sub_dir / old_ses
    new_dir = sub_dir / new_ses

    if not old_dir.exists():
        return False

    if new_dir.exists():
        log.warning(f"  SKIP: {new_dir.relative_to(DATA_ROOT)} already exists "
                    "(already renamed?)")
        return False

    log.info(f"  DIR:  {sub_dir.name}/{old_ses} -> {sub_dir.name}/{new_ses}")

    n_files = rename_files_in_dir(old_dir, old_ses, new_ses, execute)

    if execute:
        import time
        for attempt in range(5):
            try:
                old_dir.rename(new_dir)
                break
            except PermissionError:
                if attempt == 4:
                    log.error(f"  FAILED to rename dir {old_dir.relative_to(DATA_ROOT)} "
                              "(WinError 5 — locked by another process?). "
                              "Files inside ARE renamed; re-run to retry the dir.")
                    return False
                time.sleep(2)

    log.info(f"  ({n_files} files renamed)")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Rename Monash BIDS sessions from numeric to condition names."
    )
    parser.add_argument(
        "--execute", action="store_true",
        help="Actually rename files/directories. Default is dry-run."
    )
    args = parser.parse_args()

    mode = "EXECUTE" if args.execute else "DRY RUN"
    log.info("=" * 60)
    log.info(f"Monash PsiConnect Session Renaming [{mode}]")
    log.info(f"Data root: {DATA_ROOT}")
    log.info("=" * 60)

    if not DATA_ROOT.exists():
        log.error(f"Data root not found: {DATA_ROOT}")
        return

    total = {"renamed": 0, "skipped": 0, "missing": 0}

    for sub_dir in sorted(DATA_ROOT.glob("sub-PC*")):
        if not sub_dir.is_dir():
            continue
        log.info(f"\n--- {sub_dir.name} ---")
        for old_ses, new_ses in RENAME_MAP.items():
            if not (sub_dir / old_ses).exists() and not (sub_dir / new_ses).exists():
                log.info(f"  (no {old_ses})")
                total["missing"] += 1
                continue
            ok = rename_session(sub_dir, old_ses, new_ses, args.execute)
            if ok:
                total["renamed"] += 1
            else:
                total["skipped"] += 1

    log.info("\n" + "=" * 60)
    log.info(f"SUMMARY [{mode}]")
    log.info(f"  Sessions renamed:  {total['renamed']}")
    log.info(f"  Sessions skipped:  {total['skipped']}")
    log.info(f"  Sessions absent:   {total['missing']}")
    log.info("=" * 60)

    if not args.execute:
        log.info("\nThis was a DRY RUN. Re-run with --execute to apply changes.")


if __name__ == "__main__":
    main()
