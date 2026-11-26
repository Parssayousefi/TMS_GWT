#!/usr/bin/env python3
"""
Replace skull-stripped LAM T1w images with original MPRAGE scans.

The Borgwardt group provided the original (non-skull-stripped) MPRAGE images
for 26 of 28 LAM subjects as a zip file. These replace the skull-stripped
versions (~2.7 MB) so fMRIPrep can run its own skull-stripping.

Source files:  data/basel_lam/LAM MPRAGE.zip -> LAM MPRAGE/lam_fmri_{num}_*.nii
Target files:  data/basel_lam/sub-{num}/anat/sub-{num}_T1w.nii.gz

Missing subjects (no MPRAGE in zip):
  - Subject 007: Dropout (received placebo twice)
  - Subject 009: Missing from zip
  Both keep their existing skull-stripped T1s and need --skull-strip-t1w skip
  in fMRIPrep, or should be excluded.

Usage:
    python replace_lam_t1s.py --dry-run     # Preview changes
    python replace_lam_t1s.py               # Execute replacement
"""

import argparse
import gzip
import logging
import re
import shutil
import zipfile
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

BIDS_ROOT = Path("D:/projects/TMS_GWT/data/basel_lam")
MPRAGE_ZIP = BIDS_ROOT / "LAM MPRAGE.zip"
MPRAGE_DIR = BIDS_ROOT / "LAM MPRAGE"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Core logic ───────────────────────────────────────────────────────────────


def extract_zip():
    """Extract LAM MPRAGE.zip if not already extracted."""
    if MPRAGE_DIR.exists() and any(MPRAGE_DIR.iterdir()):
        log.info(f"MPRAGE directory already exists: {MPRAGE_DIR}")
        return

    if not MPRAGE_ZIP.exists():
        log.error(f"MPRAGE zip not found: {MPRAGE_ZIP}")
        raise FileNotFoundError(f"MPRAGE zip not found: {MPRAGE_ZIP}")

    log.info(f"Extracting {MPRAGE_ZIP} ...")
    with zipfile.ZipFile(MPRAGE_ZIP, "r") as zf:
        zf.extractall(BIDS_ROOT)
    log.info(f"Extracted to {MPRAGE_DIR}")


def discover_mprage_files() -> dict[str, Path]:
    """Find all MPRAGE .nii files and map subject number -> file path.

    Returns dict like {'001': Path('lam_fmri_001_...nii'), ...}
    """
    if not MPRAGE_DIR.exists():
        log.error(f"MPRAGE directory not found: {MPRAGE_DIR}")
        return {}

    pattern = re.compile(r"^lam_fmri_(\d{3})_.*\.nii$")
    files: dict[str, Path] = {}

    for fpath in sorted(MPRAGE_DIR.iterdir()):
        if not fpath.is_file():
            continue
        m = pattern.match(fpath.name)
        if not m:
            log.warning(f"  Unexpected file in MPRAGE dir: {fpath.name}")
            continue
        sub_num = m.group(1)
        files[sub_num] = fpath

    return files


def replace_t1(sub_num: str, source: Path, dry_run: bool) -> bool:
    """Replace a single subject's T1w with the gzipped MPRAGE source.

    Reads the uncompressed .nii and writes a gzip-compressed .nii.gz.
    Returns True on success.
    """
    target_dir = BIDS_ROOT / f"sub-{sub_num}" / "anat"
    target_nii = target_dir / f"sub-{sub_num}_T1w.nii.gz"

    if not target_dir.exists():
        log.warning(f"  sub-{sub_num}: BIDS anat dir not found, skipping")
        return False

    src_size_mb = source.stat().st_size / (1024 * 1024)

    if target_nii.exists():
        old_size_mb = target_nii.stat().st_size / (1024 * 1024)
        log.info(
            f"  sub-{sub_num}: {source.name} ({src_size_mb:.1f} MB) -> "
            f"sub-{sub_num}_T1w.nii.gz (replacing {old_size_mb:.1f} MB)"
        )
    else:
        log.info(
            f"  sub-{sub_num}: {source.name} ({src_size_mb:.1f} MB) -> "
            f"sub-{sub_num}_T1w.nii.gz (new)"
        )

    if not dry_run:
        # Read uncompressed .nii and write as .nii.gz
        with open(source, "rb") as f_in:
            with gzip.open(target_nii, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Replace skull-stripped LAM T1w images with original MPRAGE scans."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview changes without modifying files.",
    )
    args = parser.parse_args()

    mode = "DRY RUN" if args.dry_run else "EXECUTE"
    log.info(f"{'=' * 60}")
    log.info(f"LAM T1w Replacement [{mode}]")
    log.info(f"BIDS root: {BIDS_ROOT}")
    log.info(f"MPRAGE zip: {MPRAGE_ZIP}")
    log.info(f"{'=' * 60}")

    # Extract zip if needed
    if not args.dry_run:
        extract_zip()
    else:
        if not MPRAGE_DIR.exists():
            log.info("(Dry run) Would extract zip first")
            # For dry run, extract anyway so we can list files
            extract_zip()

    # Discover source files
    mprage_files = discover_mprage_files()
    if not mprage_files:
        log.error("No MPRAGE files found. Aborting.")
        return

    log.info(f"Found {len(mprage_files)} MPRAGE files")

    # Replace each subject's T1w
    replaced = 0
    errors = 0

    for sub_num in sorted(mprage_files.keys()):
        source = mprage_files[sub_num]
        success = replace_t1(sub_num, source, args.dry_run)
        if success:
            replaced += 1
        else:
            errors += 1

    # Check for subjects with BIDS dirs but no MPRAGE
    bids_subjects = sorted(
        d.name.replace("sub-", "")
        for d in BIDS_ROOT.iterdir()
        if d.is_dir() and d.name.startswith("sub-")
    )
    missing = [s for s in bids_subjects if s not in mprage_files]
    skipped = len(missing)
    if missing:
        log.warning(
            f"BIDS subjects without MPRAGE source (keep skull-stripped): "
            f"{', '.join(f'sub-{s}' for s in missing)}"
        )
        log.warning(
            "  These subjects need --skull-strip-t1w skip in fMRIPrep, "
            "or should be excluded."
        )

    # Summary
    log.info(f"\n{'=' * 60}")
    log.info(f"SUMMARY [{mode}]")
    log.info(f"  Subjects replaced: {replaced}")
    log.info(f"  Subjects skipped:  {skipped}")
    log.info(f"  Errors:            {errors}")
    log.info(f"{'=' * 60}")

    if args.dry_run:
        log.info("\nThis was a DRY RUN. Re-run without --dry-run to apply changes.")


if __name__ == "__main__":
    main()
