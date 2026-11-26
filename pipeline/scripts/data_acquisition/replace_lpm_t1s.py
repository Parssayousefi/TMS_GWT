#!/usr/bin/env python3
"""
Replace skull-stripped LPM T1w images with original MPRAGE scans.

The Borgwardt group provided the original (non-skull-stripped) MPRAGE images
for all 33 LPM subjects. These replace the skull-stripped versions that were
used during initial BIDS conversion, so fMRIPrep can run its own skull-stripping.

Source files:  data/basel_lpm_bids/LPM MPRAGE/LPM MPRAGE/LPM^{num}_MPRAGE_*.nii.gz
Target files:  data/basel_lpm_bids/sub-{num}/anat/sub-{num}_T1w.nii.gz

Special cases:
  - Subject 023: Does not exist (expected)
  - Subject 021: Series 0009 instead of 0005 (different series, same acquisition)
  - Subject 031: Two files (0005 and 0005a); use 0005a (re-scan, better quality)

Usage:
    python replace_lpm_t1s.py --dry-run     # Preview changes
    python replace_lpm_t1s.py               # Execute replacement
"""

import argparse
import json
import logging
import re
import shutil
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────

BIDS_ROOT = Path("D:/projects/TMS_GWT/data/basel_lpm_bids")
MPRAGE_DIR = BIDS_ROOT / "LPM MPRAGE" / "LPM MPRAGE"

# For subject 031, prefer the re-scan (0005a) over the original (0005)
SUBJECT_031_PREFERRED = "0005a"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Core logic ───────────────────────────────────────────────────────────────


def discover_mprage_files() -> dict[str, Path]:
    """Find all MPRAGE files and map subject number -> file path.

    Returns dict like {'001': Path('LPM^001_MPRAGE_...nii.gz'), ...}
    """
    if not MPRAGE_DIR.exists():
        log.error(f"MPRAGE directory not found: {MPRAGE_DIR}")
        return {}

    pattern = re.compile(r"^LPM\^(\d{3})_MPRAGE_.*\.nii\.gz$")
    files: dict[str, list[Path]] = {}

    for fpath in sorted(MPRAGE_DIR.iterdir()):
        if not fpath.is_file():
            continue
        m = pattern.match(fpath.name)
        if not m:
            log.warning(f"  Unexpected file in MPRAGE dir: {fpath.name}")
            continue
        sub_num = m.group(1)
        files.setdefault(sub_num, []).append(fpath)

    # Resolve subjects with multiple files
    resolved: dict[str, Path] = {}
    for sub_num, paths in files.items():
        if len(paths) == 1:
            resolved[sub_num] = paths[0]
        else:
            # Subject 031: pick the preferred re-scan
            preferred = [p for p in paths if SUBJECT_031_PREFERRED in p.name]
            if preferred:
                resolved[sub_num] = preferred[0]
                skipped = [p for p in paths if p != preferred[0]]
                log.info(
                    f"  sub-{sub_num}: Multiple files, using {preferred[0].name} "
                    f"(skipping {', '.join(p.name for p in skipped)})"
                )
            else:
                # Fallback: use first file
                resolved[sub_num] = paths[0]
                log.warning(
                    f"  sub-{sub_num}: Multiple files, no preferred match. "
                    f"Using {paths[0].name}"
                )

    return resolved


def update_json_sidecar(json_path: Path, dry_run: bool) -> bool:
    """Remove SkullStripped field from JSON sidecar.

    Returns True if the file was modified (or would be in dry-run).
    """
    if not json_path.exists():
        log.warning(f"  JSON sidecar not found: {json_path.name}")
        return False

    with open(json_path, "r") as f:
        data = json.load(f)

    if "SkullStripped" not in data:
        log.info(f"  JSON already clean: {json_path.name}")
        return False

    del data["SkullStripped"]

    if not dry_run:
        with open(json_path, "w") as f:
            json.dump(data, f, indent=2)
            f.write("\n")

    return True


def replace_t1(sub_num: str, source: Path, dry_run: bool) -> bool:
    """Replace a single subject's T1w with the MPRAGE source.

    Returns True on success.
    """
    target_dir = BIDS_ROOT / f"sub-{sub_num}" / "anat"
    target_nii = target_dir / f"sub-{sub_num}_T1w.nii.gz"
    target_json = target_dir / f"sub-{sub_num}_T1w.json"

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
        shutil.copy2(str(source), str(target_nii))

    # Update JSON sidecar
    json_updated = update_json_sidecar(target_json, dry_run)
    if json_updated:
        log.info(f"  sub-{sub_num}: Removed SkullStripped from JSON sidecar")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Replace skull-stripped LPM T1w images with original MPRAGE scans."
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Preview changes without modifying files."
    )
    args = parser.parse_args()

    mode = "DRY RUN" if args.dry_run else "EXECUTE"
    log.info(f"{'='*60}")
    log.info(f"LPM T1w Replacement [{mode}]")
    log.info(f"BIDS root: {BIDS_ROOT}")
    log.info(f"MPRAGE source: {MPRAGE_DIR}")
    log.info(f"{'='*60}")

    # Discover source files
    mprage_files = discover_mprage_files()
    if not mprage_files:
        log.error("No MPRAGE files found. Aborting.")
        return

    log.info(f"Found {len(mprage_files)} MPRAGE files")

    # Replace each subject's T1w
    replaced = 0
    skipped = 0
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
    if missing:
        log.warning(f"BIDS subjects without MPRAGE source: {', '.join(f'sub-{s}' for s in missing)}")
        skipped = len(missing)

    # Summary
    log.info(f"\n{'='*60}")
    log.info(f"SUMMARY [{mode}]")
    log.info(f"  Subjects replaced: {replaced}")
    log.info(f"  Subjects skipped:  {skipped}")
    log.info(f"  Errors:            {errors}")
    log.info(f"{'='*60}")

    if args.dry_run:
        log.info("\nThis was a DRY RUN. Re-run without --dry-run to apply changes.")


if __name__ == "__main__":
    main()
