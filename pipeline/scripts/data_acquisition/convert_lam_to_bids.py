#!/usr/bin/env python3
"""
Convert Basel LAM raw NIfTI data to BIDS format for fMRIPrep processing.

The LAM study (LSD, Amphetamine, MDMA, Placebo) has raw NIfTI files with
JSON sidecars available for LSD and Placebo conditions only.

Source structure:
    data/basel_lam/LAM fMRI/lam_fmri_XXX/{Condition}/
        NNN_EPI_Rest_3_5mm_iso_{Condition}/  -> functional .nii
        NNN_B0_field_map_{Condition}/         -> magnitude .nii (4D, 2 echoes)
        NNN+1_B0_field_map_{Condition}/       -> phasediff .nii (3D)
    data/basel_lam/LAM fMRI/lam_fmri_XXX/T1.nii.gz  -> anatomical
    data/basel_lam/LAM fMRI/LAM_json/                -> JSON sidecars (LSD+Placebo)

BIDS output:
    data/basel_lam_bids/
        sub-XXX/anat/sub-XXX_T1w.nii.gz
        sub-XXX/ses-{cond}/func/sub-XXX_ses-{cond}_task-rest_bold.nii.gz
        sub-XXX/ses-{cond}/fmap/sub-XXX_ses-{cond}_magnitude1.nii.gz
        sub-XXX/ses-{cond}/fmap/sub-XXX_ses-{cond}_magnitude2.nii.gz
        sub-XXX/ses-{cond}/fmap/sub-XXX_ses-{cond}_phasediff.nii.gz

Usage:
    python convert_lam_to_bids.py
    python convert_lam_to_bids.py --dry-run
    python convert_lam_to_bids.py --subjects lam_fmri_001 lam_fmri_002
"""

import os
import sys
import json
import gzip
import shutil
import re
import argparse
import csv
from pathlib import Path
from typing import Optional, List, Tuple

import nibabel as nib
import numpy as np


# =============================================================================
# Constants
# =============================================================================

SOURCE_DIR = Path("D:/projects/TMS_GWT/data/basel_lam/LAM fMRI")
JSON_DIR = SOURCE_DIR / "LAM_json"
OUTPUT_DIR = Path("D:/projects/TMS_GWT/data/basel_lam_bids")

# Condition -> BIDS session mapping
CONDITION_MAP = {
    "LSD": "ses-lsd",
    "Amph": "ses-amph",
    "MDMA": "ses-mdma",
    "Placebo": "ses-placebo",
}

# Valid conditions (skip dropout sessions)
VALID_CONDITIONS = set(CONDITION_MAP.keys())

# B0 fieldmap parameters (from MRI protocol)
B0_PHASEDIFF_PARAMS = {
    "EchoTime1": 0.00492,
    "EchoTime2": 0.00738,
}


# =============================================================================
# Helper functions
# =============================================================================

def find_json_sidecar(subject_num: str, condition: str) -> Optional[Path]:
    """
    Find the JSON sidecar for a given subject and condition.

    JSON files are named like:
        LAM-001-LSD_EPI_Rest_3.5mm_iso_YYYYMMDDHHMMSS.json
        LAM-029_-_LSD_EPI_Rest_3.5mm_iso_YYYYMMDDHHMMSS.json  (special case)

    Only LSD and Placebo JSONs exist. For Amph/MDMA, returns None.
    """
    if condition not in ("LSD", "Placebo"):
        return None

    # Try standard pattern first
    pattern = re.compile(
        rf"LAM-0*{subject_num}[-_].*{condition}_EPI_Rest_3\.5mm_iso_\d+\.json"
    )
    for f in JSON_DIR.iterdir():
        if pattern.match(f.name):
            return f

    return None


def find_epi_folder(condition_dir: Path) -> Optional[Path]:
    """Find the EPI functional data folder within a condition directory."""
    for d in sorted(condition_dir.iterdir()):
        if d.is_dir() and "EPI_Rest" in d.name:
            return d
    return None


def find_b0_folders(condition_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Find B0 fieldmap folders. Returns (magnitude_dir, phasediff_dir).

    Convention: the first B0 folder is magnitude (4D), second is phasediff (3D).
    """
    b0_dirs = sorted([
        d for d in condition_dir.iterdir()
        if d.is_dir() and "B0_field_map" in d.name
    ])
    if len(b0_dirs) >= 2:
        return b0_dirs[0], b0_dirs[1]
    elif len(b0_dirs) == 1:
        return b0_dirs[0], None
    return None, None


def find_nifti_in_dir(directory: Path) -> Optional[Path]:
    """Find the .nii or .nii.gz file in a directory."""
    for f in directory.iterdir():
        if f.name.endswith(".nii") or f.name.endswith(".nii.gz"):
            return f
    return None


def compress_nifti(src: Path, dst: Path, dry_run: bool = False):
    """Copy a NIfTI file, compressing if needed (.nii -> .nii.gz)."""
    if dry_run:
        return

    dst.parent.mkdir(parents=True, exist_ok=True)

    if src.name.endswith(".nii.gz"):
        shutil.copy2(src, dst)
    else:
        # Compress .nii to .nii.gz
        with open(src, "rb") as f_in:
            with gzip.open(dst, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)


def split_4d_magnitude(src: Path, dst1: Path, dst2: Path, dry_run: bool = False):
    """Split a 4D magnitude NIfTI (x,y,z,2) into magnitude1 and magnitude2."""
    if dry_run:
        return

    dst1.parent.mkdir(parents=True, exist_ok=True)

    img = nib.load(src)
    data = img.get_fdata()

    if data.ndim == 4 and data.shape[3] >= 2:
        # Split into two 3D volumes
        for i, dst in enumerate([dst1, dst2]):
            vol = nib.Nifti1Image(data[:, :, :, i], img.affine, img.header)
            nib.save(vol, dst)
    elif data.ndim == 3:
        # Already 3D - just copy as magnitude1, skip magnitude2
        nib.save(img, dst1)
        print(f"    WARNING: B0 magnitude is 3D, saved as magnitude1 only")
    else:
        print(f"    WARNING: Unexpected B0 magnitude shape: {data.shape}")


def write_json(filepath: Path, data: dict, dry_run: bool = False):
    """Write a JSON sidecar file."""
    if dry_run:
        return
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(data, f, indent=2)


def load_json(filepath: Path) -> dict:
    """Load a JSON file."""
    with open(filepath, "r") as f:
        return json.load(f)


# =============================================================================
# Subject conversion
# =============================================================================

def convert_subject(subj_folder: str, dry_run: bool = False) -> dict:
    """
    Convert a single LAM subject to BIDS format.

    Parameters
    ----------
    subj_folder : str
        Subject folder name (e.g., 'lam_fmri_001')
    dry_run : bool
        If True, only print what would be done

    Returns
    -------
    dict with conversion status
    """
    # Extract subject number
    match = re.match(r"lam_fmri_(\d+)", subj_folder)
    if not match:
        return {"subject": subj_folder, "status": "error",
                "issues": [f"Cannot parse subject number from '{subj_folder}'"]}

    subj_num = match.group(1)  # e.g., "001"
    bids_subj = f"sub-{subj_num}"
    subj_source = SOURCE_DIR / subj_folder

    result = {
        "subject": subj_num, "status": "success",
        "files_copied": 0, "issues": [], "sessions": []
    }

    # ---- T1w anatomical (one per subject, outside sessions) ----
    t1_source = subj_source / "T1.nii.gz"
    t1_dest = OUTPUT_DIR / bids_subj / "anat" / f"{bids_subj}_T1w.nii.gz"
    t1_json_dest = OUTPUT_DIR / bids_subj / "anat" / f"{bids_subj}_T1w.json"

    if t1_source.exists():
        if dry_run:
            print(f"  [DRY RUN] T1.nii.gz -> {t1_dest.relative_to(OUTPUT_DIR)}")
        else:
            t1_dest.parent.mkdir(parents=True, exist_ok=True)
            if not t1_dest.exists():
                shutil.copy2(t1_source, t1_dest)
            # Minimal T1 JSON
            write_json(t1_json_dest, {
                "Manufacturer": "Siemens",
                "ManufacturersModelName": "Prisma",
                "MagneticFieldStrength": 3,
            })
        result["files_copied"] += 1
    else:
        result["issues"].append("Missing T1.nii.gz")

    # ---- Process each condition ----
    # Store a reference JSON from LSD/Placebo to clone for Amph/MDMA
    reference_json_data = None

    for condition in ["LSD", "Placebo", "Amph", "MDMA"]:
        condition_dir = subj_source / condition
        if not condition_dir.exists():
            result["issues"].append(f"Missing condition: {condition}")
            continue

        ses = CONDITION_MAP[condition]
        result["sessions"].append(ses)

        # --- Functional EPI ---
        epi_dir = find_epi_folder(condition_dir)
        if epi_dir is None:
            result["issues"].append(f"No EPI folder in {condition}")
            continue

        epi_nifti = find_nifti_in_dir(epi_dir)
        if epi_nifti is None:
            result["issues"].append(f"No NIfTI in EPI folder for {condition}")
            continue

        bold_dest = (OUTPUT_DIR / bids_subj / ses / "func" /
                     f"{bids_subj}_{ses}_task-rest_bold.nii.gz")
        bold_json_dest = bold_dest.with_suffix("").with_suffix(".json")

        if dry_run:
            print(f"  [DRY RUN] {epi_nifti.name} -> {bold_dest.relative_to(OUTPUT_DIR)}")
        else:
            if not bold_dest.exists():
                compress_nifti(epi_nifti, bold_dest)
        result["files_copied"] += 1

        # --- BOLD JSON sidecar ---
        json_sidecar = find_json_sidecar(subj_num, condition)
        if json_sidecar is not None:
            # Use the real JSON
            json_data = load_json(json_sidecar)
            json_data["TaskName"] = "rest"
            reference_json_data = json_data  # Save for Amph/MDMA
            write_json(bold_json_dest, json_data, dry_run)
            if dry_run:
                print(f"  [DRY RUN] JSON: {json_sidecar.name}")
        elif reference_json_data is not None:
            # Clone from LSD/Placebo JSON (identical protocol)
            cloned = dict(reference_json_data)
            cloned["TaskName"] = "rest"
            write_json(bold_json_dest, cloned, dry_run)
            if dry_run:
                print(f"  [DRY RUN] JSON: cloned from reference for {condition}")
        else:
            result["issues"].append(f"No JSON sidecar for {condition} (and no reference)")

        # --- B0 Fieldmaps ---
        mag_dir, phase_dir = find_b0_folders(condition_dir)

        if mag_dir is not None:
            mag_nifti = find_nifti_in_dir(mag_dir)
            if mag_nifti is not None:
                mag1_dest = (OUTPUT_DIR / bids_subj / ses / "fmap" /
                             f"{bids_subj}_{ses}_magnitude1.nii.gz")
                mag2_dest = (OUTPUT_DIR / bids_subj / ses / "fmap" /
                             f"{bids_subj}_{ses}_magnitude2.nii.gz")

                if dry_run:
                    print(f"  [DRY RUN] B0 mag -> magnitude1 + magnitude2")
                else:
                    if not mag1_dest.exists():
                        split_4d_magnitude(mag_nifti, mag1_dest, mag2_dest)
                result["files_copied"] += 2

        if phase_dir is not None:
            phase_nifti = find_nifti_in_dir(phase_dir)
            if phase_nifti is not None:
                phase_dest = (OUTPUT_DIR / bids_subj / ses / "fmap" /
                              f"{bids_subj}_{ses}_phasediff.nii.gz")
                phase_json_dest = phase_dest.with_suffix("").with_suffix(".json")

                if dry_run:
                    print(f"  [DRY RUN] B0 phase -> phasediff")
                else:
                    if not phase_dest.exists():
                        compress_nifti(phase_nifti, phase_dest)

                # Phasediff JSON with IntendedFor
                intended_for = (f"bids::{bids_subj}/{ses}/func/"
                                f"{bids_subj}_{ses}_task-rest_bold.nii.gz")
                phase_json = {
                    **B0_PHASEDIFF_PARAMS,
                    "IntendedFor": intended_for,
                }
                write_json(phase_json_dest, phase_json, dry_run)
                result["files_copied"] += 1

    if result["issues"]:
        result["status"] = "partial"

    return result


# =============================================================================
# Dataset-level BIDS files
# =============================================================================

def create_dataset_description(output_dir: Path):
    """Create BIDS dataset_description.json."""
    desc = {
        "Name": "Basel LAM Study (LSD, Amphetamine, MDMA, Placebo)",
        "BIDSVersion": "1.8.0",
        "DatasetType": "raw",
        "License": "Restricted",
        "Authors": ["Borgwardt, S.", "Preller, K.H.", "Vollenweider, F.X."],
        "GeneratedBy": [{
            "Name": "convert_lam_to_bids.py",
            "Description": "LAM raw NIfTI to BIDS conversion for GWT pipeline"
        }]
    }
    with open(output_dir / "dataset_description.json", "w") as f:
        json.dump(desc, f, indent=2)


def create_participants_tsv(output_dir: Path, subjects: list):
    """Create BIDS participants.tsv."""
    with open(output_dir / "participants.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["participant_id"])
        for subj in sorted(subjects):
            writer.writerow([f"sub-{subj}"])


# =============================================================================
# Main
# =============================================================================

def discover_subjects() -> List[str]:
    """Find all subject folders in the LAM source directory."""
    subjects = []
    for d in sorted(SOURCE_DIR.iterdir()):
        if d.is_dir() and d.name.startswith("lam_fmri_"):
            subjects.append(d.name)
    return subjects


def main():
    parser = argparse.ArgumentParser(
        description="Convert Basel LAM raw NIfTI data to BIDS format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python convert_lam_to_bids.py --dry-run
  python convert_lam_to_bids.py --subjects lam_fmri_001
  python convert_lam_to_bids.py
        """
    )
    parser.add_argument("--subjects", nargs="+",
                        help="Specific subject folder names (default: all)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview without copying files")

    args = parser.parse_args()

    print("=" * 60)
    print("Basel LAM -> BIDS Conversion")
    print("=" * 60)
    print(f"Source: {SOURCE_DIR}")
    print(f"Output: {OUTPUT_DIR}")
    print()

    # Discover subjects
    if args.subjects:
        subjects = args.subjects
    else:
        subjects = discover_subjects()

    if not subjects:
        print("ERROR: No subjects found.")
        sys.exit(1)

    print(f"Found {len(subjects)} subjects")
    print()

    if not args.dry_run:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        create_dataset_description(OUTPUT_DIR)

    # Convert each subject
    results = []
    for i, subj in enumerate(subjects, 1):
        print(f"[{i}/{len(subjects)}] {subj}...")
        result = convert_subject(subj, dry_run=args.dry_run)
        results.append(result)

        if result["issues"]:
            for issue in result["issues"]:
                print(f"  WARNING: {issue}")

    # Create participants.tsv
    successful = [r["subject"] for r in results
                  if r["status"] in ("success", "partial")]
    if not args.dry_run and successful:
        create_participants_tsv(OUTPUT_DIR, successful)

    # Summary
    print()
    print("=" * 60)
    print("CONVERSION SUMMARY")
    print("=" * 60)
    success = sum(1 for r in results if r["status"] == "success")
    partial = sum(1 for r in results if r["status"] == "partial")
    errors = sum(1 for r in results if r["status"] == "error")
    total_files = sum(r["files_copied"] for r in results)
    print(f"  Successful: {success}")
    print(f"  Partial: {partial}")
    print(f"  Errors: {errors}")
    print(f"  Total files {'to copy' if args.dry_run else 'copied'}: {total_files}")

    if partial > 0:
        print("\nSubjects with issues:")
        for r in results:
            if r["status"] == "partial":
                print(f"  sub-{r['subject']}: {'; '.join(r['issues'])}")

    print()
    if args.dry_run:
        print("DRY RUN complete. No files were copied.")
    else:
        print(f"BIDS dataset created at: {OUTPUT_DIR}")
        print("Next step: Run bids-validator, then fMRIPrep")


if __name__ == "__main__":
    main()
