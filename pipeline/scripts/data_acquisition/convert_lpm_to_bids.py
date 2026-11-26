#!/usr/bin/env python3
"""
Convert Basel LPM raw DICOM data to BIDS format for fMRIPrep processing.

The LPM study (LSD, Psilocybin, Mescaline, Placebo) has raw DICOM files
that must be converted with dcm2niix before organizing into BIDS.

NOTE: Only skull-stripped T1s are available. fMRIPrep should use
      --skull-strip-t1w skip for this dataset.

Source structure:
    data/basel_lpm/LPM_fMRI/LPM_XXX/{Condition}/
        EPI_REST_3_5MM_ISO_NNNN/   -> EPI DICOMs (300 per run)
        B0_FIELD_MAP_NNNN/          -> B0 magnitude DICOMs
        B0_FIELD_MAP_NNNN+1/        -> B0 phasediff DICOMs
    data/basel_lpm/LPM_fMRI/LPM_XXX/T1_skull_stripped.nii.gz

BIDS output:
    data/basel_lpm_bids/
        sub-XXX/anat/sub-XXX_T1w.nii.gz
        sub-XXX/ses-{cond}/func/sub-XXX_ses-{cond}_task-rest[_run-N]_bold.nii.gz
        sub-XXX/ses-{cond}/fmap/sub-XXX_ses-{cond}_magnitude1.nii.gz
        sub-XXX/ses-{cond}/fmap/sub-XXX_ses-{cond}_magnitude2.nii.gz
        sub-XXX/ses-{cond}/fmap/sub-XXX_ses-{cond}_phasediff.nii.gz

Usage:
    python convert_lpm_to_bids.py --dry-run
    python convert_lpm_to_bids.py --subjects LPM_001
    python convert_lpm_to_bids.py
"""

import os
import sys
import json
import shutil
import re
import argparse
import csv
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List

import nibabel as nib
import numpy as np


# =============================================================================
# Constants
# =============================================================================

SOURCE_DIR = Path("D:/projects/TMS_GWT/data/basel_lpm/LPM_fMRI")
OUTPUT_DIR = Path("D:/projects/TMS_GWT/data/basel_lpm_bids")

# Condition folder name -> BIDS session mapping
# Handles typos and alternative naming found in the raw data
CONDITION_MAP = {
    "LSD": "ses-lsd",
    "Psilocybin": "ses-psilocybin",
    "Mescaline": "ses-mescaline",
    "Mecaline": "ses-mescaline",       # Typo in LPM_030
    "Placebo": "ses-placebo",
}

# Folders that should be treated as LSD condition
LSD_ALIASES = {"LSD"}

# B0 fieldmap parameters (from MRI protocol - same as LAM)
B0_PHASEDIFF_PARAMS = {
    "EchoTime1": 0.00492,
    "EchoTime2": 0.00738,
}

# Minimum DICOM count to consider a valid EPI run (300 vols * 35 slices = 10500,
# but DICOMs may be organized differently - use volume count as proxy)
MIN_DICOM_COUNT_FOR_VALID_RUN = 100


# =============================================================================
# Helper functions
# =============================================================================

def check_dcm2niix():
    """Check if dcm2niix is available."""
    try:
        result = subprocess.run(["dcm2niix", "-v"], capture_output=True, text=True)
        return True
    except FileNotFoundError:
        return False


def run_dcm2niix(dicom_dir: Path, output_dir: Path, prefix: str = "out") -> List[Path]:
    """
    Run dcm2niix on a DICOM directory.

    Returns list of generated NIfTI files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "dcm2niix",
        "-z", "y",       # Compress output (.nii.gz)
        "-b", "y",       # Generate BIDS JSON sidecars
        "-f", prefix,    # Output filename prefix
        "-o", str(output_dir),
        str(dicom_dir)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    WARNING: dcm2niix returned {result.returncode}")
        if result.stderr:
            print(f"    STDERR: {result.stderr[:200]}")

    # Collect generated NIfTI files
    niftis = sorted(output_dir.glob(f"{prefix}*.nii.gz"))
    return niftis


def classify_condition_folder(folder_name: str) -> Optional[str]:
    """
    Map a condition folder name to a standard condition name.

    Handles typos and alternative naming conventions.
    """
    # Direct match
    if folder_name in CONDITION_MAP:
        return folder_name

    # LPM_020 has a weird folder name for LSD
    if folder_name.startswith("MRI_SCH"):
        return "LSD"

    return None


def find_epi_dirs(condition_dir: Path) -> List[Path]:
    """Find EPI DICOM directories within a condition folder, sorted by series number."""
    epi_dirs = []
    for d in sorted(condition_dir.iterdir()):
        if d.is_dir() and "EPI" in d.name.upper():
            epi_dirs.append(d)
    return epi_dirs


def find_b0_dirs(condition_dir: Path) -> List[Path]:
    """Find B0 fieldmap DICOM directories."""
    b0_dirs = []
    for d in sorted(condition_dir.iterdir()):
        if d.is_dir() and "B0" in d.name.upper():
            b0_dirs.append(d)
    return b0_dirs


def count_dicoms(directory: Path) -> int:
    """Count DICOM files in a directory."""
    return len([f for f in directory.iterdir()
                if f.is_file() and (f.suffix == ".IMA" or f.suffix == ".dcm"
                                    or f.suffix == "")])


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
    Convert a single LPM subject to BIDS format.

    Parameters
    ----------
    subj_folder : str
        Subject folder name (e.g., 'LPM_001')
    dry_run : bool
        If True, only print what would be done

    Returns
    -------
    dict with conversion status
    """
    match = re.match(r"LPM_(\d+)", subj_folder)
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

    # Create a temp directory for dcm2niix output
    tmp_base = OUTPUT_DIR / ".tmp_dcm2niix" / bids_subj
    if not dry_run:
        tmp_base.mkdir(parents=True, exist_ok=True)

    # ---- T1w anatomical (skull-stripped) ----
    t1_source = subj_source / "T1_skull_stripped.nii.gz"
    if not t1_source.exists():
        # Some subjects have .nii instead of .nii.gz
        t1_source = subj_source / "T1_skull_stripped.nii"

    t1_dest = OUTPUT_DIR / bids_subj / "anat" / f"{bids_subj}_T1w.nii.gz"
    t1_json_dest = t1_dest.with_suffix("").with_suffix(".json")

    if t1_source.exists():
        if dry_run:
            print(f"  [DRY RUN] {t1_source.name} -> {t1_dest.relative_to(OUTPUT_DIR)}")
        else:
            t1_dest.parent.mkdir(parents=True, exist_ok=True)
            if not t1_dest.exists():
                if t1_source.name.endswith(".nii.gz"):
                    shutil.copy2(t1_source, t1_dest)
                else:
                    # Compress .nii to .nii.gz
                    import gzip
                    with open(t1_source, "rb") as f_in:
                        with gzip.open(t1_dest, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
            write_json(t1_json_dest, {
                "SkullStripped": True,
                "Manufacturer": "Siemens",
                "ManufacturersModelName": "Prisma",
                "MagneticFieldStrength": 3,
            })
        result["files_copied"] += 1
    else:
        result["issues"].append("Missing T1_skull_stripped")

    # ---- Process each condition ----
    for cond_dir in sorted(subj_source.iterdir()):
        if not cond_dir.is_dir():
            continue

        condition = classify_condition_folder(cond_dir.name)
        if condition is None:
            continue

        if condition not in CONDITION_MAP:
            continue

        ses = CONDITION_MAP[condition]
        result["sessions"].append(ses)

        # --- Functional EPI ---
        epi_dirs = find_epi_dirs(cond_dir)
        if not epi_dirs:
            result["issues"].append(f"No EPI folders in {cond_dir.name}")
            continue

        # Filter to valid runs (enough DICOMs)
        valid_epis = []
        for ed in epi_dirs:
            n = count_dicoms(ed)
            if n >= MIN_DICOM_COUNT_FOR_VALID_RUN:
                valid_epis.append(ed)
            else:
                if dry_run:
                    print(f"  [DRY RUN] Skipping {ed.name} ({n} DICOMs - too few)")

        if not valid_epis:
            result["issues"].append(f"No valid EPI runs in {cond_dir.name}")
            continue

        # Determine if we need run labels
        use_runs = len(valid_epis) > 1

        for run_idx, epi_dir in enumerate(valid_epis, 1):
            run_suffix = f"_run-{run_idx}" if use_runs else ""

            bold_dest = (OUTPUT_DIR / bids_subj / ses / "func" /
                         f"{bids_subj}_{ses}_task-rest{run_suffix}_bold.nii.gz")
            bold_json_dest = bold_dest.with_suffix("").with_suffix(".json")

            if dry_run:
                print(f"  [DRY RUN] dcm2niix {epi_dir.name} -> "
                      f"{bold_dest.relative_to(OUTPUT_DIR)}")
                result["files_copied"] += 1
                continue

            if bold_dest.exists():
                result["files_copied"] += 1
                continue

            # Run dcm2niix
            tmp_out = tmp_base / f"{ses}_func{run_suffix}"
            niftis = run_dcm2niix(epi_dir, tmp_out, prefix="epi")

            if niftis:
                bold_dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(niftis[0]), str(bold_dest))

                # Handle JSON sidecar
                json_src = niftis[0].with_suffix("").with_suffix(".json")
                if json_src.exists():
                    json_data = load_json(json_src)
                    json_data["TaskName"] = "rest"
                    write_json(bold_json_dest, json_data)
                else:
                    # Create minimal JSON
                    write_json(bold_json_dest, {
                        "TaskName": "rest",
                        "RepetitionTime": 1.8,
                        "EchoTime": 0.028,
                        "Manufacturer": "Siemens",
                        "ManufacturersModelName": "Prisma",
                        "MagneticFieldStrength": 3,
                    })

                result["files_copied"] += 1
            else:
                result["issues"].append(
                    f"dcm2niix produced no output for {epi_dir.name}")

        # --- B0 Fieldmaps ---
        b0_dirs = find_b0_dirs(cond_dir)
        if len(b0_dirs) < 2:
            result["issues"].append(
                f"Expected 2 B0 folders in {cond_dir.name}, found {len(b0_dirs)}")
            continue

        # First B0 dir = magnitude, second = phasediff
        mag_dicom_dir = b0_dirs[0]
        phase_dicom_dir = b0_dirs[1]

        # -- Magnitude --
        mag1_dest = (OUTPUT_DIR / bids_subj / ses / "fmap" /
                     f"{bids_subj}_{ses}_magnitude1.nii.gz")
        mag2_dest = (OUTPUT_DIR / bids_subj / ses / "fmap" /
                     f"{bids_subj}_{ses}_magnitude2.nii.gz")

        if dry_run:
            print(f"  [DRY RUN] dcm2niix {mag_dicom_dir.name} -> magnitude1 + magnitude2")
        elif not mag1_dest.exists():
            tmp_out = tmp_base / f"{ses}_mag"
            niftis = run_dcm2niix(mag_dicom_dir, tmp_out, prefix="mag")

            mag1_dest.parent.mkdir(parents=True, exist_ok=True)
            if niftis:
                img = nib.load(niftis[0])
                data = img.get_fdata()
                if data.ndim == 4 and data.shape[3] >= 2:
                    for i, dst in enumerate([mag1_dest, mag2_dest]):
                        vol = nib.Nifti1Image(
                            data[:, :, :, i], img.affine, img.header)
                        nib.save(vol, dst)
                elif data.ndim == 3:
                    shutil.move(str(niftis[0]), str(mag1_dest))
                    if len(niftis) > 1:
                        shutil.move(str(niftis[1]), str(mag2_dest))
                    else:
                        print(f"    WARNING: Only 1 magnitude volume for "
                              f"{bids_subj}/{ses}")
                # If dcm2niix produced 2 separate files (echo split)
                elif len(niftis) >= 2:
                    shutil.move(str(niftis[0]), str(mag1_dest))
                    shutil.move(str(niftis[1]), str(mag2_dest))
        result["files_copied"] += 2

        # -- Phasediff --
        phase_dest = (OUTPUT_DIR / bids_subj / ses / "fmap" /
                      f"{bids_subj}_{ses}_phasediff.nii.gz")
        phase_json_dest = phase_dest.with_suffix("").with_suffix(".json")

        # Determine IntendedFor (first/only run)
        if use_runs:
            intended_bold = (f"bids::{bids_subj}/{ses}/func/"
                             f"{bids_subj}_{ses}_task-rest_run-1_bold.nii.gz")
        else:
            intended_bold = (f"bids::{bids_subj}/{ses}/func/"
                             f"{bids_subj}_{ses}_task-rest_bold.nii.gz")

        if dry_run:
            print(f"  [DRY RUN] dcm2niix {phase_dicom_dir.name} -> phasediff")
        elif not phase_dest.exists():
            tmp_out = tmp_base / f"{ses}_phase"
            niftis = run_dcm2niix(phase_dicom_dir, tmp_out, prefix="phase")

            phase_dest.parent.mkdir(parents=True, exist_ok=True)
            if niftis:
                shutil.move(str(niftis[0]), str(phase_dest))

        # Phasediff JSON
        phase_json = {
            **B0_PHASEDIFF_PARAMS,
            "IntendedFor": intended_bold,
        }
        write_json(phase_json_dest, phase_json, dry_run)
        result["files_copied"] += 1

    # Cleanup temp directory
    if not dry_run and tmp_base.exists():
        shutil.rmtree(tmp_base, ignore_errors=True)

    if result["issues"]:
        result["status"] = "partial"

    return result


# =============================================================================
# Dataset-level BIDS files
# =============================================================================

def create_dataset_description(output_dir: Path):
    """Create BIDS dataset_description.json."""
    desc = {
        "Name": "Basel LPM Study (LSD, Psilocybin, Mescaline, Placebo)",
        "BIDSVersion": "1.8.0",
        "DatasetType": "raw",
        "License": "Restricted",
        "Authors": ["Borgwardt, S.", "Preller, K.H.", "Vollenweider, F.X."],
        "GeneratedBy": [{
            "Name": "convert_lpm_to_bids.py",
            "Description": "LPM raw DICOM to BIDS conversion for GWT pipeline"
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
    """Find all subject folders in the LPM source directory."""
    subjects = []
    for d in sorted(SOURCE_DIR.iterdir()):
        if d.is_dir() and d.name.startswith("LPM_"):
            subjects.append(d.name)
    return subjects


def main():
    parser = argparse.ArgumentParser(
        description="Convert Basel LPM raw DICOM data to BIDS format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python convert_lpm_to_bids.py --dry-run
  python convert_lpm_to_bids.py --subjects LPM_001
  python convert_lpm_to_bids.py

Requirements:
  dcm2niix must be installed and in PATH.
  Install via: conda install -c conda-forge dcm2niix
        """
    )
    parser.add_argument("--subjects", nargs="+",
                        help="Specific subject folder names (default: all)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Preview without converting")

    args = parser.parse_args()

    print("=" * 60)
    print("Basel LPM -> BIDS Conversion (DICOM via dcm2niix)")
    print("=" * 60)
    print(f"Source: {SOURCE_DIR}")
    print(f"Output: {OUTPUT_DIR}")
    print()

    # Check dcm2niix
    if not args.dry_run and not check_dcm2niix():
        print("ERROR: dcm2niix not found in PATH.")
        print("Install with: conda install -c conda-forge dcm2niix")
        sys.exit(1)

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

    # Cleanup global temp dir
    tmp_dir = OUTPUT_DIR / ".tmp_dcm2niix"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir, ignore_errors=True)

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
    print(f"  Total files {'to process' if args.dry_run else 'created'}: {total_files}")

    if partial > 0:
        print("\nSubjects with issues:")
        for r in results:
            if r["status"] == "partial":
                print(f"  sub-{r['subject']}: {'; '.join(r['issues'])}")

    print()
    if args.dry_run:
        print("DRY RUN complete. No files were created.")
    else:
        print(f"BIDS dataset created at: {OUTPUT_DIR}")
        print("NOTE: This dataset uses skull-stripped T1s.")
        print("      Use fMRIPrep with: --skull-strip-t1w skip")


if __name__ == "__main__":
    main()
