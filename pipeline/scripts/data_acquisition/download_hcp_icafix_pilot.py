#!/usr/bin/env python3
"""
Download HCP minimally preprocessed (ICA-FIX cleaned) data for the
20-subject pilot — exact match to our existing 100-subject HCP cohort
(first 20 of `hcp_rest_workspace_metadata.yaml`).

Why: Luppi 2024 SAPHIRE feeds CONN the HCP-pipeline-cleaned (ICA-FIX) volumes,
not raw fMRI. Our pipeline currently feeds CONN fMRIPrep-aligned raw volumes.
This pilot lets us test whether the synergy reproduction gap (r=0.37 vs Luppi's
matrices, redundancy r=0.90) closes when we match Luppi's input more carefully.

Per subject (~1 GB), files pulled:
  - MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_hp2000_clean.nii.gz
      The ICA-FIX cleaned BOLD in MNI152NLin6Asym 2mm.
  - MNINonLinear/Results/rfMRI_REST1_LR/Movement_Regressors.txt
      12-column motion file (first 6 = translations + rotations).
  - MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_SBRef.nii.gz
      Single-band reference (small, useful for QC).
  - MNINonLinear/T1w_restore_brain.nii.gz
      Skull-stripped MNI-aligned T1 (input for CONN structural_segment).
  - MNINonLinear/brainmask_fs.nii.gz
      Brain mask aligned to BOLD.

Total disk: ~20 GB for 20 subjects.

PREREQUISITES
=============
HCP S3 access is NOT anonymous. You need:
  1) A registered ConnectomeDB account (https://db.humanconnectome.org/)
  2) Accepted "WU-Minn HCP Data - 1200 Subjects, ..." Open Access Data Use Terms.
  3) Personal AWS keys generated via ConnectomeDB > Activate Amazon Bucket.
  4) Configured locally either as env vars
       set AWS_ACCESS_KEY_ID=...
       set AWS_SECRET_ACCESS_KEY=...
     or as a named profile (~/.aws/credentials):
       [hcp]
       aws_access_key_id = ...
       aws_secret_access_key = ...
     and pass --profile hcp.

USAGE
=====
  python download_hcp_icafix_pilot.py                          # default 20-subj pilot
  python download_hcp_icafix_pilot.py --profile hcp            # use named AWS profile
  python download_hcp_icafix_pilot.py --subjects 100307 100408 # explicit subset
  python download_hcp_icafix_pilot.py --dry-run                # print plan, don't download

Output layout (BIDS-derivative-flavour, kept separate from raw HCP):

  data/hcp_icafix/
      sub-{subject}/
          func/
              sub-{subject}_task-rest_dir-LR_desc-icafix_bold.nii.gz
              sub-{subject}_task-rest_dir-LR_desc-icafix_motion.txt
              sub-{subject}_task-rest_dir-LR_desc-icafix_sbref.nii.gz
          anat/
              sub-{subject}_T1w_brain.nii.gz
              sub-{subject}_brainmask.nii.gz
      participants.tsv
      manifest.json
"""

import argparse
import json
import logging
import subprocess
import sys
import time
from pathlib import Path

S3_BUCKET = "s3://hcp-openaccess/HCP_1200"
LOCAL_ROOT = Path("D:/projects/TMS_GWT/data/hcp_icafix")

# First 20 of the 100-subject HCP cohort already used by the existing pipeline
# (`results/group_analysis/hcp/rest/hcp_rest_workspace_metadata.yaml`).
PILOT_SUBJECTS = [
    "100307", "100408", "101107", "101309", "101915",
    "103111", "103414", "103818", "105014", "105115",
    "106016", "108828", "110411", "111312", "111716",
    "113619", "113922", "114419", "115320", "116524",
]

REST_RUN = "rfMRI_REST1_LR"

# (s3 relative key, local relative path under sub-{subject}/)
PER_SUBJECT_FILES = [
    (
        f"MNINonLinear/Results/{REST_RUN}/{REST_RUN}_hp2000_clean.nii.gz",
        "func/sub-{subject}_task-rest_dir-LR_desc-icafix_bold.nii.gz",
    ),
    (
        f"MNINonLinear/Results/{REST_RUN}/Movement_Regressors.txt",
        "func/sub-{subject}_task-rest_dir-LR_desc-icafix_motion.txt",
    ),
    (
        f"MNINonLinear/Results/{REST_RUN}/{REST_RUN}_SBRef.nii.gz",
        "func/sub-{subject}_task-rest_dir-LR_desc-icafix_sbref.nii.gz",
    ),
    (
        "MNINonLinear/T1w_restore_brain.nii.gz",
        "anat/sub-{subject}_T1w_brain.nii.gz",
    ),
    (
        "MNINonLinear/brainmask_fs.nii.gz",
        "anat/sub-{subject}_brainmask.nii.gz",
    ),
]


def setup_logging(log_dir: Path) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"download_hcp_icafix_{time.strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def aws_check(profile: str | None) -> bool:
    """Quick credential sanity check before launching downloads."""
    cmd = ["aws", "sts", "get-caller-identity"]
    if profile:
        cmd += ["--profile", profile]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    except FileNotFoundError:
        logging.error("`aws` CLI not found. Install AWS CLI v2 first.")
        return False
    except subprocess.TimeoutExpired:
        logging.error("aws sts call timed out")
        return False
    if result.returncode != 0:
        logging.error("AWS credentials not configured or invalid:")
        logging.error(result.stderr.strip())
        return False
    logging.info("AWS credentials OK: %s", result.stdout.strip())
    return True


def s3_cp(s3_url: str, local_path: Path, profile: str | None, dry_run: bool) -> bool:
    """Download single S3 object via aws s3 cp; skip if file already exists."""
    if local_path.exists() and local_path.stat().st_size > 0:
        logging.info("  SKIP (already exists): %s", local_path.name)
        return True

    local_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["aws", "s3", "cp", s3_url, str(local_path)]
    if profile:
        cmd += ["--profile", profile]

    if dry_run:
        logging.info("  DRY-RUN: would run: %s", " ".join(cmd))
        return True

    logging.info("  GET %s -> %s", s3_url, local_path.name)
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error("    FAILED (%s): %s", result.returncode, result.stderr.strip())
        return False
    size_mb = local_path.stat().st_size / 1e6 if local_path.exists() else 0
    logging.info("    OK (%.1f MB, %.1fs)", size_mb, time.time() - t0)
    return True


def download_subject(subject: str, profile: str | None, dry_run: bool) -> dict:
    """Download all needed files for one HCP subject."""
    logging.info("=== sub-%s ===", subject)
    sub_root = LOCAL_ROOT / f"sub-{subject}"
    record = {"subject": subject, "files": {}, "n_ok": 0, "n_fail": 0}

    for s3_rel, local_template in PER_SUBJECT_FILES:
        s3_url = f"{S3_BUCKET}/{subject}/{s3_rel}"
        local_path = sub_root / local_template.format(subject=subject)
        ok = s3_cp(s3_url, local_path, profile, dry_run)
        record["files"][s3_rel] = {
            "local": str(local_path.relative_to(LOCAL_ROOT)),
            "ok": ok,
        }
        if ok:
            record["n_ok"] += 1
        else:
            record["n_fail"] += 1
    return record


def write_participants_tsv(subjects: list[str]) -> None:
    """Minimal participants.tsv, mirrors HCP-1200 'unrelated' subset."""
    p = LOCAL_ROOT / "participants.tsv"
    if p.exists():
        return
    with p.open("w") as f:
        f.write("participant_id\tdataset_source\trun\n")
        for s in subjects:
            f.write(f"sub-{s}\tHCP_1200_minimal_preproc\tREST1_LR\n")
    logging.info("Wrote %s", p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--subjects",
        nargs="+",
        default=PILOT_SUBJECTS,
        help=f"Subject IDs to download. Defaults to first-20 pilot ({len(PILOT_SUBJECTS)} subjects).",
    )
    ap.add_argument("--profile", default=None, help="AWS named profile (e.g. 'hcp').")
    ap.add_argument("--dry-run", action="store_true", help="Print plan only, no downloads.")
    ap.add_argument("--skip-cred-check", action="store_true", help="Skip aws sts check.")
    args = ap.parse_args()

    LOCAL_ROOT.mkdir(parents=True, exist_ok=True)
    log_path = setup_logging(LOCAL_ROOT / "logs")
    logging.info("Pilot HCP ICA-FIX download")
    logging.info("  Subjects (%d): %s", len(args.subjects), ", ".join(args.subjects))
    logging.info("  Output root: %s", LOCAL_ROOT)
    logging.info("  Profile: %s", args.profile or "(default)")
    logging.info("  Dry run: %s", args.dry_run)
    logging.info("  Log: %s", log_path)

    if not args.skip_cred_check and not args.dry_run:
        if not aws_check(args.profile):
            logging.error(
                "Set up HCP S3 keys via ConnectomeDB > Activate Amazon Bucket, then "
                "configure as env vars or aws profile. See script docstring."
            )
            sys.exit(2)

    manifest = {
        "subjects": args.subjects,
        "s3_bucket": S3_BUCKET,
        "rest_run": REST_RUN,
        "files_per_subject": [s3_rel for s3_rel, _ in PER_SUBJECT_FILES],
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "records": [],
    }

    for subject in args.subjects:
        rec = download_subject(subject, args.profile, args.dry_run)
        manifest["records"].append(rec)

    write_participants_tsv(args.subjects)
    manifest["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    manifest_path = LOCAL_ROOT / "manifest.json"
    with manifest_path.open("w") as f:
        json.dump(manifest, f, indent=2)
    logging.info("Manifest written: %s", manifest_path)

    n_ok = sum(r["n_ok"] for r in manifest["records"])
    n_fail = sum(r["n_fail"] for r in manifest["records"])
    expected = len(args.subjects) * len(PER_SUBJECT_FILES)
    logging.info("=" * 60)
    logging.info("Summary: %d OK, %d FAIL (expected %d)", n_ok, n_fail, expected)
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
