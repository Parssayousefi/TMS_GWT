"""Map deposited Aachen targets using independently audited corrected anatomy.

Writes only output/aachen_alignment_review/<label>.  No production files,
timeseries, effect matrices or manuscript results are read or changed.
An anatomy audit from aachen_alignment_anatomical_review.py must exist first.
"""
from argparse import ArgumentParser
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

import nibabel as nib
import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
from aachen_alignment_anatomical_review import (
    ROOT, REVIEW, SESSIONS, array, correlation, fraction, mask_image,
    selected_masks, sitk,
)


def main():
    parser = ArgumentParser()
    parser.add_argument("--derivatives-root", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--subjects", nargs="*")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.label):
        raise ValueError("Label must be a plain directory name")
    destination = REVIEW / args.label
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f"Use a new label; refusing to overwrite {destination}")
    audit_provenance = json.loads((args.audit_dir / "provenance.json").read_text())
    if Path(audit_provenance["derivatives_root"]).resolve() != args.derivatives_root.resolve():
        raise ValueError("Audit belongs to a different derivative directory")
    anatomy = pd.read_csv(args.audit_dir / "anatomical_support.csv", dtype={"subject": str}).set_index("subject")
    settings = pd.read_csv(args.audit_dir / "effective_skullstrip_settings.csv", dtype={"subject": str})
    audited_sessions = pd.read_csv(args.audit_dir / "functional_anatomical_agreement.csv", dtype={"subject": str})
    source = selected_masks()
    subjects = sorted({r["subject"] for r in source})
    if args.subjects:
        requested = {s.zfill(2) for s in args.subjects}
        if not requested <= set(subjects):
            raise ValueError("Requested participant has no usable deposited target")
        subjects = sorted(requested)
    # Fail the whole requested mapping if anatomy has not been corrected; do not
    # quietly create a favorable subset of subjects after a quality failure.
    for subject in subjects:
        if subject not in anatomy.index:
            raise ValueError(f"Missing independent anatomy audit for {subject}")
        if anatomy.loc[subject, "fraction_native_support_retained"] < .95:
            raise ValueError(f"Native anatomy support is still truncated for {subject}")
        modes = set(settings.loc[settings.subject.eq(subject), "skull_strip_t1w"])
        if modes != {"skip"}:
            raise ValueError(f"Expected audited skull-strip skip, found {modes} for {subject}")
        needed = {SESSIONS[r["session_number"]] for r in source if r["subject"] == subject}
        have = set(audited_sessions.loc[audited_sessions.subject.eq(subject), "session"])
        if not needed <= have:
            raise ValueError(f"Missing corrected functional audit for {subject}: {needed - have}")
    for relative, expected in audit_provenance["input_sha256"].items():
        path = ROOT / relative
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Audited input changed: {path}")
    destination.mkdir(parents=True, exist_ok=True)
    hashes = {}

    def track(path):
        path = Path(path)
        hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    atlas_path = track(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    atlas = nib.load(atlas_path)
    labels = np.asarray(atlas.dataobj).astype(int)
    cortical = (labels >= 1) & (labels <= 400)
    reference = sitk.ReadImage(str(atlas_path))
    qc, coverage, weights, transformations = [], [], [], []
    for subject in subjects:
        base = args.derivatives_root / f"sub-{subject}"
        anat = base / "ses-anat/anat"
        prefix = f"sub-{subject}_ses-anat"
        t1 = sitk.ReadImage(str(anat / f"{prefix}_desc-preproc_T1w.nii.gz"))
        t1mask = array(sitk.ReadImage(str(anat / f"{prefix}_desc-brain_mask.nii.gz"))) > 0
        a2m = sitk.ReadTransform(str(track(anat / f"{prefix}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5")))
        m2a = sitk.ReadTransform(str(track(anat / f"{prefix}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5")))
        block = [r for r in source if r["subject"] == subject]
        chosen = sorted({r["session_number"] for r in block})
        masks = {"l": [], "r": []}
        for session in chosen:
            condition = SESSIONS[session]
            stem = f"sub-{subject}_ses-{condition}_task-rest"
            func = base / f"ses-{condition}/func"
            native_path = track(func / f"{stem}_desc-coreg_boldref.nii.gz")
            native = sitk.ReadImage(str(native_path))
            native_nib = nib.load(native_path)
            raw = nib.load(ROOT / f"data/aachen_itbs/sub-{subject}/ses-{condition}/func/{stem}_bold.nii.gz")
            np.testing.assert_allclose(native_nib.affine, raw.affine, atol=1e-4)
            coreg = sitk.ReadTransform(str(track(func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt")))
            composite = sitk.CompositeTransform([coreg, a2m])
            standard = sitk.ReadImage(str(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz"))
            reconstructed = array(sitk.Resample(native, standard, composite, sitk.sitkLinear, 0, sitk.sitkFloat32))
            saved = array(standard)
            common = ((reconstructed > np.percentile(reconstructed[reconstructed > 0], 15)) &
                      (saved > np.percentile(saved[saved > 0], 15)))
            agreement = correlation(reconstructed, saved, common)
            if agreement < .93:
                raise ValueError(f"BOLD transform reproduction failed: {subject} {condition} r={agreement}")
            transformations.append({"subject": subject, "session": condition,
                                    "bold_reference_correlation": agreement})
            for side in ["l", "r"]:
                source_path = next(r["path"] for r in block if r["session_number"] == session and r["side"] == side)
                source_image = nib.load(track(source_path))
                values = np.asarray(source_image.dataobj) > 0
                assert values.shape == native_nib.shape
                corners = np.array(np.meshgrid(*[[-(d-1)/2, (d-1)/2] for d in values.shape])).reshape(3, -1).T
                grid_error = float(np.linalg.norm(corners @ (source_image.affine[:3, :3] - native_nib.affine[:3, :3]).T, axis=1).max())
                if grid_error >= .5:
                    raise ValueError(f"Native mask lattice mismatch: {subject} {condition} {side}")
                restored = mask_image(values, native)
                mapped = sitk.Resample(restored, reference, composite, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
                data = array(mapped).astype(float)
                native_t1_target = array(sitk.Resample(restored, t1, coreg, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8))
                if not data.sum():
                    raise ValueError(f"Empty mapped target: {subject} {condition} {side}")
                center = nib.affines.apply_affine(atlas.affine, np.argwhere(data).mean(0))
                source_point = nib.affines.apply_affine(native_nib.affine, np.argwhere(values).mean(0)) * [-1, -1, 1]
                t1point = coreg.GetInverse().TransformPoint(source_point.tolist())
                mnipoint_lps = np.asarray(m2a.TransformPoint(t1point))
                roundtrip = np.linalg.norm(np.asarray(a2m.TransformPoint(mnipoint_lps.tolist())) - t1point)
                centroid_difference = np.linalg.norm(center - mnipoint_lps * [-1, -1, 1])
                if roundtrip >= .5 or centroid_difference >= 4:
                    raise ValueError(f"Transform centroid consistency failed for {subject} {condition} {side}")
                assert center[0] < 0 if side == "l" else center[0] > 0
                masks[side].append(data)
                sitk.WriteImage(mapped, str(destination / f"sub-{subject}_{condition}_{side}_mapped.nii.gz"))
                qc.append({"subject": subject, "session": condition, "side": side,
                           "native_voxels": int(values.sum()), "mni_voxels": int(data.sum()),
                           "mni_x": center[0], "mni_y": center[1], "mni_z": center[2],
                           "native_t1_brainmask_fraction": fraction(native_t1_target, t1mask),
                           "cortical_fraction": fraction(data, cortical),
                           "linear_grid_corner_error_mm": grid_error,
                           "inverse_transform_error_mm": float(roundtrip),
                           "centroid_agreement_mm": float(centroid_difference)})
        for side, images in masks.items():
            mean = np.mean(images, axis=0)
            cortical_fraction = fraction(mean, cortical)
            eligible = cortical_fraction > .1  # Preserve the original anatomy-only rule.
            coverage.append({"subject_id": f"aachen_itbs_{subject}", "side": side,
                             "cortical_fraction": cortical_fraction, "eligible": eligible})
            nib.save(nib.Nifti1Image(mean.astype("float32"), atlas.affine),
                     destination / f"sub-{subject}_{side}_fixed_mask.nii.gz")
            if eligible:
                counts = np.bincount(labels.ravel(), weights=mean.ravel(), minlength=455)[1:455]
                counts[400:] = 0
                counts /= counts.sum()
                for roi, weight in enumerate(counts, 1):
                    weights.append({"subject_id": f"aachen_itbs_{subject}", "side": side,
                                    "roi": roi, "weight": weight, "cortical_fraction": cortical_fraction,
                                    "source_sessions": ",".join(map(str, chosen))})
        print("Mapped", subject, "from corrected derivatives", flush=True)
    for name, data in [("mask_mapping_qc", qc), ("transform_validation", transformations),
                       ("subject_region_weights", weights), ("parcel_coverage", coverage)]:
        pd.DataFrame(data).to_csv(destination / f"{name}.csv", index=False)
    provenance = {"created_utc": datetime.now(timezone.utc).isoformat(),
                  "derivatives_root": str(args.derivatives_root), "independent_audit": str(args.audit_dir),
                  "subjects": subjects, "all_deposited_usable_subjects_mapped": len(subjects) == 8,
                  "input_sha256": hashes, "audit_provenance_sha256": hashlib.sha256((args.audit_dir / "provenance.json").read_bytes()).hexdigest(),
                  "functional_outcomes_read": False, "production_files_changed": False,
                  "selection": "Available baseline masks averaged within subject; earliest available source only if neither baseline was deposited.",
                  "eligibility": "Original >10% cortical atlas overlap rule retained. Low coverage is recorded and never relocated or replaced.",
                  "geometry": "Deposited mean-functional voxel lattice restored to matching raw/fMRIPrep native lattice, then mapped by corrected fMRIPrep transforms.",
                  "limitation": "Mean-functional versus reference residual motion and ROI partial volume remain; these are seed masks, not electric-field models.",
                  "required_for_effect_analysis": "These weights must be used with timeseries and outcomes reprocessed through the same corrected anatomical chain."}
    (destination / "anatomy_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print("Saved", destination, flush=True)


if __name__ == "__main__":
    main()
