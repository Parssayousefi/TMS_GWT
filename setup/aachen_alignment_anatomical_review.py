"""Read-only Aachen anatomical audit; write only a new review directory.

This examines anatomical support and saved transforms, not functional outcomes.
It neither repairs production derivatives nor computes target effect statistics.
Example pilot comparison: --derivatives-root <clean-pilot-fmriprep> --subjects 06
                         --label pilot_skip
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
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tmp/target_roi_deps"))
import SimpleITK as sitk
sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(4)

REVIEW = ROOT / "output/aachen_alignment_review"
SESSIONS = {1: "pre90", 2: "post90", 3: "pre120", 4: "post120"}


def array(image):
    return sitk.GetArrayFromImage(image).transpose(2, 1, 0)


def mask_image(data, reference):
    image = sitk.GetImageFromArray(np.asarray(data, dtype="uint8").transpose(2, 1, 0))
    image.CopyInformation(reference)
    return image


def dice(a, b):
    return float(2 * np.count_nonzero(a & b) / (np.count_nonzero(a) + np.count_nonzero(b)))


def fraction(mask, support):
    return float(np.sum(mask[support]) / np.sum(mask))


def correlation(a, b, support):
    return float(np.corrcoef(a[support], b[support])[0, 1])


def geometry(image):
    return {
        "shape": list(image.shape), "zooms": [float(x) for x in image.header.get_zooms()],
        "affine": image.affine.tolist(), "qform": image.get_qform().tolist(),
        "sform": image.get_sform().tolist(),
        "qform_code": int(image.header["qform_code"]),
        "sform_code": int(image.header["sform_code"]),
    }


def display(ax, data, mask=None, title="", index=None, axis=0):
    if index is None:
        index = data.shape[axis] // 2
    positive = data[data > 0]
    ax.imshow(np.take(data, index, axis=axis).T, origin="lower", cmap="gray",
              vmin=0, vmax=np.percentile(positive, 98) if len(positive) else 1)
    if mask is not None:
        plane = np.take(mask, index, axis=axis).T
        if plane.max() > 0:
            ax.contour(plane, levels=[0.5], colors=["#ff6339"], linewidths=1)
    ax.set_title(title, fontsize=9)
    ax.axis("off")


def selected_masks():
    files = []
    for path in sorted((ROOT / "data/aachen_target_rois").glob("*rFC*.nii")):
        match = re.search(r"sub-(\d+)_ses-(\d+)_rFC([lr])", path.name)
        files.append({"subject": match[1], "session_number": int(match[2]),
                      "side": match[3], "path": path})
    all_files = pd.DataFrame(files)
    result = []
    for subject, group in all_files.groupby("subject"):
        if subject == "02":
            continue  # Existing functional exclusion, unrelated to anatomy/target findings.
        sessions = sorted(set(group.session_number) & {1, 3}) or [int(group.session_number.min())]
        result.extend(group[group.session_number.isin(sessions)].to_dict("records"))
    return result


def main():
    parser = ArgumentParser()
    parser.add_argument("--derivatives-root", type=Path,
                        default=ROOT / "data/aachen_itbs/derivatives/fmriprep")
    parser.add_argument("--subjects", nargs="*")
    parser.add_argument("--label", default="original")
    parser.add_argument("--save-targets", action="store_true",
                        help="Save diagnostic target NIfTIs inside the review output only")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.label):
        raise ValueError("Label must be a plain directory name")
    out = REVIEW / args.label
    out.mkdir(parents=True, exist_ok=True)
    subject_ids = args.subjects or [f"{i:02d}" for i in range(1, 17)]
    source_masks = selected_masks()
    anatomy_rows, session_rows, target_rows, configs, headers, hashes = [], [], [], [], [], {}
    montage, montage_axes = plt.subplots(4, 4, figsize=(14, 13), layout="constrained")
    for ax in montage_axes.ravel():
        ax.axis("off")

    def track(path):
        path = Path(path)
        try:
            key = path.relative_to(ROOT).as_posix()
        except ValueError:
            key = str(path)
        hashes[key] = hashlib.sha256(path.read_bytes()).hexdigest()
        return path

    atlas = nib.load(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    atlas_labels = np.asarray(atlas.dataobj)
    atlas_cortex = (atlas_labels >= 1) & (atlas_labels <= 400)
    for num, subject in enumerate(subject_ids):
        base = args.derivatives_root / f"sub-{subject}"
        anat = base / "ses-anat/anat"
        prefix = f"sub-{subject}_ses-anat"
        native_path = anat / f"{prefix}_desc-preproc_T1w.nii.gz"
        if not native_path.exists():
            print("No native anatomy:", subject, flush=True)
            continue
        raw_path = ROOT / f"data/aachen_itbs/sub-{subject}/ses-anat/anat/{prefix}_T1w.nii.gz"
        raw_original = nib.load(track(raw_path))
        raw = nib.as_closest_canonical(raw_original)
        native_nib = nib.load(native_path)
        t1 = np.asarray(native_nib.dataobj, dtype=np.float32)
        raw_data = np.asarray(raw.dataobj)
        native = sitk.ReadImage(str(native_path))
        mask_path = anat / f"{prefix}_desc-brain_mask.nii.gz"
        mask = np.asarray(nib.load(track(mask_path)).dataobj) > 0
        original_support = t1 > 0
        # This is anatomical image support, not a substitute segmentation/GM mask.
        filled = ndimage.binary_fill_holes(original_support)
        distance = ndimage.distance_transform_edt(filled, sampling=native_nib.header.get_zooms())
        surface_shell = original_support & (distance <= 8)
        intensity_threshold = np.percentile(t1[original_support], 20)
        above_p20 = t1 > intensity_threshold
        voxel_mm3 = abs(np.linalg.det(native_nib.affine[:3, :3]))
        config_files = sorted((base / "log").glob("*/fmriprep.toml"))
        for config_path in config_files:
            text = config_path.read_text(encoding="utf-8")
            match = re.search(r'^skull_strip_t1w\s*=\s*"([^\"]+)"', text, re.M)
            configs.append({"subject": subject, "config_path": str(config_path),
                            "skull_strip_t1w": match[1] if match else "missing"})
        headers.append({"subject": subject, "raw": geometry(raw_original),
                        "raw_canonical": geometry(raw), "native_preproc": geometry(native_nib)})
        row = {
            "subject": subject, "included_functional_cohort": subject != "02",
            "raw_nonzero_mL": float(np.count_nonzero(raw_data) * abs(np.linalg.det(raw.affine[:3, :3])) / 1000),
            "native_nonzero_mL": float(original_support.sum() * voxel_mm3 / 1000),
            "fmriprep_native_mask_mL": float(mask.sum() * voxel_mm3 / 1000),
            "fraction_native_support_retained": fraction(original_support, mask),
            "fraction_above_native_p20_retained": fraction(above_p20, mask),
            "fraction_outer_8mm_support_retained": fraction(surface_shell, mask),
            "raw_native_array_support_dice": dice(raw_data > 0, original_support) if raw_data.shape == t1.shape else np.nan,
            "raw_native_affine_max_abs_diff": float(np.abs(raw.affine - native_nib.affine).max()),
        }
        center = np.rint(np.argwhere(original_support).mean(0)).astype(int)
        if num < 16:
            display(montage_axes.ravel()[num], t1, mask,
                    f"sub-{subject}: support retained {row['fraction_native_support_retained']:.1%}",
                    index=int(center[0]))
        mni_path = anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz"
        mni_mask_path = anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz"
        transform_path = anat / f"{prefix}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"
        if not (mni_path.exists() and transform_path.exists()):
            anatomy_rows.append(row)
            continue
        standard = sitk.ReadImage(str(mni_path))
        standard_nib = nib.load(mni_path)
        std_t1 = array(standard)
        std_mask = np.asarray(nib.load(mni_mask_path).dataobj) > 0
        a2m = sitk.ReadTransform(str(track(transform_path)))
        full_mni = array(sitk.Resample(native, standard, a2m, sitk.sitkLinear, 0, sitk.sitkFloat32))
        reconstructed_mask = array(sitk.Resample(mask_image(mask, native), standard, a2m,
                                                sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)) > 0
        full_mni_support = array(sitk.Resample(mask_image(original_support, native), standard, a2m,
                                              sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)) > 0
        row.update({"native_to_saved_mni_mask_dice": dice(reconstructed_mask, std_mask),
                    "native_to_saved_mni_t1_r_within_saved_mask": correlation(full_mni, std_t1, std_mask),
                    "mni_t1_atlas_affine_max_abs_diff": float(np.abs(standard_nib.affine - atlas.affine).max())})
        anatomy_rows.append(row)
        subject_masks = [x for x in source_masks if x["subject"] == subject]
        sessions = sorted({SESSIONS[x["session_number"]] for x in subject_masks})
        # Also inspect every available saved BOLD brainmask, including subjects with no deposited target.
        sessions = sorted(set(sessions) | {p.parent.parent.name.replace("ses-", "") for p in base.glob("ses-*/func/*space-MNI*desc-brain_mask.nii.gz")})
        for condition in sessions:
            func = base / f"ses-{condition}/func"
            stem = f"sub-{subject}_ses-{condition}_task-rest"
            bold_mni_mask_path = func / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz"
            bold_mni_path = func / f"{stem}_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz"
            coreg_path = func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"
            bold_native_path = func / f"{stem}_desc-coreg_boldref.nii.gz"
            if not all(p.exists() for p in [bold_mni_mask_path, bold_mni_path, coreg_path, bold_native_path]):
                continue
            bold_mni_nib = nib.load(bold_mni_path)
            bold_mni = np.asarray(bold_mni_nib.dataobj)
            bold_mni_mask = np.asarray(nib.load(bold_mni_mask_path).dataobj) > 0
            assert bold_mni.shape == std_t1.shape
            coreg = sitk.ReadTransform(str(track(coreg_path)))
            bold_native = sitk.ReadImage(str(bold_native_path))
            bold_native_nib = nib.load(bold_native_path)
            composite = sitk.CompositeTransform([coreg, a2m])
            reconstructed_bold = array(sitk.Resample(bold_native, standard, composite,
                                                    sitk.sitkLinear, 0, sitk.sitkFloat32))
            common = ((reconstructed_bold > np.percentile(reconstructed_bold[reconstructed_bold > 0], 15)) &
                      (bold_mni > np.percentile(bold_mni[bold_mni > 0], 15)))
            session_rows.append({
                "subject": subject, "session": condition,
                "bold_to_masked_t1_mni_dice": dice(bold_mni_mask, std_mask),
                "bold_to_full_t1_support_mni_dice": dice(bold_mni_mask, full_mni_support),
                "bold_t1_mask_volume_ratio": float(bold_mni_mask.sum() / std_mask.sum()),
                "bold_transform_reproduction_r": correlation(reconstructed_bold, bold_mni, common),
                "bold_t1_mni_affine_max_abs_diff": float(np.abs(bold_mni_nib.affine - standard_nib.affine).max()),
            })
            for source in [x for x in subject_masks if SESSIONS[x["session_number"]] == condition]:
                source_nib = nib.load(track(source["path"]))
                source_data = np.asarray(source_nib.dataobj) > 0
                assert source_data.shape == bold_native_nib.shape
                corners = np.array(np.meshgrid(*[[-(d-1)/2, (d-1)/2] for d in source_data.shape])).reshape(3, -1).T
                grid_error = float(np.linalg.norm(corners @ (source_nib.affine[:3, :3] - bold_native_nib.affine[:3, :3]).T, axis=1).max())
                assert grid_error < .5
                native_target = mask_image(source_data, bold_native)
                t1_target = array(sitk.Resample(native_target, native, coreg, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8))
                mni_target = array(sitk.Resample(native_target, standard, composite, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8))
                mni_centroid = nib.affines.apply_affine(standard_nib.affine, np.argwhere(mni_target).mean(0))
                same_atlas_grid = (mni_target.shape == atlas_cortex.shape and
                                   np.allclose(standard_nib.affine, atlas.affine, atol=1e-5))
                if args.save_targets:
                    nib.save(nib.Nifti1Image(mni_target.astype("uint8"), standard_nib.affine),
                             out / f"sub-{subject}_{condition}_{source['side']}_diagnostic_target_mni.nii.gz")
                target_rows.append({
                    "subject": subject, "session": condition, "side": source["side"],
                    "linear_grid_corner_error_mm": grid_error,
                    "fraction_target_in_native_full_t1_support": fraction(t1_target, original_support),
                    "fraction_target_in_native_fmriprep_t1_mask": fraction(t1_target, mask),
                    "fraction_target_in_mni_full_t1_support": fraction(mni_target, full_mni_support),
                    "fraction_target_in_mni_fmriprep_t1_mask": fraction(mni_target, std_mask),
                    "fraction_target_in_mni_bold_mask": fraction(mni_target, bold_mni_mask),
                    "fraction_target_in_existing_cortical_atlas": fraction(mni_target, atlas_cortex) if same_atlas_grid else np.nan,
                    "mni_centroid_x": float(mni_centroid[0]),
                    "mni_centroid_y": float(mni_centroid[1]),
                    "mni_centroid_z": float(mni_centroid[2]),
                })
                if source["side"] == "l":
                    target_center = np.rint(np.argwhere(t1_target).mean(0)).astype(int)
                    mni_center = np.rint(np.argwhere(mni_target).mean(0)).astype(int)
                    fig, axes = plt.subplots(1, 4, figsize=(15, 4), layout="constrained")
                    display(axes[0], t1, t1_target, "Full native T1 + target", index=target_center[0])
                    display(axes[1], full_mni, mni_target, "Full T1 via current MNI transform", index=mni_center[0])
                    display(axes[2], std_t1, mni_target, "Saved skull-stripped MNI T1", index=mni_center[0])
                    display(axes[3], bold_mni, mni_target, "Saved MNI BOLD", index=mni_center[0])
                    fig.suptitle(f"sub-{subject} {condition}: deposited left target; no target relocation")
                    fig.savefig(out / f"sub-{subject}_{condition}_target_chain.png", dpi=140)
                    plt.close(fig)
        print(f"Audited sub-{subject}: native support retained {row['fraction_native_support_retained']:.3f}", flush=True)
    montage.suptitle("Aachen T1 anatomy: existing fMRIPrep brainmask in orange\nFull native T1 intensity support shown; masks visibly omit cortex/cerebellum")
    montage.savefig(out / "all_subjects_native_t1_brainmask.png", dpi=150)
    plt.close(montage)
    for name, data in [("anatomical_support", anatomy_rows), ("functional_anatomical_agreement", session_rows),
                       ("target_support", target_rows), ("effective_skullstrip_settings", configs)]:
        pd.DataFrame(data).to_csv(out / f"{name}.csv", index=False)
    (out / "headers.json").write_text(json.dumps(headers, indent=2), encoding="utf-8")
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "derivatives_root": str(args.derivatives_root),
        "software": {"SimpleITK": sitk.Version_VersionString(), "nibabel": nib.__version__, "numpy": np.__version__},
        "functional_outcomes_read": False, "production_inputs_modified": False,
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": hashes,
        "support_definition": "Nonzero native T1 image support; not an independent gray-matter segmentation. Outer shell is original support within 8 mm of the hole-filled external boundary.",
        "target_selection": "Available baseline masks averaged nowhere in this audit; earliest available session only when neither baseline was deposited. Existing sub-02 functional exclusion retained.",
        "purpose": "Diagnose anatomy/mapping, not validate a repaired pipeline or release new target weights.",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print("Saved", out, flush=True)


if __name__ == "__main__":
    main()
