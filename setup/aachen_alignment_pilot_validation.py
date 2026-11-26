"""Compare an audited corrected Aachen pilot with the original spatial chain.

Outputs tables and anatomical figures only, under aachen_alignment_review.
Run after the corrected audit and mapper have completed.
"""
from argparse import ArgumentParser
from pathlib import Path
import json
import re
import sys

import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
from aachen_alignment_anatomical_review import ROOT, REVIEW, array, mask_image, sitk


def show(ax, data, axis, index, title, atlas=None, target=None):
    ax.imshow(np.take(data, index, axis=axis).T, origin="lower", cmap="gray", vmin=0,
              vmax=np.percentile(data[data > 0], 98))
    for overlay, color, width in [(atlas, "#39c5dc", .8), (target, "#ff6339", 1.3)]:
        if overlay is not None:
            plane = np.take(overlay, index, axis=axis).T
            if plane.min() < .5 < plane.max():
                ax.contour(plane, levels=[.5], colors=[color], linewidths=width)
    ax.axis("off")
    ax.set_title(title, fontsize=10)


def main():
    parser = ArgumentParser()
    parser.add_argument("--corrected-root", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--mapping-dir", type=Path, required=True)
    parser.add_argument("--subject", default="06")
    parser.add_argument("--label", default="pilot_validation")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.label):
        raise ValueError("Unsafe label")
    subject = args.subject.zfill(2)
    out = REVIEW / args.label
    out.mkdir(parents=True, exist_ok=True)
    original_root = ROOT / "data/aachen_itbs/derivatives/fmriprep"
    original_audit = REVIEW / f"original_sub{subject}"
    if not original_audit.exists():
        raise ValueError("Run representative original audit with current audit script first")
    summary = {"subject": subject, "functional_outcomes_read": False,
               "production_inputs_modified": False, "corrected_root": str(args.corrected_root)}
    for table in ["anatomical_support", "functional_anatomical_agreement", "target_support"]:
        old = pd.read_csv(original_audit / f"{table}.csv", dtype={"subject": str})
        new = pd.read_csv(args.audit_dir / f"{table}.csv", dtype={"subject": str})
        old = old[old.subject.eq(subject)]
        new = new[new.subject.eq(subject)]
        keys = ["subject"] + (["session"] if "session" in old else []) + (["side"] if "side" in old else [])
        merged = old.merge(new, on=keys, suffixes=("_original", "_corrected"), validate="one_to_one")
        merged.to_csv(out / f"{table}_comparison.csv", index=False)
        summary[table + "_rows"] = len(merged)

    atlas_path = ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz"
    atlas_nib = nib.load(atlas_path)
    atlas_sitk = sitk.ReadImage(str(atlas_path))
    cortex = mask_image((np.asarray(atlas_nib.dataobj) >= 1) & (np.asarray(atlas_nib.dataobj) <= 400), atlas_sitk)
    template_dir = ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym"
    template_path = template_dir / "tpl-MNI152NLin2009cAsym_res-02_desc-brain_T1w.nii.gz"
    if not template_path.exists() or template_path.stat().st_size == 0:
        # TemplateFlow may leave an unmaterialized zero-byte placeholder.
        # The cached full T1 and its official mask define the same reference.
        template_path = template_dir / "tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz"
    template = nib.load(template_path)
    template_mask = np.asarray(nib.load(template_dir / "tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz").dataobj) > 0
    np.testing.assert_allclose(template.affine, atlas_nib.affine)
    template_data = np.asarray(template.dataobj) * template_mask
    summary["template_T1_path"] = str(template_path)
    prefix = f"sub-{subject}_ses-anat"
    corrected_anat = args.corrected_root / f"sub-{subject}/ses-anat/anat"
    native_path = corrected_anat / f"{prefix}_desc-preproc_T1w.nii.gz"
    native_nib = nib.load(native_path)
    native = sitk.ReadImage(str(native_path))
    native_data = np.asarray(native_nib.dataobj)
    gm = np.asarray(nib.load(corrected_anat / f"{prefix}_label-GM_probseg.nii.gz").dataobj)
    stem = f"sub-{subject}_ses-pre90_task-rest"
    func = args.corrected_root / f"sub-{subject}/ses-pre90/func"
    bold_native = sitk.ReadImage(str(func / f"{stem}_desc-coreg_boldref.nii.gz"))
    coreg = sitk.ReadTransform(str(func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"))
    source_path = next((ROOT / "data/aachen_target_rois").glob(f"sub-{subject}_ses-1_rFCl*.nii"))
    source = np.asarray(nib.load(source_path).dataobj) > 0
    native_target = array(sitk.Resample(mask_image(source, bold_native), native, coreg,
                                       sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)) > 0
    native_center = np.rint(np.argwhere(native_target).mean(0)).astype(int)
    corrected_target = np.asarray(nib.load(args.mapping_dir / f"sub-{subject}_pre90_l_mapped.nii.gz").dataobj) > 0
    mni_center = np.rint(np.argwhere(corrected_target).mean(0)).astype(int)
    mni_world = nib.affines.apply_affine(atlas_nib.affine, mni_center)
    native_fig, native_axes = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    standard_fig, standard_axes = plt.subplots(3, 3, figsize=(12, 10), layout="constrained")
    for row, (name, root) in enumerate([("Original force", original_root), ("Corrected skip", args.corrected_root)]):
        anat = root / f"sub-{subject}/ses-anat/anat"
        old_native = nib.load(anat / f"{prefix}_desc-preproc_T1w.nii.gz")
        np.testing.assert_allclose(old_native.affine, native_nib.affine, atol=1e-4)
        m2a = sitk.ReadTransform(str(anat / f"{prefix}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"))
        native_atlas = array(sitk.Resample(cortex, native, m2a, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)) > 0
        summary[name + "_native_atlas_mean_corrected_GM_probability"] = float(gm[native_atlas].mean())
        summary[name + "_native_atlas_voxels"] = int(native_atlas.sum())
        standardized = nib.load(anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz")
        np.testing.assert_allclose(standardized.affine, template.affine, atol=1e-4)
        standard_data = np.asarray(standardized.dataobj)
        brainmask = np.asarray(nib.load(anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz").dataobj) > 0
        summary[name + "_template_mask_dice"] = float(2*(brainmask & template_mask).sum()/(brainmask.sum()+template_mask.sum()))
        summary[name + "_T1_template_intensity_r"] = float(np.corrcoef(standard_data[template_mask], template_data[template_mask])[0, 1])
        for axis in range(3):
            show(native_axes[row, axis], native_data, axis, native_center[axis],
                 f"{name}: atlas on same native T1", native_atlas, native_target)
            show(standard_axes[row, axis], standard_data, axis, mni_center[axis],
                 f"{name}; {'xyz'[axis]}={mni_world[axis]:.1f} mm", template_mask,
                 corrected_target if row == 1 else None)
    for axis in range(3):
        show(standard_axes[2, axis], template_data, axis, mni_center[axis],
             f"MNI template; {'xyz'[axis]}={mni_world[axis]:.1f} mm", None, corrected_target)
    native_fig.suptitle(f"sub-{subject}: cortical atlas in cyan; unchanged deposited left target in orange\nBoth rows use the same corrected native T1 image and same native target")
    native_fig.savefig(out / "native_atlas_placement_comparison.png", dpi=150)
    plt.close(native_fig)
    standard_fig.suptitle(f"sub-{subject}: same MNI slices through corrected left target\nTemplate brain boundary in cyan; corrected target in orange")
    standard_fig.savefig(out / "standard_t1_template_comparison.png", dpi=150)
    plt.close(standard_fig)
    (out / "validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
