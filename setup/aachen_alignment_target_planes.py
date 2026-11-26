"""Bilateral three-plane target evidence for completed Aachen anatomy audits.

This preserves deposited voxel masks and uses the same corrected transforms.
It adds anatomical views, not alternative target locations or functional tests.
"""
from argparse import ArgumentParser
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
import aachen_alignment_anatomical_review as audit

ROOT = audit.ROOT
OUT = audit.REVIEW / "corrected_cohort"
DERIV = ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel(ax, data, target, axis, index, title, cortex=None):
    ax.imshow(np.take(data, index, axis=axis).T, origin="lower", cmap="gray",
              vmin=0, vmax=np.percentile(data[data > 0], 98))
    if cortex is not None:
        plane = np.take(cortex, index, axis=axis).T
        if plane.min() < .5 < plane.max():
            ax.contour(plane, levels=[.5], colors=["#39c5dc"], linewidths=.5)
    plane = np.take(target, index, axis=axis).T
    if plane.max() > 0:
        ax.contour(plane, levels=[.5], colors=["#ff6339"], linewidths=1.1)
    ax.set_title(title, fontsize=9)
    ax.axis("off")


def main():
    parser = ArgumentParser()
    parser.add_argument("--subjects", nargs="*")
    args = parser.parse_args()
    selected = audit.selected_masks()
    subjects = sorted({r["subject"] for r in selected})
    if args.subjects:
        subjects = [s for s in subjects if s in {x.zfill(2) for x in args.subjects}]
    queue = json.loads((ROOT / "logs/aachen_corrected_20261125/queue_status.json").read_text())
    atlas_path = ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz"
    atlas = nib.load(atlas_path)
    labels = np.asarray(atlas.dataobj)
    cortex = (labels >= 1) & (labels <= 400)
    for subject in subjects:
        out = OUT / f"sub-{subject}"
        done = queue.get("finished", {}).get(subject, {})
        if not ((out / "metrics.json").exists() and done.get("success") is True
                and done.get("state", {}).get("ExitCode") == 0
                and done.get("inventory", {}).get("complete") is True):
            continue
        receipt_path = out / "bilateral_target_planes_provenance.json"
        if receipt_path.exists():
            continue
        inputs, outputs = {}, {}

        def track(path):
            inputs[str(path)] = digest(path)
            return path

        anat = DERIV / f"sub-{subject}/ses-anat/anat"
        prefix = f"sub-{subject}_ses-anat"
        native = audit.sitk.ReadImage(str(track(anat / f"{prefix}_desc-preproc_T1w.nii.gz")))
        native_data = audit.array(native)
        standard_image = nib.load(track(anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz"))
        standard_data = np.asarray(standard_image.dataobj)
        np.testing.assert_allclose(standard_image.affine, atlas.affine, atol=1e-4)
        track(atlas_path)
        block = [r for r in selected if r["subject"] == subject]
        for session_number in sorted({r["session_number"] for r in block}):
            session = audit.SESSIONS[session_number]
            func = DERIV / f"sub-{subject}/ses-{session}/func"
            stem = f"sub-{subject}_ses-{session}_task-rest"
            boldref = audit.sitk.ReadImage(str(track(func / f"{stem}_desc-coreg_boldref.nii.gz")))
            coreg = audit.sitk.ReadTransform(str(track(func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt")))
            fig, axes = plt.subplots(4, 3, figsize=(14, 15), layout="constrained")
            for side_number, side in enumerate(["l", "r"]):
                source = next(r["path"] for r in block if r["session_number"] == session_number and r["side"] == side)
                values = np.asarray(nib.load(track(source)).dataobj) > 0
                restored = audit.mask_image(values, boldref)
                native_target = audit.array(audit.sitk.Resample(restored, native, coreg,
                    audit.sitk.sitkNearestNeighbor, 0, audit.sitk.sitkUInt8)) > 0
                standard_target = np.asarray(nib.load(track(out / f"sub-{subject}_{session}_{side}_diagnostic_target_mni.nii.gz")).dataobj) > 0
                for space_number, (data, target) in enumerate([(native_data, native_target), (standard_data, standard_target)]):
                    center = np.rint(np.argwhere(target).mean(0)).astype(int)
                    for axis, plane in enumerate(["sagittal", "coronal", "axial"]):
                        title = f"{'Left' if side == 'l' else 'Right'} | {'native T1' if space_number == 0 else 'MNI T1 / cortical atlas'} | {plane}"
                        panel(axes[side_number*2+space_number, axis], data, target, axis, center[axis], title,
                              cortex if space_number == 1 else None)
            fig.suptitle(f"sub-{subject} {session}: deposited bilateral analysis ROIs, unchanged source voxels\nOrange: ROI; cyan: cortical atlas. Geometry evidence does not independently identify coil position.")
            path = out / f"sub-{subject}_{session}_bilateral_target_planes.png"
            fig.savefig(path, dpi=120)
            plt.close(fig)
            outputs[str(path)] = digest(path)
        receipt_path.write_text(json.dumps({"subject": subject, "created_utc": datetime.now(timezone.utc).isoformat(),
            "script_sha256": digest(Path(__file__)), "derivative_root": str(DERIV),
            "input_sha256": inputs, "output_sha256": outputs,
            "scope": "Three-plane anatomical localization of unchanged deposited left/right analysis ROI masks."}, indent=2), encoding="utf-8")
        print(json.dumps({"bilateral_target_planes_ready": subject, "figures": len(outputs)}), flush=True)


if __name__ == "__main__":
    main()
