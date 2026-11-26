"""Create anatomical/coverage review material after a complete ds001927 pilot.

Does not automatically approve images or run treatment statistics.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DERIV = ROOT / "data/ds001927/derivatives/fmriprep_23_2_1"
SPACE = "MNI152NLin2009cAsym_res-2"


def show(ax, data, axis, index, title, contours=()):
    plane = np.take(data, index, axis=axis).T
    positive = data[data > 0]
    ax.imshow(plane, origin="lower", cmap="gray", vmin=0,
              vmax=np.percentile(positive, 98) if positive.size else 1)
    for mask, color in contours:
        sl = np.take(mask, index, axis=axis).T
        if sl.min() < .5 < sl.max():
            ax.contour(sl, levels=[.5], colors=[color], linewidths=.6)
    ax.set_title(title, fontsize=8)
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="s001", choices=[f"s{i:03d}" for i in range(1, 24)])
    args = ap.parse_args()
    sub = "sub-" + args.subject
    out = ROOT / "output/ds001927_pilot_review_20261126" / sub
    out.mkdir(parents=True, exist_ok=True)
    template_dir = ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym"
    atlas_path = ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz"
    atlas = nib.load(atlas_path)
    labels = np.asarray(atlas.dataobj).astype(int)
    counts = np.bincount(labels.ravel(), minlength=455)[1:455]
    assert len(counts) == 454 and np.all(counts > 0)
    anat = DERIV / sub / "anat"
    source_paths = []

    def read(path, standard=True):
        im = nib.load(path)
        source_paths.append(path)
        if standard:
            np.testing.assert_allclose(im.affine, atlas.affine, atol=1e-4)
            assert im.shape == atlas.shape, (path, im.shape, atlas.shape)
        a = im.get_fdata(dtype=np.float32)
        assert np.isfinite(a).all(), path
        return a

    template = read(template_dir / "tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz")
    template_mask = read(template_dir / "tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz") > 0
    t1 = read(anat / f"{sub}_space-{SPACE}_desc-preproc_T1w.nii.gz")
    t1mask = read(anat / f"{sub}_space-{SPACE}_desc-brain_mask.nii.gz") > 0
    gm = read(anat / f"{sub}_space-{SPACE}_label-GM_probseg.nii.gz")
    native = read(anat / f"{sub}_desc-preproc_T1w.nii.gz", False)
    native_mask = read(anat / f"{sub}_desc-brain_mask.nii.gz", False) > 0
    native_gm = read(anat / f"{sub}_label-GM_probseg.nii.gz", False)
    summary = dict(subject=sub, created_utc=datetime.now(timezone.utc).isoformat(),
                   template_mask_dice=float(2 * (t1mask & template_mask).sum() / (t1mask.sum() + template_mask.sum())),
                   t1_template_intensity_r=float(np.corrcoef(t1[template_mask], template[template_mask])[0, 1]),
                   cortical_atlas_mean_gm_probability=float(gm[(labels > 0) & (labels <= 400)].mean()),
                   analysis_release=False,
                   status="review material generated; requires visual inspection and independent motion gate",
                   limitation="Mask/reference support does not measure temporal reliability or prove registration accuracy")
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    center = np.rint(np.argwhere(native_mask).mean(axis=0)).astype(int)
    for axis in range(3):
        for row in range(2):
            contour = [(native_mask if row == 0 else native_gm > .5, "#ffac54")]
            show(axes[row, axis], native, axis, center[axis],
                 f"{sub}: {'brain boundary' if row == 0 else 'GM >0.5'}", contour)
    fig.suptitle("Extraction and tissue segmentation on common native T1")
    fig.tight_layout()
    fig.savefig(out / "native_extraction_and_segmentation.png", dpi=140)
    plt.close(fig)
    zvals = [-20, 0, 20, 40, 60]
    indices = [int(round(nib.affines.apply_affine(np.linalg.inv(atlas.affine), [0, 0, z])[2])) for z in zvals]
    fig, axes = plt.subplots(3, 5, figsize=(17, 10))
    for row, (a, title) in enumerate([(template * template_mask, "MNI template"),
                                    (t1, "Subject T1"), (t1, "GM segmentation")]):
        for col, (z, index) in enumerate(zip(zvals, indices)):
            overlays = [(template_mask if row < 2 else gm > .5, "#22bdce")]
            show(axes[row, col], a, 2, index, f"{title}, z={z} mm", overlays)
    fig.tight_layout()
    fig.savefig(out / "standard_alignment.png", dpi=140)
    plt.close(fig)
    fig, axes = plt.subplots(6, 5, figsize=(17, 17))
    rows = []
    for row, (ses, task) in enumerate((s, t) for s in ["ctr", "fro", "occ"] for t in ["restpre", "restpost"]):
        stem = f"{sub}_ses-{ses}_task-{task}"
        func = DERIV / sub / f"ses-{ses}" / "func"
        mask = read(func / f"{stem}_space-{SPACE}_desc-brain_mask.nii.gz") > 0
        ref = read(func / f"{stem}_space-{SPACE}_boldref.nii.gz")
        inmask = np.bincount(labels[mask], minlength=455)[1:455]
        nonzero = np.bincount(labels[ref > 0], minlength=455)[1:455]
        for roi in range(454):
            rows.append(dict(subject=sub, session=ses, task=task, roi=roi + 1,
                             atlas_voxels=int(counts[roi]), mask_coverage=float(inmask[roi] / counts[roi]),
                             nonzero_reference_coverage=float(nonzero[roi] / counts[roi])))
        for col, (z, index) in enumerate(zip(zvals, indices)):
            show(axes[row, col], ref, 2, index, f"{ses} {task}, z={z} mm",
                 [(template_mask, "#22bdce"), (mask, "#ffac54")])
    fig.suptitle("BOLD references: MNI template boundary cyan; functional mask orange")
    fig.tight_layout(rect=(0, 0, 1, .98))
    fig.savefig(out / "six_run_bold_alignment_and_coverage.png", dpi=140)
    plt.close(fig)
    df = pd.DataFrame(rows)
    df.to_csv(out / "atlas_coverage.csv", index=False)
    summary["zero_mask_parcel_runs"] = int((df.mask_coverage == 0).sum())
    summary["zero_reference_parcel_runs"] = int((df.nonzero_reference_coverage == 0).sum())
    summary["parcel_runs_below_half_mask_coverage"] = int((df.mask_coverage < .5).sum())
    summary["mask_coverage_minimum"] = float(df.mask_coverage.min())
    summary["input_sha256"] = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    (out / "review_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "input_sha256"}, indent=2))


if __name__ == "__main__":
    main()
