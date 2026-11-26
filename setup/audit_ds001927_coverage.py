"""Trace pilot coverage from acquired BOLD through the saved spatial transforms.

This generates QC evidence only: no outcome estimates, new registration, masking
of production images, or automatic anatomical approval. Run prepare/summarize in
the gwt environment and warp inside the pinned fMRIPrep image with read-only
project input and a separately writable audit directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd

SPACE = "MNI152NLin2009cAsym_res-2"
ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def run_paths(root, subject):
    for session in ("ctr", "fro", "occ"):
        for task in ("restpre", "restpost"):
            stem = f"sub-{subject}_ses-{session}_task-{task}"
            raw = root / "data/ds001927" / f"sub-{subject}" / f"ses-{session}" / "func" / f"{stem}_bold.nii.gz"
            deriv = root / "data/ds001927/derivatives/fmriprep_23_2_1" / f"sub-{subject}"
            func = deriv / f"ses-{session}" / "func"
            yield session, task, stem, raw, deriv, func


def save_nifti(path, values, reference):
    nib.save(nib.Nifti1Image(values, reference.affine), path)


def prepare(root, out, subjects):
    records = []
    for subject in subjects:
        for session, task, stem, raw, deriv, func in run_paths(root, subject):
            im = nib.load(raw)
            assert im.shape == (64, 64, 40, 600)
            ref = nib.load(func / f"{stem}_desc-coreg_boldref.nii.gz")
            np.testing.assert_allclose(im.affine, ref.affine, atol=1e-5)
            assert ref.shape == im.shape[:3]
            # Read the full raw payload, then retain the declared initial trim.
            values = im.get_fdata(dtype=np.float32)
            assert np.isfinite(values).all(), raw
            retained = values[..., 8:]
            mean = retained.mean(axis=3, dtype=np.float64).astype(np.float32)
            sd = retained.std(axis=3, dtype=np.float64).astype(np.float32)
            save_nifti(out / f"{stem}_raw_mean.nii.gz", mean, ref)
            save_nifti(out / f"{stem}_raw_sd.nii.gz", sd, ref)
            save_nifti(out / f"{stem}_acquired_fov.nii.gz", np.ones(ref.shape, np.uint8), ref)
            records.append(dict(subject=subject, session=session, task=task,
                                raw_path=raw.relative_to(root).as_posix(), raw_sha256=sha(raw),
                                raw_shape=list(im.shape), raw_payload_all_finite=True,
                                reference_grid_matches_raw=True, retained_frames=592,
                                zero_temporal_variance_voxels=int((sd == 0).sum())))
            print(f"Prepared raw QC: {stem}", flush=True)
            del values, retained
    (out / "raw_payload_audit.json").write_text(json.dumps(records, indent=2), encoding="utf-8")


def warp(root, out, subjects):
    commands = []
    atlas = root / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz"
    for subject in subjects:
        for session, task, stem, raw, deriv, func in run_paths(root, subject):
            transforms = [deriv / "anat" / f"sub-{subject}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5",
                          func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt"]
            items = [(out / f"{stem}_raw_mean.nii.gz", "raw_mean", "Linear"),
                     (out / f"{stem}_acquired_fov.nii.gz", "acquired_fov", "NearestNeighbor"),
                     (func / f"{stem}_desc-coreg_boldref.nii.gz", "ref_parity", "LanczosWindowedSinc")]
            for source, label, interpolation in items:
                destination = out / f"{stem}_{label}_MNI.nii.gz"
                command = ["antsApplyTransforms", "--dimensionality", "3", "--float", "1",
                           "--input", str(source), "--reference-image", str(atlas),
                           "--output", str(destination), "--interpolation", interpolation,
                           "--default-value", "0"]
                for transform in transforms:
                    command.extend(["--transform", str(transform)])
                subprocess.run(command, check=True)
                commands.append(command)
            print(f"Applied saved transforms: {stem}", flush=True)
    (out / "transform_commands.json").write_text(json.dumps(commands, indent=2), encoding="utf-8")


def image_panel(ax, values, index, title, contours=()):
    positive = values[values > 0]
    ax.imshow(values[:, :, index].T, origin="lower", cmap="gray", vmin=0,
              vmax=np.percentile(positive, 98) if len(positive) else 1)
    for mask, color in contours:
        plane = mask[:, :, index].T
        if plane.min() < .5 < plane.max():
            ax.contour(plane, levels=[.5], colors=[color], linewidths=.6)
    ax.set_title(title, fontsize=9)
    ax.axis("off")


def summarize(root, out, subjects):
    atlas_im = nib.load(root / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    atlas = np.asarray(atlas_im.dataobj, dtype=int)
    counts = np.bincount(atlas.ravel(), minlength=455)[1:455]
    assert np.all(counts > 0) and atlas.max() == 454
    centroid_z = np.array([nib.affines.apply_affine(atlas_im.affine, np.argwhere(atlas == r).mean(axis=0))[2]
                           for r in range(1, 455)])
    template = nib.load(root / "tmp/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz").get_fdata() > 0
    zlow = int(round(nib.affines.apply_affine(np.linalg.inv(atlas_im.affine), [0, 0, -20])[2]))
    rows, run_records = [], []
    for subject in subjects:
        fig, axes = plt.subplots(6, 4, figsize=(14, 16))
        for row, (session, task, stem, raw, deriv, func) in enumerate(run_paths(root, subject)):
            def load(path):
                im = nib.load(path)
                np.testing.assert_allclose(im.affine, atlas_im.affine, atol=1e-5)
                assert im.shape == atlas_im.shape
                values = im.get_fdata(dtype=np.float32)
                assert np.isfinite(values).all()
                return values
            raw_mean = load(out / f"{stem}_raw_mean_MNI.nii.gz")
            fov = load(out / f"{stem}_acquired_fov_MNI.nii.gz") > .5
            ref_parity = load(out / f"{stem}_ref_parity_MNI.nii.gz")
            ref = load(func / f"{stem}_space-{SPACE}_boldref.nii.gz")
            mask = load(func / f"{stem}_space-{SPACE}_desc-brain_mask.nii.gz") > .5
            parity_error = float(abs(ref_parity - ref).max())
            assert parity_error < 1e-4, (stem, parity_error)
            scale = float(np.median(raw_mean[mask & (atlas > 0)]))
            assert scale > 0
            relative_raw = raw_mean / scale
            # Full processed payload check. Mean/std and support are QC only.
            bold_im = nib.load(func / f"{stem}_space-{SPACE}_desc-preproc_bold.nii.gz")
            np.testing.assert_allclose(bold_im.affine, atlas_im.affine, atol=1e-5)
            assert bold_im.shape == (*atlas_im.shape, 600)
            bold = bold_im.get_fdata(dtype=np.float32)
            assert np.isfinite(bold).all(), stem
            mean = bold[..., 8:].mean(axis=3, dtype=np.float64).astype(np.float32)
            sd = bold[..., 8:].std(axis=3, dtype=np.float64).astype(np.float32)
            del bold
            save_nifti(out / f"{stem}_preprocessed_mean.nii.gz", mean, atlas_im)
            save_nifti(out / f"{stem}_preprocessed_sd.nii.gz", sd, atlas_im)
            fractions = {}
            for label, flags in (("fov_coverage", fov), ("mask_coverage", mask),
                                 ("raw_signal_above_20pct", relative_raw > .2),
                                 ("nonzero_temporal_variance", sd > 0)):
                fractions[label] = np.bincount(atlas[flags], minlength=455)[1:455] / counts
            for roi in range(454):
                inroi = atlas == roi + 1
                rows.append(dict(subject=subject, session=session, task=task, roi=roi + 1,
                                 centroid_z_mm=centroid_z[roi],
                                 raw_relative_median=float(np.median(relative_raw[inroi])),
                                 **{key: value[roi] for key, value in fractions.items()}))
            missing = (atlas > 0) & ~mask
            record = dict(subject=subject, session=session, task=task, ref_transform_max_error=parity_error,
                          full_preprocessed_payload_finite=True,
                          geometric_fov_atlas_fraction=float(fov[atlas > 0].mean()),
                          atlas_voxels_outside_functional_mask=int(missing.sum()),
                          fraction_outside_mask_still_in_acquired_fov=float(fov[missing].mean()),
                          outside_mask_raw_relative_median=float(np.median(relative_raw[missing])),
                          outside_mask_processed_relative_median=float(np.median(mean[missing]) / np.median(mean[mask & (atlas > 0)])),
                          zero_temporal_variance_atlas_voxels=int(((sd == 0) & (atlas > 0)).sum()),
                          zero_variance_parcels=int((fractions["nonzero_temporal_variance"] == 0).sum()))
            run_records.append(record)
            raw_native = nib.load(out / f"{stem}_raw_mean.nii.gz").get_fdata(dtype=np.float32)
            # Native acquisition slice shown as acquired, without registration.
            image_panel(axes[row, 0], raw_native, 5, f"{session} {task}: raw slice 5")
            image_panel(axes[row, 1], raw_mean, zlow, "Raw mean, saved warp only", [(template, "#26bfcc")])
            image_panel(axes[row, 2], mean, zlow, "Preprocessed mean", [(template, "#26bfcc"), (mask, "#ffa64c")])
            image_panel(axes[row, 3], ref, zlow, "Reference + acquired FOV", [(fov, "#df68e3"), (mask, "#ffa64c")])
            print(f"Full payload / coverage checked: {stem}", flush=True)
        fig.suptitle(f"sub-{subject}: tracing inferior signal before and after preprocessing\nNative raw slice at left; remaining panels at MNI z=-20 mm. Raw mean is not motion-corrected.", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, .96))
        fig.savefig(out / f"sub-{subject}_raw_to_preprocessed_coverage.png", dpi=130)
        plt.close(fig)
    df = pd.DataFrame(rows)
    df.to_csv(out / "parcel_coverage_trace.csv", index=False)
    pd.DataFrame(run_records).to_csv(out / "run_coverage_trace.csv", index=False)
    low = df[df.mask_coverage < .5]
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), subjects=subjects,
                   run_count=len(run_records), below_half_mask_parcel_runs=len(low),
                   below_half_mask_parcel_runs_with_fov_above_95pct=int((low.fov_coverage > .95).sum()),
                   below_half_mask_raw_relative_median=float(low.raw_relative_median.median()),
                   below_half_mask_centroid_z_median=float(low.centroid_z_mm.median()),
                   zero_variance_parcel_runs=sum(x["zero_variance_parcels"] for x in run_records),
                   max_transform_parity_error=max(x["ref_transform_max_error"] for x in run_records),
                   raw_payload_receipt_sha256=sha(out / "raw_payload_audit.json"),
                   script_sha256=sha(Path(__file__)), analysis_release=False,
                   scope="Pilot QC; weak raw signal and finite variance do not prove valid neuronal time series. Raw mean uses all 592 retained frames without realignment; warped mean is a visual diagnostic, not an aligned voxelwise treatment comparison. The 20% signal display threshold is diagnostic only.")
    (out / "coverage_trace_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prepare", "warp", "summarize"])
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--subjects", nargs="+", choices=[f"s{i:03}" for i in range(1, 24)], default=["s001", "s022"])
    args = ap.parse_args()
    out = args.out or args.root / "output/ds001927_coverage_trace_20261126"
    out.mkdir(parents=True, exist_ok=True)
    globals()[args.stage](args.root, out, args.subjects)


if __name__ == "__main__":
    main()
