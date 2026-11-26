"""Inspect pinned ds001927 inputs without changing source images.

Full T1 intensities plus the first 24 frames of each BOLD are inspected.
Contact sheets support visual skull/defacing review, not automated segmentation
or motion QC. No exclusions or treatment statistics are derived here.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ("RepetitionTime", "EchoTime", "SliceTiming", "PhaseEncodingDirection",
          "TotalReadoutTime", "EffectiveEchoSpacing", "SliceEncodingDirection",
          "MultibandAccelerationFactor", "NonSteadyStateOutlierIndices")


def panel(ax, array, voxels, axis, title):
    index = array.shape[axis] // 2
    sl = np.take(array, index, axis=axis).T
    plane = [i for i in range(3) if i != axis]
    values = array[array > 0]
    vmax = float(np.percentile(values, 99)) if values.size else 1
    ax.imshow(sl, origin="lower", cmap="gray", vmin=0, vmax=vmax,
              aspect=voxels[plane[1]] / voxels[plane[0]])
    ax.set_title(title, fontsize=8)
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "output/ds001927_input_audit_20261126")
    args = ap.parse_args()
    out = args.out.resolve()
    out.relative_to(ROOT)
    out.mkdir(parents=True, exist_ok=True)
    raw = ROOT / "data/ds001927"
    subjects = sorted(p.name for p in raw.glob("sub-*") if p.is_dir())
    assert len(subjects) == 23
    records, metadata = [], []
    anat_thumbs = {}
    bold_thumbs = {}
    detail_subjects = {subjects[0], subjects[len(subjects) // 2], subjects[-1]}
    details = []
    for subject in subjects:
        for session in ("ctr", "fro", "occ"):
            files = [raw / subject / f"ses-{session}" / "anat" / f"{subject}_ses-{session}_T1w.nii.gz"]
            files += [raw / subject / f"ses-{session}" / "func" /
                      f"{subject}_ses-{session}_task-rest{time}_bold.nii.gz"
                      for time in ("pre", "post")]
            for path in files:
                img = nib.load(path)
                is_bold = len(img.shape) == 4
                js = json.loads(path.with_name(path.name.replace(".nii.gz", ".json")).read_text())
                orientation = "".join(nib.aff2axcodes(img.affine))
                rec = dict(path=path.relative_to(raw).as_posix(), subject=subject,
                           session=session, modality="bold" if is_bold else "T1w",
                           shape=list(img.shape), voxel_sizes=list(map(float, img.header.get_zooms())),
                           orientation=orientation, qform_code=int(img.header["qform_code"]),
                           sform_code=int(img.header["sform_code"]))
                if int(img.header["qform_code"]) and int(img.header["sform_code"]):
                    rec["qform_sform_max_abs_mm"] = float(np.max(np.abs(img.get_qform() - img.get_sform())))
                if is_bold:
                    assert img.shape == (64, 64, 40, 600), (path, img.shape)
                    assert abs(float(img.header.get_zooms()[3]) - 1.25) < 1e-6
                    data = np.asarray(img.dataobj[..., :24], dtype=np.float32)
                    reference = np.median(data[..., 16:24], axis=-1)
                    rec["inspected_frames"] = [0, 24]
                    rec["early_global_means"] = np.mean(data, axis=(0, 1, 2)).astype(float).tolist()
                    metadata.append({"path": rec["path"], **{k: js.get(k) for k in FIELDS}})
                else:
                    data = img.get_fdata(dtype=np.float32)
                    reference = data
                    rec["inspected_frames"] = "complete structural volume"
                rec["finite_fraction_inspected"] = float(np.isfinite(data).mean())
                rec["zero_fraction_inspected"] = float((data == 0).mean())
                assert np.isfinite(data).all(), path
                canon = nib.as_closest_canonical(nib.Nifti1Image(reference, img.affine))
                vol = canon.get_fdata(dtype=np.float32)
                zooms = canon.header.get_zooms()
                # Keep only a central sagittal slice for whole-cohort contact sheets.
                thumb = (vol[vol.shape[0] // 2, :, :].T.copy(), float(zooms[2] / zooms[1]),
                         float(np.percentile(vol[vol > 0], 99)))
                if is_bold:
                    time = "pre" if "restpre" in path.name else "post"
                    bold_thumbs[(subject, session, time)] = thumb
                else:
                    anat_thumbs[(subject, session)] = thumb
                    if subject in detail_subjects:
                        details.append((subject, session, vol, zooms))
                records.append(rec)
                del data, reference, vol
        print(f"Inspected {subject}: 3 T1w and 6 BOLD early-frame samples", flush=True)

    for start in range(0, len(subjects), 8):
        subset = subjects[start:start + 8]
        for kind in ("T1w", "BOLD"):
            cols = 3 if kind == "T1w" else 6
            fig, axes = plt.subplots(len(subset), cols, figsize=(3.1 * cols, 2.7 * len(subset)), squeeze=False)
            for i, subject in enumerate(subset):
                for j, session in enumerate(("ctr", "fro", "occ")):
                    for t, time in enumerate(("pre", "post") if kind == "BOLD" else ("",)):
                        col = j * 2 + t if kind == "BOLD" else j
                        key = (subject, session, time) if kind == "BOLD" else (subject, session)
                        sl, aspect, vmax = (bold_thumbs if kind == "BOLD" else anat_thumbs)[key]
                        ax = axes[i, col]
                        ax.imshow(sl, origin="lower", cmap="gray", vmin=0, vmax=vmax, aspect=aspect)
                        ax.set_title(f"{subject} {session} {time}", fontsize=8)
                        ax.axis("off")
            fig.suptitle(f"ds001927 raw {kind}: central sagittal slices in RAS orientation\n"
                         "Visual input check; BOLD is median of frames 16–23, not motion-corrected", fontsize=12)
            fig.tight_layout(rect=(0, 0, 1, .97))
            fig.savefig(out / f"{kind.lower()}_contact_{start // 8 + 1}.png", dpi=100)
            plt.close(fig)
    fig, axes = plt.subplots(len(details), 3, figsize=(12, 3.2 * len(details)))
    for row, (subject, session, vol, zooms) in enumerate(details):
        for axis in range(3):
            panel(axes[row, axis], vol, zooms, axis, f"{subject} {session} {('sagittal', 'coronal', 'axial')[axis]}")
    fig.suptitle("Raw T1w images: orthogonal central sections; no extraction applied", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .98))
    fig.savefig(out / "t1w_orthogonal_examples.png", dpi=110)
    plt.close(fig)

    # Digest groups detect deposited duplicate anatomies without rehashing 23 GB.
    verified = json.loads((ROOT / "output/ds001927_download_20261126/verified_files.json").read_text())
    anat_hashes = [v["sha256"] for v in verified if v["path"].endswith("_T1w.nii.gz")]
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), subjects=len(subjects),
                   structural_images=len(anat_thumbs), bold_runs=len(bold_thumbs),
                   unique_structural_sha256=len(set(anat_hashes)),
                   missing_metadata_counts={k: sum(m[k] is None for m in metadata) for k in FIELDS},
                   fieldmap_files=len(list(raw.glob("sub-*/ses-*/fmap/*"))),
                   all_inspected_voxels_finite=all(r["finite_fraction_inspected"] == 1 for r in records),
                   orientations=dict(Counter(r["orientation"] for r in records)),
                   limitations=["BOLD intensities sampled only in first 24 frames; whole-file checksums were verified separately.",
                                "No head-motion, registration, atlas coverage or treatment-effect analysis performed.",
                                "Zero fraction and contact sheets are not an automated skull-stripping classifier."])
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (out / "image_checks.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    (out / "bold_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    pd.DataFrame([{k: v for k, v in r.items() if k != "early_global_means"} for r in records]).to_csv(out / "image_checks.csv", index=False)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
