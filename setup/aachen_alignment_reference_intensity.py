"""Bounded sub-10 intensity check: two raw volumes and native reference images.

No temporal variance, connectivity, effects, or exclusion changes are computed.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
import aachen_alignment_anatomical_review as audit

ROOT = audit.ROOT
OLD = ROOT / "data/aachen_itbs/derivatives/fmriprep"
NEW = ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125"
OUT = audit.REVIEW / "corrected_cohort/sub-10/reference_intensity"


def sha(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    inputs, rows, agreement = {}, [], []

    def track(path):
        if str(path) not in inputs:
            inputs[str(path)] = sha(path)
        return path

    for session in ["pre90", "post90"]:
        stem = f"sub-10_ses-{session}_task-rest"
        raw_path = track(ROOT / f"data/aachen_itbs/sub-10/ses-{session}/func/{stem}_bold.nii.gz")
        raw = nib.load(raw_path)
        stages = [("Raw volume 1", np.asarray(raw.dataobj[..., 0], dtype=np.float32)),
                  ("Raw volume 150", np.asarray(raw.dataobj[..., 149], dtype=np.float32))]
        reference_paths = []
        for root, label, suffix in [(OLD, "Original HMC reference", "hmc"),
                                    (NEW, "Corrected HMC reference", "hmc"),
                                    (OLD, "Original coreg reference", "coreg"),
                                    (NEW, "Corrected coreg reference", "coreg")]:
            path = track(root / f"sub-10/ses-{session}/func/{stem}_desc-{suffix}_boldref.nii.gz")
            image = nib.load(path)
            assert image.shape == raw.shape[:3]
            np.testing.assert_allclose(image.affine, raw.affine, atol=1e-4)
            stages.append((label, np.asarray(image.dataobj, dtype=np.float32)))
            reference_paths.append(path)
        native = audit.sitk.ReadImage(str(reference_paths[-1]))
        t1mask_path = track(NEW / "sub-10/ses-anat/anat/sub-10_ses-anat_desc-brain_mask.nii.gz")
        coreg_path = track(NEW / f"sub-10/ses-{session}/func/{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt")
        transform = audit.sitk.ReadTransform(str(coreg_path))
        support = audit.array(audit.sitk.Resample(audit.sitk.ReadImage(str(t1mask_path)), native,
            transform.GetInverse(), audit.sitk.sitkNearestNeighbor, 0, audit.sitk.sitkUInt8)) > 0
        interior = ndimage.binary_erosion(support, iterations=2)
        grid = np.stack(np.indices(interior.shape), axis=-1)
        y = nib.affines.apply_affine(raw.affine, grid)[..., 1]
        lo, hi = np.quantile(y[interior], [1/3, 2/3])
        anterior, posterior = interior & (y >= hi), interior & (y <= lo)
        center = np.rint(np.argwhere(interior).mean(0)).astype(int)
        fig, axes = plt.subplots(6, 3, figsize=(11, 18), layout="constrained")
        for row, (label, data) in enumerate(stages):
            a, p = float(np.median(data[anterior])), float(np.median(data[posterior]))
            median = float(np.median(data[interior]))
            ratio = a/p
            rows.append({"subject": "10", "session": session, "image": label,
                         "anterior_median": a, "posterior_median": p,
                         "anterior_posterior_median_ratio": ratio,
                         "interior_median": median})
            for axis in range(3):
                other = [i for i in range(3) if i != axis]
                aspect = raw.header.get_zooms()[other[1]] / raw.header.get_zooms()[other[0]]
                axes[row, axis].imshow(np.take(data/median, center[axis], axis=axis).T,
                    origin="lower", cmap="gray", vmin=0, vmax=2.5, aspect=aspect)
                axes[row, axis].set_title(f"{label} | axis {axis}\nanterior/posterior median = {ratio:.2f}", fontsize=9)
                axes[row, axis].axis("off")
        fig.suptitle(f"sub-10 {session}: raw and original/corrected native references\nSame native voxel slices; each image scaled by its interior-brain median; fixed display range 0–2.5")
        fig.savefig(OUT / f"sub-10_{session}_raw_reference_comparison.png", dpi=115)
        plt.close(fig)
        agreement.append({"session": session,
            "original_corrected_hmc_r": audit.correlation(stages[2][1], stages[3][1], interior),
            "original_corrected_coreg_r": audit.correlation(stages[4][1], stages[5][1], interior),
            "raw_volume150_corrected_hmc_r": audit.correlation(stages[1][1], stages[3][1], interior),
            "original_corrected_hmc_arrays_identical": bool(np.array_equal(stages[2][1], stages[3][1])),
            "original_corrected_coreg_arrays_identical": bool(np.array_equal(stages[4][1], stages[5][1]))})
    pd.DataFrame(rows).to_csv(OUT / "regional_intensity_ratios.csv", index=False)
    pd.DataFrame(agreement).to_csv(OUT / "reference_agreement.csv", index=False)
    paths = sorted(OUT.glob("*.png")) + sorted(OUT.glob("*.csv"))
    provenance = {"created_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": sha(Path(__file__)),
        "input_sha256": inputs, "output_sha256": {str(p):sha(p) for p in paths},
        "raw_volumes_read": {"pre90": [1,150], "post90": [1,150]},
        "method": "Anterior and posterior thirds of physical Y within a two-voxel-eroded corrected T1 support mask mapped onto the native BOLD grid. Per-image median scaling for display. No raw-image realignment was added.",
        "limitations": "A bounded spatial intensity comparison, not a temporal signal-quality test or diagnosis of the physical acquisition cause.",
        "exclusions_or_processing_changed": False, "functional_effects_read": False}
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(pd.DataFrame(rows)[["session","image","anterior_posterior_median_ratio"]].to_string(index=False), flush=True)
    print(pd.DataFrame(agreement).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
