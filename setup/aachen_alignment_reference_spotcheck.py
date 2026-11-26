"""Bounded raw/reference/processed spatial check for a flagged cohort subject.

No temporal statistics, effects, preprocessing or exclusion changes. Reads only
volumes 1 and 150 of pre90/post90; writes within the anatomical review folder.
"""
from argparse import ArgumentParser
from datetime import datetime, timezone
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
from aachen_alignment_reference_intensity import NEW, OLD, sha


def main():
    parser = ArgumentParser()
    parser.add_argument("--subject", required=True)
    subject = parser.parse_args().subject.zfill(2)
    queue = json.loads((audit.ROOT / "logs/aachen_corrected_20261125/queue_status.json").read_text())
    item = queue.get("finished", {}).get(subject, {})
    if not (item.get("success") is True and item.get("state", {}).get("ExitCode") == 0):
        raise ValueError("Subject processing has not completed successfully")
    folder = audit.REVIEW / f"corrected_cohort/sub-{subject}"
    if not (folder / "metrics.json").exists():
        raise ValueError("Main anatomical audit is required first")
    out = folder / "reference_intensity"
    out.mkdir(exist_ok=True)
    if (out / "reference_spotcheck_provenance.json").exists():
        raise FileExistsError("Refusing to replace a completed spot-check")
    inputs, native_rows, processed_rows = {}, [], []

    def track(path):
        path = Path(path)
        if str(path) not in inputs:
            inputs[str(path)] = sha(path)
        return path

    anat = NEW / f"sub-{subject}/ses-anat/anat"
    prefix = f"sub-{subject}_ses-anat"
    t1_native = audit.sitk.ReadImage(str(track(anat / f"{prefix}_desc-brain_mask.nii.gz")))
    t1_mni = np.asarray(nib.load(track(anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz")).dataobj) > 0
    a2m = audit.sitk.ReadTransform(str(track(anat / f"{prefix}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5")))
    for session in ["pre90", "post90"]:
        stem = f"sub-{subject}_ses-{session}_task-rest"
        func = NEW / f"sub-{subject}/ses-{session}/func"
        raw = nib.load(track(audit.ROOT / f"data/aachen_itbs/sub-{subject}/ses-{session}/func/{stem}_bold.nii.gz"))
        arrays = {f"Raw volume {v}": np.asarray(raw.dataobj[..., v-1], dtype=np.float32) for v in [1, 150]}
        for root, label in [(OLD, "Original"), (NEW, "Corrected")]:
            for kind in ["hmc", "coreg"]:
                image = nib.load(track(root / f"sub-{subject}/ses-{session}/func/{stem}_desc-{kind}_boldref.nii.gz"))
                np.testing.assert_allclose(image.affine, raw.affine, atol=1e-4)
                arrays[f"{label} {kind} reference"] = np.asarray(image.dataobj, dtype=np.float32)
        native = audit.sitk.ReadImage(str(func / f"{stem}_desc-coreg_boldref.nii.gz"))
        coreg = audit.sitk.ReadTransform(str(track(func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt")))
        support = ndimage.binary_erosion(audit.array(audit.sitk.Resample(t1_native, native, coreg.GetInverse(),
            audit.sitk.sitkNearestNeighbor, 0, audit.sitk.sitkUInt8)) > 0, iterations=2)
        center = np.rint(np.argwhere(support).mean(0)).astype(int)
        fig, axes = plt.subplots(6, 3, figsize=(11, 18), layout="constrained")
        for row, (label, data) in enumerate(arrays.items()):
            normalized = data / np.median(data[support])
            for axis in range(3):
                other = [i for i in range(3) if i != axis]
                aspect = raw.header.get_zooms()[other[1]] / raw.header.get_zooms()[other[0]]
                axes[row, axis].imshow(np.take(normalized, center[axis], axis=axis).T,
                    origin="lower", cmap="gray", vmin=0, vmax=2.5, aspect=aspect)
                axes[row, axis].set_title(f"{label} | axis {axis}", fontsize=9)
                axes[row, axis].axis("off")
        fig.suptitle(f"sub-{subject} {session}: raw and original/corrected references\nSame native cuts; interior-median scaled; fixed range 0–2.5")
        fig.savefig(out / f"sub-{subject}_{session}_raw_reference_comparison.png", dpi=115)
        plt.close(fig)
        native_rows.append({"subject": subject, "session": session,
            "original_corrected_hmc_arrays_identical": bool(np.array_equal(arrays["Original hmc reference"], arrays["Corrected hmc reference"])),
            "original_corrected_coreg_r": audit.correlation(arrays["Original coreg reference"], arrays["Corrected coreg reference"], support),
            "raw_volume1_corrected_hmc_r": audit.correlation(arrays["Raw volume 1"], arrays["Corrected hmc reference"], support),
            "raw_volume150_corrected_hmc_r": audit.correlation(arrays["Raw volume 150"], arrays["Corrected hmc reference"], support)})

        reference_path = track(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz")
        reference = audit.sitk.ReadImage(str(reference_path))
        ref_nib = nib.load(reference_path)
        enhanced = np.asarray(ref_nib.dataobj)
        hmc_mni = audit.array(audit.sitk.Resample(audit.sitk.ReadImage(str(func / f"{stem}_desc-hmc_boldref.nii.gz")),
            reference, audit.sitk.CompositeTransform([coreg, a2m]), audit.sitk.sitkLinear, 0, audit.sitk.sitkFloat32))
        brain = np.asarray(nib.load(track(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz")).dataobj) > 0
        support = ndimage.binary_erosion(brain & t1_mni, iterations=1)
        support &= hmc_mni > np.percentile(hmc_mni[support], 10)
        bold = nib.load(track(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz"))
        stages = [("HMC reference through corrected geometry", hmc_mni), ("Saved enhanced MNI reference", enhanced)]
        stages += [(f"Actual processed volume {v}", np.asarray(bold.dataobj[..., v-1], dtype=np.float32)) for v in [1, 150]]
        center = np.rint(nib.affines.apply_affine(np.linalg.inv(ref_nib.affine), [0, 0, 25])).astype(int)
        fig, axes = plt.subplots(4, 3, figsize=(12, 13), layout="constrained")
        for row, (label, data) in enumerate(stages):
            normalized = data / np.median(data[support])
            gain = normalized[support] / (hmc_mni[support] / np.median(hmc_mni[support]))
            processed_rows.append({"subject": subject, "session": session, "image": label,
                "r_to_unenhanced_hmc_through_corrected_geometry": audit.correlation(data, hmc_mni, support),
                "r_to_saved_enhanced_reference": audit.correlation(data, enhanced, support),
                "relative_gain_p05": float(np.percentile(gain, 5)), "relative_gain_p95": float(np.percentile(gain, 95))})
            for axis in range(3):
                axes[row, axis].imshow(np.take(normalized, center[axis], axis=axis).T, origin="lower", cmap="gray", vmin=0, vmax=2.5)
                axes[row, axis].contour(np.take(t1_mni, center[axis], axis=axis).T, levels=[.5], colors=["#39c5dc"], linewidths=.7)
                axes[row, axis].set_title(f"{label}\naxis {axis}", fontsize=9)
                axes[row, axis].axis("off")
        fig.suptitle(f"sub-{subject} {session}: actual BOLD samples versus enhanced reference\nx=0, y=0, z=25 mm; median-scaled range 0–2.5; cyan = corrected T1 boundary")
        fig.savefig(out / f"sub-{subject}_{session}_sampled_processed_bold.png", dpi=115)
        plt.close(fig)

    work = audit.ROOT / f"tmp/fmriprep_aachen_skullstrip_skip_20261125/sub-{subject}/fmriprep_23_2_wf/sub_{subject}_wf/bold_ses_pre90_task_rest_wf"
    reports = {}
    for relative in ["bold_fit_wf/enhance_and_skullstrip_bold_wf/n4_correct/_report/report.rst",
                     "bold_native_wf/bold_stc_wf/slice_timing_correction/_report/report.rst",
                     "bold_std_wf/_in_tuple_MNI152NLin2009cAsym.res2/resample/_report/report.rst",
                     "ds_bold_std_wf/_in_tuple_MNI152NLin2009cAsym.res2/ds_bold/_report/report.rst"]:
        path = track(work / relative)
        reports[relative] = {"path": str(path), "sha256": inputs[str(path)], "relevant_lines": [
            line for line in path.read_text().splitlines() if line.startswith("* ") and any(
                key in line for key in ["in_file :", "out_file :", "input_image :", "output_image :", "transforms :", "weight_image :", "fieldmap :", "cmdline :"])]}
    n4 = nib.load(track(work / f"bold_fit_wf/enhance_and_skullstrip_bold_wf/n4_correct/sub-{subject}_ses-pre90_task-rest_bold_average_corrected.nii.gz"))
    saved = nib.load(NEW / f"sub-{subject}/ses-pre90/func/sub-{subject}_ses-pre90_task-rest_desc-coreg_boldref.nii.gz")
    metadata = {"reports": reports, "native_saved_coreg_to_n4_output_r": float(np.corrcoef(np.asarray(n4.dataobj).ravel(), np.asarray(saved.dataobj).ravel())[0, 1])}
    (out / "enhancement_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    pd.DataFrame(native_rows).to_csv(out / "reference_agreement.csv", index=False)
    pd.DataFrame(processed_rows).to_csv(out / "sampled_bold_reference_agreement.csv", index=False)
    outputs = list(out.glob("*.png")) + [out / "enhancement_metadata.json", out / "reference_agreement.csv", out / "sampled_bold_reference_agreement.csv"]
    provenance = {"subject": subject, "created_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": sha(Path(__file__)),
        "input_sha256": inputs, "output_sha256": {str(p): sha(p) for p in outputs},
        "raw_and_processed_volumes_read": {"pre90": [1, 150], "post90": [1, 150]},
        "scope": "Bounded reference-versus-4D intensity and anatomical sample check, not temporal quality or effects.",
        "production_inputs_modified": False, "exclusions_changed": False}
    (out / "reference_spotcheck_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(pd.DataFrame(native_rows).to_string(index=False), flush=True)
    print(pd.DataFrame(processed_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
