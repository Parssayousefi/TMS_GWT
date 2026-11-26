"""Compare four sampled sub-10 BOLD volumes with enhanced/uniform references.

Only volumes 1 and 150 of pre90/post90 are read. This distinguishes reference
enhancement from 4D intensity processing; it is not a temporal-quality audit.
"""
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
from aachen_alignment_reference_intensity import NEW, OUT, sha


def main():
    inputs, rows = {}, []

    def track(path):
        inputs[str(path)] = sha(path)
        return path

    anat = NEW / "sub-10/ses-anat/anat"
    a2m_path = track(anat / "sub-10_ses-anat_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5")
    a2m = audit.sitk.ReadTransform(str(a2m_path))
    t1mask_path = track(anat / "sub-10_ses-anat_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz")
    t1mask = np.asarray(nib.load(t1mask_path).dataobj) > 0
    for session in ["pre90", "post90"]:
        stem = f"sub-10_ses-{session}_task-rest"
        func = NEW / f"sub-10/ses-{session}/func"
        hmc_path = track(func / f"{stem}_desc-hmc_boldref.nii.gz")
        coreg_path = track(func / f"{stem}_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt")
        standard_path = track(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz")
        mask_path = track(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz")
        bold_path = track(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz")
        standard = audit.sitk.ReadImage(str(standard_path))
        ref = nib.load(standard_path)
        enhanced = np.asarray(ref.dataobj, dtype=np.float32)
        transform = audit.sitk.CompositeTransform([audit.sitk.ReadTransform(str(coreg_path)), a2m])
        hmc_mni = audit.array(audit.sitk.Resample(audit.sitk.ReadImage(str(hmc_path)), standard,
            transform, audit.sitk.sitkLinear, 0, audit.sitk.sitkFloat32))
        support = (np.asarray(nib.load(mask_path).dataobj) > 0) & t1mask
        support = ndimage.binary_erosion(support, iterations=1)
        support &= hmc_mni > np.percentile(hmc_mni[support], 10)
        proxy = nib.load(bold_path)
        stages = [("HMC reference through corrected geometry", hmc_mni),
                  ("Saved enhanced MNI reference", enhanced),
                  ("Actual preprocessed BOLD volume 1", np.asarray(proxy.dataobj[..., 0], dtype=np.float32)),
                  ("Actual preprocessed BOLD volume 150", np.asarray(proxy.dataobj[..., 149], dtype=np.float32))]
        center = np.rint(nib.affines.apply_affine(np.linalg.inv(ref.affine), [0,0,25])).astype(int)
        fig, axes = plt.subplots(4, 3, figsize=(12, 13), layout="constrained")
        for row, (label, data) in enumerate(stages):
            median = float(np.median(data[support]))
            gain = (data[support]/median) / (hmc_mni[support]/np.median(hmc_mni[support]))
            rows.append({"subject": "10", "session": session, "image": label,
                "r_to_unenhanced_hmc_through_corrected_geometry": audit.correlation(data, hmc_mni, support),
                "r_to_saved_enhanced_reference": audit.correlation(data, enhanced, support),
                "relative_gain_p05": float(np.percentile(gain,5)),
                "relative_gain_p50": float(np.percentile(gain,50)),
                "relative_gain_p95": float(np.percentile(gain,95)),
                "support_voxels": int(support.sum())})
            for axis in range(3):
                axes[row,axis].imshow(np.take(data/median, center[axis], axis=axis).T,
                    origin="lower", cmap="gray", vmin=0, vmax=2.5)
                plane = np.take(t1mask, center[axis], axis=axis).T
                axes[row,axis].contour(plane, levels=[.5], colors=["#39c5dc"], linewidths=.7)
                axes[row,axis].set_title(f"{label}\naxis {axis}", fontsize=9)
                axes[row,axis].axis("off")
        fig.suptitle(f"sub-10 {session}: actual 4D samples versus reference enhancement\nFixed cuts x=0, y=0, z=25 mm; median-scaled range 0–2.5; cyan = corrected T1 boundary")
        fig.savefig(OUT / f"sub-10_{session}_sampled_processed_bold.png", dpi=115)
        plt.close(fig)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "sampled_bold_reference_agreement.csv", index=False)
    work = audit.ROOT / "tmp/fmriprep_aachen_skullstrip_skip_20261125/sub-10/fmriprep_23_2_wf/sub_10_wf/bold_ses_pre90_task_rest_wf"
    reports = {}
    for relative in ["bold_fit_wf/enhance_and_skullstrip_bold_wf/n4_correct/_report/report.rst",
                     "bold_native_wf/bold_stc_wf/slice_timing_correction/_report/report.rst",
                     "bold_std_wf/_in_tuple_MNI152NLin2009cAsym.res2/resample/_report/report.rst",
                     "ds_bold_std_wf/_in_tuple_MNI152NLin2009cAsym.res2/ds_bold/_report/report.rst"]:
        path = track(work / relative)
        text = path.read_text()
        reports[relative] = {"path": str(path), "sha256": inputs[str(path)],
            "relevant_lines": [line for line in text.splitlines() if line.startswith("* ") and
                any(field in line for field in ["in_file :", "out_file :", "input_image :", "output_image :",
                    "transforms :", "weight_image :", "fieldmap :", "cmdline :"])]}
    n4_path = track(work / "bold_fit_wf/enhance_and_skullstrip_bold_wf/n4_correct/sub-10_ses-pre90_task-rest_bold_average_corrected.nii.gz")
    native_coreg_path = track(NEW / "sub-10/ses-pre90/func/sub-10_ses-pre90_task-rest_desc-coreg_boldref.nii.gz")
    n4, coreg = np.asarray(nib.load(n4_path).dataobj), np.asarray(nib.load(native_coreg_path).dataobj)
    n4_agreement = float(np.corrcoef(n4.ravel(), coreg.ravel())[0,1])
    (OUT / "enhancement_metadata.json").write_text(json.dumps({"reports": reports,
        "native_saved_coreg_to_n4_output_r": n4_agreement,
        "dataflow_observation": "The 4D resampling in_file is slice-time-corrected raw BOLD. Its transforms are HMC, coregistration and anatomical normalization; no N4-corrected intensity image is used as 4D in_file."}, indent=2), encoding="utf-8")
    outputs = list(OUT.glob("*sampled_processed_bold.png")) + [OUT / "sampled_bold_reference_agreement.csv", OUT / "enhancement_metadata.json"]
    (OUT / "sampled_bold_provenance.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": sha(Path(__file__)),
        "input_sha256": inputs, "output_sha256": {str(p):sha(p) for p in outputs},
        "processed_volumes_read": {"pre90":[1,150],"post90":[1,150]},
        "scope": "Reference-versus-4D intensity comparison on four sampled volumes only; no full temporal statistics or effects.",
        "method": "Resample the unchanged HMC reference through corrected fMRIPrep geometry. Compare its intensity pattern with saved enhanced reference and actual processed samples, within one-voxel-eroded common T1/BOLD support excluding the lowest HMC intensity decile.",
        "production_inputs_modified": False, "exclusions_changed": False}, indent=2), encoding="utf-8")
    print(frame.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
