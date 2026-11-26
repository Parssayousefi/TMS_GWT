"""Localize saved SUDMEX mask disagreements; report spatial coverage, not effects."""
import json
import re
import nibabel as nib
import numpy as np
import pandas as pd
from audit_sudmex_registration import ROOT, OUT, DERIV, SPACE, one, on_grid


def main():
    runs = pd.read_csv(OUT / "alignment_metrics.csv", dtype={"subject": str})
    motion = pd.read_csv(ROOT / "pipeline/results/qc/sudmex_tms/exclusions.csv", dtype={"subject": str})
    motion["subject"] = motion.subject.str.zfill(3)
    motion = motion.rename(columns={"status": "existing_motion_status"})
    atlas = nib.load(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    rows, parcel_rows, report_rows = [], [], []
    for subject, group in runs.groupby("subject", sort=True):
        tm_img = nib.load(one(DERIV / f"sub-{subject}/ses-t0/anat", f"*space-{SPACE}*desc-brain_mask.nii.gz"))
        for run in group.to_dict("records"):
            b_img = nib.load(ROOT / run["MNI_BOLD_mask_path"])
            b = np.asarray(b_img.dataobj) > 0
            tm = on_grid(tm_img, b_img) > 0
            ref = on_grid(nib.load(ROOT / run["MNI_BOLD_reference_path"]), b_img, 1)
            labels = np.rint(on_grid(atlas, b_img)).astype(np.int16)
            cortex = (labels > 0) & (labels <= 400) & tm
            subcortex = (labels > 400) & tm
            z = (b_img.affine[2, 0] * np.arange(b.shape[0])[:, None, None]
                 + b_img.affine[2, 1] * np.arange(b.shape[1])[None, :, None]
                 + b_img.affine[2, 2] * np.arange(b.shape[2])[None, None, :]
                 + b_img.affine[2, 3])
            low, high = tm & (z < -30), tm & (z >= 0)
            missing = tm & ~b
            denom = np.bincount(labels[tm], minlength=455)[1:455]
            numerator = np.bincount(labels[tm & b], minlength=455)[1:455]
            coverage = np.divide(numerator, denom, out=np.full(454, np.nan), where=denom > 0)
            rows.append({
                "subject": subject, "session": run["session"],
                "anatomical_cortical_atlas_in_bold_mask_fraction": float((cortex & b).sum() / cortex.sum()),
                "anatomical_subcortical_atlas_in_bold_mask_fraction": float((subcortex & b).sum() / subcortex.sum()),
                "anatomical_mask_below_z_minus30_in_bold_mask_fraction": float((low & b).sum() / low.sum()),
                "anatomical_mask_at_or_above_z0_in_bold_mask_fraction": float((high & b).sum() / high.sum()),
                "anatomical_mask_in_positive_boldref_fraction": float((tm & (ref > 0)).sum() / tm.sum()),
                "missing_anatomical_mask_fraction_with_zero_boldref": float((missing & (ref <= 0)).sum() / missing.sum()),
                "cortical_parcels_with_less_than_half_anatomical_support_in_bold_mask": int((coverage[:400] < .5).sum()),
                "subcortical_parcels_with_less_than_half_anatomical_support_in_bold_mask": int((coverage[400:] < .5).sum()),
            })
            parcel_rows.extend({"subject": subject, "session": run["session"], "parcel_id": i + 1,
                                "anatomical_support_voxels": int(denom[i]),
                                "bold_mask_overlap_voxels": int(numerator[i]),
                                "fraction_of_anatomical_parcel_in_bold_mask": float(coverage[i])} for i in range(454))
            summary_path = DERIV / f"sub-{subject}/figures/sub-{subject}_ses-{run['session']}_task-rest_desc-summary_bold.html"
            summary = summary_path.read_text(encoding="utf-8")
            def extract(label):
                match = re.search(r"<li>" + re.escape(label) + r": (.*?)</li>", summary)
                return re.sub(r"<[^>]+>", "", match[1]) if match else "NOT_FOUND"
            report_rows.append({"subject": subject, "session": run["session"],
                                "susceptibility_distortion_correction": extract("Susceptibility distortion correction"),
                                "registration_method": extract("Registration"),
                                "summary_path": str(summary_path.relative_to(ROOT))})
    coverage = pd.DataFrame(rows)
    reports = pd.DataFrame(report_rows)
    merged = runs.merge(coverage, on=["subject", "session"], validate="one_to_one").merge(
        motion[["subject", "session", "existing_motion_status", "reason", "mean_fd_mm", "max_trans_mm", "max_rot_deg"]],
        on=["subject", "session"], validate="one_to_one").merge(reports, on=["subject", "session"], validate="one_to_one")
    assert len(merged) == 142
    merged.to_csv(OUT / "run_qc_with_coverage_and_motion.csv", index=False)
    parcels = pd.DataFrame(parcel_rows)
    parcels.to_csv(OUT / "parcel_mask_coverage.csv", index=False)
    names = pd.read_csv(ROOT / "pipeline/atlases/roi_to_network_mapping.csv")
    named = parcels.merge(names[["ROI", "ROI_Label", "Network_Name", "Hemisphere"]],
                          left_on="parcel_id", right_on="ROI", validate="many_to_one")
    named[named.fraction_of_anatomical_parcel_in_bold_mask < .5].to_csv(
        OUT / "poorly_covered_parcels_with_labels.csv", index=False)
    summary = {"motion_status_counts": merged.existing_motion_status.value_counts().to_dict(),
               "distortion_correction_counts": merged.susceptibility_distortion_correction.value_counts().to_dict(),
               "registration_method_counts": merged.registration_method.value_counts().to_dict(),
               "coverage_metrics": coverage.drop(columns=["subject", "session"]).describe().to_dict(),
               "limits": "Coverage describes saved masks; it is not signal quality, neural activity, or a validated exclusion threshold. z slabs are geometric, not tissue segmentations."}
    (OUT / "coverage_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "coverage_metrics"}, indent=2))
    print(coverage.drop(columns=["subject", "session"]).describe().to_string())


if __name__ == "__main__":
    main()
