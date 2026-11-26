"""Read-only saved SUDMEX MNI anatomy/BOLD audit for the entire configured cohort."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import nibabel as nib
from nibabel.processing import resample_from_to
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/sudmex_registration_audit"
RAW = ROOT / "data/sudmex_tms"
DERIV = RAW / "derivatives/fmriprep"
SPACE = "MNI152NLin2009cAsym"
TEMPLATE = ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym"


def configured_cohort():
    cfg = yaml.safe_load((ROOT / "pipeline/configs/datasets/sudmex_tms.yaml").read_text())
    excluded = set(cfg["subjects"].get("exclude", [])) | set(cfg["qc"].get("exclude_subjects", []))
    included = cfg["subjects"].get("include", [])
    subjects = sorted(included or [p.name[4:] for p in RAW.glob("sub-*") if p.is_dir()])
    return [s for s in subjects if s not in excluded], cfg["qc"]["include_sessions"]


def one(folder, pattern):
    paths = sorted(folder.glob(pattern))
    if len(paths) != 1:
        raise ValueError(f"Expected one {folder / pattern}; found {len(paths)}")
    return paths[0]


def on_grid(image, reference, order=0):
    if image.shape == reference.shape and np.allclose(image.affine, reference.affine, atol=1e-4):
        return np.asarray(image.dataobj)
    return np.asarray(resample_from_to(image, reference, order=order).dataobj)


def dice(a, b):
    return float(2 * (a & b).sum() / (a.sum() + b.sum()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    subjects, sessions = configured_cohort()
    atlas = nib.load(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    template_mask = nib.load(TEMPLATE / "tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz")
    template_t1 = nib.load(TEMPLATE / "tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz")
    rois = {s: nib.load(ROOT / f"output/sudmex_target_anatomy/common_{s}_target_neighbourhood.nii.gz") for s in ["l", "r"]}
    rows, anatomical_rows, expected, display_rows = [], [], [], []
    for subject in subjects:
        folder = DERIV / f"sub-{subject}"
        anat = folder / "ses-t0/anat"
        t1mask_path = one(anat, f"*space-{SPACE}*desc-brain_mask.nii.gz")
        t1_path = one(anat, f"*space-{SPACE}*desc-preproc_T1w.nii.gz")
        t1mask, t1 = nib.load(t1mask_path), nib.load(t1_path)
        t = np.asarray(t1mask.dataobj) > 0
        tm = on_grid(template_mask, t1mask) > 0
        labels = on_grid(atlas, t1mask)
        cortex = (labels > 0) & (labels <= 400)
        whole_atlas = labels > 0
        t1data = on_grid(t1, t1mask, 1)
        template_data = on_grid(template_t1, t1mask, 1)
        anatomical_rows.append({
            "subject": subject, "MNI_anatomy_template_mask_dice": dice(t, tm),
            "MNI_anatomy_to_template_mask_volume_ratio": float(t.sum() / tm.sum()),
            "MNI_template_brain_in_anatomy_mask_fraction": float((t & tm).sum() / tm.sum()),
            "MNI_cortical_atlas_in_anatomy_mask_fraction": float((t & cortex).sum() / cortex.sum()),
            "MNI_whole_atlas_in_anatomy_mask_fraction": float((t & whole_atlas).sum() / whole_atlas.sum()),
            "MNI_anatomy_template_intensity_pearson_r": float(np.corrcoef(t1data[t & tm], template_data[t & tm])[0, 1]),
            "MNI_anatomy_path": str(t1_path.relative_to(ROOT)),
        })
        for session in sessions:
            raw_bolds = sorted((RAW / f"sub-{subject}/ses-{session}/func").glob("*_task-rest_bold.nii.gz"))
            if not raw_bolds:
                continue
            if len(raw_bolds) != 1:
                raise ValueError(f"Ambiguous raw BOLD: {raw_bolds}")
            expected.append((subject, session))
            func = folder / f"ses-{session}/func"
            boldmask_path = one(func, f"*space-{SPACE}*desc-brain_mask.nii.gz")
            boldref_path = one(func, f"*space-{SPACE}*boldref.nii.gz")
            bold_path = one(func, f"*space-{SPACE}*desc-preproc_bold.nii.gz")
            confounds_path = one(func, "*_desc-confounds_timeseries.tsv")
            boldmask, boldref = nib.load(boldmask_path), nib.load(boldref_path)
            b = np.asarray(boldmask.dataobj) > 0
            tb = on_grid(t1mask, boldmask) > 0
            template_b = on_grid(template_mask, boldmask) > 0
            atlas_b = on_grid(atlas, boldmask)
            cortex_b = (atlas_b > 0) & (atlas_b <= 400)
            atlas_support = atlas_b > 0
            bold_img, raw_bold_img = nib.load(bold_path), nib.load(raw_bolds[0])
            confounds = pd.read_csv(confounds_path, sep="\t")
            row = {
                "subject_id": f"sudmex_tms_{subject}", "subject": subject, "session": session,
                "bold_to_T1_mask_volume_ratio": float(b.sum() / tb.sum()), "dice": dice(b, tb),
                "T1_mask_in_bold_mask_fraction": float((b & tb).sum() / tb.sum()),
                "bold_mask_in_T1_mask_fraction": float((b & tb).sum() / b.sum()),
                "template_mask_dice": dice(b, template_b),
                "template_brain_in_bold_mask_fraction": float((b & template_b).sum() / template_b.sum()),
                "cortical_atlas_in_bold_brain_fraction": float((b & cortex_b).sum() / cortex_b.sum()),
                "whole_atlas_in_bold_brain_fraction": float((b & atlas_support).sum() / atlas_support.sum()),
                "bold_volumes": bold_img.shape[3], "raw_bold_volumes": raw_bold_img.shape[3],
                "confounds_rows": len(confounds),
                "boldref_nonfinite_voxels": int(np.count_nonzero(~np.isfinite(np.asarray(boldref.dataobj)))),
                "MNI_BOLD_path": str(bold_path.relative_to(ROOT)),
                "MNI_BOLD_mask_path": str(boldmask_path.relative_to(ROOT)),
                "MNI_BOLD_reference_path": str(boldref_path.relative_to(ROOT)),
                "confounds_path": str(confounds_path.relative_to(ROOT)),
            }
            for side, roi in rois.items():
                r = on_grid(roi, boldmask) > 0
                row[side + "_target_in_bold_brain_fraction"] = float((r & b).sum() / r.sum())
            rows.append(row)
            bref, bm = on_grid(boldref, t1mask, 1), on_grid(boldmask, t1mask) > 0
            indices = [(0, int(round(nib.affines.apply_affine(np.linalg.inv(t1mask.affine), [-32, 0, 0])[0]))),
                       (2, int(round(nib.affines.apply_affine(np.linalg.inv(t1mask.affine), [0, 0, 36])[2])))]
            planes = [(np.take(bref, idx, axis=ax), np.take(bm, idx, axis=ax),
                       np.take(t, idx, axis=ax), np.take(t1data, idx, axis=ax)) for ax, idx in indices]
            display_rows.append((subject, session, row["dice"], planes, float(np.percentile(bref[bref > 0], 98)),
                                 float(np.percentile(t1data[t1data > 0], 98))))
        print(f"MNI audit sub-{subject}: {len(rows)} runs", flush=True)
    d, a = pd.DataFrame(rows), pd.DataFrame(anatomical_rows)
    assert len(subjects) == 44, f"Configured participant count changed: {len(subjects)}"
    assert len(expected) == len(d) == 142, f"Expected/saved run count: {len(expected)}/{len(d)}"
    assert not d.duplicated(["subject", "session"]).any()
    assert (d.bold_volumes == d.raw_bold_volumes).all()
    assert (d.bold_volumes == d.confounds_rows).all()
    d.to_csv(OUT / "alignment_metrics.csv", index=False)
    a.to_csv(OUT / "standard_anatomy_metrics.csv", index=False)
    plot_dir = OUT / "standard_contact_sheets"
    plot_dir.mkdir(exist_ok=True)
    pages = []
    for start in range(0, len(display_rows), 6):
        part = display_rows[start:start + 6]
        fig, axes = plt.subplots(len(part), 4, figsize=(12, 2.7 * len(part)), squeeze=False, layout="constrained")
        for i, (subject, session, score, planes, bmax, tmax) in enumerate(part):
            for j, (bplane, bm, tm, tplane) in enumerate(planes):
                for k, (plane, vmax, title) in enumerate([(bplane, bmax, "BOLD"), (tplane, tmax, "T1")]):
                    ax = axes[i, j * 2 + k]
                    ax.imshow(plane.T, origin="lower", cmap="gray", vmin=0, vmax=vmax)
                    ax.contour(tm.T, levels=[.5], colors=["#50dbe8"], linewidths=.7)
                    ax.contour(bm.T, levels=[.5], colors=["#ff9f40"], linewidths=.7)
                    ax.set_title(f"{subject} {session} {title} {'x=-32' if j == 0 else 'z=36'} | Dice {score:.3f}", fontsize=8)
                    ax.axis("off")
        fig.suptitle("SUDMEX all saved runs | cyan: anatomical mask; orange: BOLD mask", fontsize=11)
        path = plot_dir / f"page_{start // 6 + 1:02d}.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        pages.append({"path": str(path.relative_to(OUT)), "runs": [f"{r[0]}_{r[1]}" for r in part]})
    summary = {"configured_subjects": len(subjects), "expected_raw_runs": len(expected), "complete_saved_runs": len(d),
               "session_counts": d.session.value_counts().to_dict(), "volume_and_confound_row_counts_match": True,
               "run_metrics": d.select_dtypes(include=np.number).describe().to_dict(),
               "standard_anatomy_metrics": a.select_dtypes(include=np.number).describe().to_dict(),
               "standard_contact_sheets": pages,
               "caution": "Mask agreement is descriptive QC, not proof of correct extraction, normalization, or undistorted signal."}
    (OUT / "alignment_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(d.select_dtypes(include=np.number).describe().to_string())
    print("Worst Dice:", d.sort_values("dice").head(10)[["subject", "session", "dice", "bold_to_T1_mask_volume_ratio"]].to_string(index=False))


if __name__ == "__main__":
    main()
