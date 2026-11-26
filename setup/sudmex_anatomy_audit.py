"""Read-only native extraction and provenance audit of all configured SUDMEX T1s.

Raw positive support is descriptive, never equated with brain tissue: scans may
contain scalp, skull, neck, background noise and defacing edits.
"""
import hashlib
import json
from datetime import datetime, timezone
import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import tomllib
except ImportError:
    from pip._vendor import tomli as tomllib
from audit_sudmex_registration import ROOT, OUT, RAW, DERIV, configured_cohort, one, on_grid


def geometry(image):
    return {"shape": list(image.shape), "zooms": [float(v) for v in image.header.get_zooms()],
            "affine": image.affine.tolist(), "orientation": list(nib.aff2axcodes(image.affine)),
            "qform_code": int(image.header["qform_code"]), "sform_code": int(image.header["sform_code"])}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plot_dir = OUT / "native_contact_sheets"
    plot_dir.mkdir(exist_ok=True)
    subjects, sessions = configured_cohort()
    rows, provenance, geometries, plots, config_rows = [], [], [], [], []
    for subject in subjects:
        prefix = f"sub-{subject}_ses-t0"
        raw_path = one(RAW / f"sub-{subject}/ses-t0/anat", "*_T1w.nii.gz")
        anat_dir = DERIV / f"sub-{subject}/ses-t0/anat"
        native_path = anat_dir / f"{prefix}_desc-preproc_T1w.nii.gz"
        mask_path = anat_dir / f"{prefix}_desc-brain_mask.nii.gz"
        raw_img, native_img, mask_img = nib.load(raw_path), nib.load(native_path), nib.load(mask_path)
        native = np.asarray(native_img.dataobj, dtype=np.float32)
        raw_original = np.asarray(raw_img.dataobj, dtype=np.float32)
        raw = on_grid(raw_img, native_img, 1).astype(np.float32)
        mask = on_grid(mask_img, native_img) > 0
        if nib.aff2axcodes(native_img.affine) != ("R", "A", "S"):
            raise ValueError(f"Unexpected orientation for display: {subject}")
        voxvol = abs(np.linalg.det(native_img.affine[:3, :3]))
        components, n_components = ndimage.label(mask)
        counts = np.bincount(components.ravel())[1:]
        native_support, raw_support = native > 0, raw_original > 0
        coords = np.nonzero(mask)
        bounds = [(int(c.min()), int(c.max())) for c in coords]
        x = int(round(np.median(coords[0])))
        y = int(round(bounds[1][0] + .62 * (bounds[1][1] - bounds[1][0])))
        z = int(round(bounds[2][0] + .70 * (bounds[2][1] - bounds[2][0])))
        sections = [(0, x, "mid-sagittal"), (1, y, "anterior coronal"), (2, z, "upper axial")]
        plots.append({"subject": subject,
                      "planes": [(np.take(raw, idx, axis=ax), np.take(native, idx, axis=ax),
                                  np.take(mask, idx, axis=ax), title, idx) for ax, idx, title in sections],
                      "raw_max": float(np.percentile(raw[raw > 0], 98)),
                      "native_max": float(np.percentile(native[native > 0], 98))})
        r = {"subject": subject, "raw_path": str(raw_path.relative_to(ROOT)),
             "native_path": str(native_path.relative_to(ROOT)), "mask_path": str(mask_path.relative_to(ROOT)),
             "native_mask_volume_ml": float(mask.sum() * voxvol / 1000),
             "raw_nonzero_support_ml_NOT_brain_volume": float(raw_support.sum() * abs(np.linalg.det(raw_img.affine[:3, :3])) / 1000),
             "raw_zero_fraction": float((raw_original == 0).mean()),
             "raw_positive_voxels_outside_native_mask_fraction": float(((raw > 0) & ~mask).sum() / (raw > 0).sum()),
             "native_nonzero_support_ml_NOT_brain_volume": float(native_support.sum() * voxvol / 1000),
             "native_mask_largest_component_fraction": float(counts.max() / mask.sum()),
             "native_mask_component_count": int(n_components),
             "native_mask_empty_intensity_fraction": float((mask & (native <= 0)).sum() / mask.sum()),
             "raw_native_grid_matches_within_0_01mm": bool(raw_img.shape == native_img.shape and np.allclose(raw_img.affine, native_img.affine, atol=.01)),
             "raw_native_inmask_intensity_pearson_r": float(np.corrcoef(raw[mask], native[mask])[0, 1]),
             "raw_nonfinite_voxels": int(np.count_nonzero(~np.isfinite(raw_original))),
             "native_nonfinite_voxels": int(np.count_nonzero(~np.isfinite(native))),
             "sagittal_index": x, "coronal_index": y, "axial_index": z,
             "native_json_skullstripped": json.loads(native_path.with_suffix("").with_suffix(".json").read_text()).get("SkullStripped"),
             "report_present": (DERIV / f"sub-{subject}.html").is_file()}
        rows.append(r)
        geometries.append({"subject": subject, "raw": geometry(raw_img), "native": geometry(native_img), "mask": geometry(mask_img)})
        configs = sorted((DERIV / f"sub-{subject}/log").glob("*/fmriprep.toml"))
        if not configs:
            raise ValueError(f"Missing saved configuration: {subject}")
        for cfg_path in configs:
            cfg = tomllib.loads(cfg_path.read_text())
            w, e = cfg["workflow"], cfg["execution"]
            crashes = sorted(cfg_path.parent.glob("crash*"))
            config_rows.append({"subject": subject, "run_uuid": e["run_uuid"],
                                "is_latest_saved_config": cfg_path == configs[-1],
                                "fmriprep_version": cfg["environment"]["version"],
                                "skull_strip_t1w": w["skull_strip_t1w"],
                                "skull_strip_template": w["skull_strip_template"],
                                "run_reconall": w["run_reconall"],
                                "bold2t1w_dof": w["bold2t1w_dof"],
                                "bold2t1w_init": w["bold2t1w_init"], "use_syn_sdc": w["use_syn_sdc"],
                                "output_spaces": e["output_spaces"],
                                "T1_filter": json.dumps(e.get("bids_filters", {}).get("t1w", {})),
                                "BOLD_filter": json.dumps(e.get("bids_filters", {}).get("bold", {})),
                                "crash_files_in_this_run_directory": len(crashes),
                                "crash_file_paths": "|".join(str(p.relative_to(ROOT)) for p in crashes),
                                "config_path": str(cfg_path.relative_to(ROOT))})
        for path in [raw_path, native_path, mask_path] + configs:
            provenance.append({"path": str(path.relative_to(ROOT)), "size_bytes": path.stat().st_size,
                               "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        print(f"Native audit sub-{subject}: mask {r['native_mask_volume_ml']:.1f} ml", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "native_anatomy_metrics.csv", index=False)
    pd.DataFrame(config_rows).to_csv(OUT / "fmriprep_config_inventory.csv", index=False)
    pd.DataFrame(provenance).to_csv(OUT / "anatomy_input_manifest.csv", index=False)
    (OUT / "native_geometry.json").write_text(json.dumps(geometries, indent=2), encoding="utf-8")
    pages = []
    for start in range(0, len(plots), 4):
        part = plots[start:start + 4]
        fig, axes = plt.subplots(len(part), 6, figsize=(18, 3.1 * len(part)), squeeze=False, layout="constrained")
        for i, item in enumerate(part):
            for j, (rp, np_, mp, label, idx) in enumerate(item["planes"]):
                for k, (plane, vmax, kind) in enumerate([(rp, item["raw_max"], "raw"), (np_, item["native_max"], "preproc")]):
                    ax = axes[i, j * 2 + k]
                    ax.imshow(plane.T, origin="lower", cmap="gray", vmin=0, vmax=vmax)
                    ax.contour(mp.T, levels=[.5], colors=["#ff793a"], linewidths=.8)
                    ax.set_title(f"{item['subject']} {kind} | {label} {idx}", fontsize=9)
                    ax.axis("off")
        fig.suptitle("SUDMEX baseline native anatomy | orange: fMRIPrep brain mask on raw and preprocessed T1", fontsize=13)
        path = plot_dir / f"page_{start // 4 + 1:02d}.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        pages.append({"path": str(path.relative_to(OUT)), "subjects": [p["subject"] for p in part]})
    latest = pd.DataFrame(config_rows).query("is_latest_saved_config")
    summary = {"created_utc": datetime.now(timezone.utc).isoformat(), "subjects": subjects,
               "n_subjects": len(subjects), "config_records": len(config_rows), "latest_config_subjects": len(latest),
               "latest_skull_strip_t1w": latest.skull_strip_t1w.value_counts().to_dict(),
               "latest_fmriprep_version": latest.fmriprep_version.value_counts().to_dict(),
               "latest_run_crash_count": int(latest.crash_files_in_this_run_directory.sum()),
               "all_run_directory_crash_count": int(sum(r["crash_files_in_this_run_directory"] for r in config_rows)),
               "metrics": pd.DataFrame(rows).select_dtypes(include=np.number).describe().to_dict(),
               "native_contact_sheets": pages,
               "interpretation_guardrail": "Raw positive support includes head/neck/noise. Mask/raw volume ratio does not measure brain loss. Review images."}
    (OUT / "native_anatomy_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ["metrics", "native_contact_sheets", "subjects"]}, indent=2))


if __name__ == "__main__":
    main()
