"""Independent anatomical QA of completed corrected Aachen participants.

Only reads production data and queue state. All writes stay in the dedicated
review directory. No processing, statistical effects, or subject exclusion is
performed. Automatic screens require explicit visual pass/flag review.
"""
from argparse import ArgumentParser
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

import nibabel as nib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.dont_write_bytecode = True
import aachen_alignment_anatomical_review as audit
import aachen_alignment_map_targets as mapper

ROOT = audit.ROOT
COHORT = audit.REVIEW / "corrected_cohort"
QUEUE = ROOT / "logs/aachen_corrected_20261125/queue_status.json"
DERIVATIVES = ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125"
SUBJECTS = [f"{i:02d}" for i in range(1, 17) if i != 2]
SESSIONS = ["pre90", "post90", "pre120", "post120"]
TABLES = ["anatomical_support", "functional_anatomical_agreement", "target_support", "effective_skullstrip_settings"]


def stamp():
    return datetime.now(timezone.utc).isoformat()


def read_table(path):
    try:
        return pd.read_csv(path, dtype={"subject": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def evidence(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size}


def queue_ready(state, subject):
    item = state.get("finished", {}).get(subject, {})
    return (item.get("success") is True and item.get("state", {}).get("ExitCode") == 0
            and item.get("inventory", {}).get("complete") is True)


def get_template():
    folder = ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym"
    t1 = nib.load(folder / "tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz")
    mask = np.asarray(nib.load(folder / "tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz").dataobj) > 0
    return t1, np.asarray(t1.dataobj) * mask, mask


def show(ax, data, index, axis, title, overlays=()):
    positive = data[data > 0]
    ax.imshow(np.take(data, index, axis=axis).T, origin="lower", cmap="gray", vmin=0,
              vmax=np.percentile(positive, 98) if len(positive) else 1)
    for mask, color in overlays:
        plane = np.take(mask, index, axis=axis).T
        if plane.min() < .5 < plane.max():
            ax.contour(plane, levels=[.5], colors=[color], linewidths=.9)
    ax.set_title(title, fontsize=9)
    ax.axis("off")


def make_visuals(subject, out):
    prefix = f"sub-{subject}_ses-anat"
    anat = DERIVATIVES / f"sub-{subject}/ses-anat/anat"
    native_path = anat / f"{prefix}_desc-preproc_T1w.nii.gz"
    native = audit.sitk.ReadImage(str(native_path))
    native_data = audit.array(native)
    native_mask = np.asarray(nib.load(anat / f"{prefix}_desc-brain_mask.nii.gz").dataobj) > 0
    center = np.rint(np.argwhere(native_data > 0).mean(0)).astype(int)
    atlas_path = ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz"
    atlas = audit.sitk.ReadImage(str(atlas_path))
    labels = audit.array(atlas)
    cortex = audit.mask_image((labels >= 1) & (labels <= 400), atlas)
    inverse = audit.sitk.ReadTransform(str(anat / f"{prefix}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"))
    native_cortex = audit.array(audit.sitk.Resample(cortex, native, inverse,
                                                  audit.sitk.sitkNearestNeighbor, 0, audit.sitk.sitkUInt8)) > 0
    lateral = center.copy()
    lateral[0] = int(np.percentile(np.argwhere(native_data > 0)[:, 0], 25))
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), layout="constrained")
    for axis in range(3):
        show(axes[0, axis], native_data, center[axis], axis, f"Native T1 mask, axis {axis}",
             [(native_mask, "#ff6339")])
        show(axes[1, axis], native_data, lateral[axis], axis, f"Cortical atlas on native T1, axis {axis}",
             [(native_cortex, "#39c5dc")])
    fig.suptitle(f"sub-{subject}: corrected native anatomy; cortex/cerebellum and atlas placement")
    fig.savefig(out / "native_anatomy_and_atlas.png", dpi=130)
    # Replace the inherited audit's original-failure caption in this new output
    # directory with a neutral corrected-anatomy view; source audit is untouched.
    fig.savefig(out / "all_subjects_native_t1_brainmask.png", dpi=130)
    plt.close(fig)

    template_nib, template, template_mask = get_template()
    standard_nib = nib.load(anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz")
    np.testing.assert_allclose(standard_nib.affine, template_nib.affine, atol=1e-4)
    standard = np.asarray(standard_nib.dataobj)
    standard_mask = np.asarray(nib.load(anat / f"{prefix}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz").dataobj) > 0
    indices = np.rint(nib.affines.apply_affine(np.linalg.inv(template_nib.affine), [0, 0, 25])).astype(int)
    fig, axes = plt.subplots(6, 3, figsize=(12, 19), layout="constrained")
    for axis in range(3):
        show(axes[0, axis], template, indices[axis], axis, "MNI template")
        show(axes[1, axis], standard, indices[axis], axis, "Corrected MNI T1 + template boundary",
             [(template_mask, "#39c5dc")])
    for row, session in enumerate(SESSIONS, 2):
        path = DERIVATIVES / f"sub-{subject}/ses-{session}/func/sub-{subject}_ses-{session}_task-rest_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz"
        image = nib.load(path)
        np.testing.assert_allclose(image.affine, template_nib.affine, atol=1e-4)
        data = np.asarray(image.dataobj)
        for axis in range(3):
            show(axes[row, axis], data, indices[axis], axis, f"{session} BOLD + corrected T1 boundary",
                 [(standard_mask, "#ff6339")])
    fig.suptitle(f"sub-{subject}: same standard-space slices in all four sessions; x=0, y=0, z=25 mm")
    fig.savefig(out / "standard_anatomy_all_sessions.png", dpi=120)
    plt.close(fig)
    return {"template_t1_intensity_r": audit.correlation(standard, template, template_mask),
            "template_brainmask_dice": audit.dice(standard_mask, template_mask)}


def audit_subject(subject, queue):
    if not queue_ready(queue, subject):
        raise ValueError(f"Queue success/output completion not established for {subject}")
    out = COHORT / f"sub-{subject}"
    if (out / "metrics.json").exists():
        return False
    if subject == "06" and not out.exists():
        shutil.copytree(audit.REVIEW / "pilot_skip", out)
    else:
        argv, previous = sys.argv, audit.REVIEW
        try:
            audit.REVIEW = COHORT
            sys.argv = ["aachen_alignment_anatomical_review.py", "--derivatives-root", str(DERIVATIVES),
                        "--subjects", subject, "--label", f"sub-{subject}", "--save-targets"]
            audit.main()
        finally:
            sys.argv, audit.REVIEW = argv, previous
    a = read_table(out / "anatomical_support.csv")
    f = read_table(out / "functional_anatomical_agreement.csv")
    t = read_table(out / "target_support.csv")
    c = read_table(out / "effective_skullstrip_settings.csv")
    flags = []
    screens = {
        "native_support_preserved_over_95pct": bool(a.fraction_native_support_retained.min() > .95),
        "outer_support_preserved_over_90pct": bool(a.fraction_outer_8mm_support_retained.min() > .90),
        "skull_strip_skip": set(c.skull_strip_t1w) == {"skip"},
        "all_four_sessions": set(f.session) == set(SESSIONS),
        "bold_t1_dice_all_over_0_85": bool(f.bold_to_masked_t1_mni_dice.min() > .85),
        "bold_transform_r_all_over_0_99": bool(f.bold_transform_reproduction_r.min() > .99),
        "mni_affines_match": bool((f.bold_t1_mni_affine_max_abs_diff == 0).all()),
        "t1_mask_transform_dice_over_0_98": bool(a.native_to_saved_mni_mask_dice.min() > .98),
    }
    shapes = []
    for session in SESSIONS:
        stem = f"sub-{subject}_ses-{session}_task-rest"
        raw = nib.load(ROOT / f"data/aachen_itbs/sub-{subject}/ses-{session}/func/{stem}_bold.nii.gz")
        processed = nib.load(DERIVATIVES / f"sub-{subject}/ses-{session}/func/{stem}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz")
        shapes.append({"session": session, "raw_volumes": int(raw.shape[3]), "processed_volumes": int(processed.shape[3]),
                       "same_volume_count": raw.shape[3] == processed.shape[3]})
    screens["bold_volume_counts_preserved"] = all(r["same_volume_count"] for r in shapes)
    flags.extend(name for name, passed in screens.items() if not passed)
    image_metrics = make_visuals(subject, out)
    if image_metrics["template_t1_intensity_r"] < .65:
        flags.append("low_template_tissue_intensity_agreement_requires_review")
    if not t.empty:
        if t.fraction_target_in_native_fmriprep_t1_mask.min() < .8:
            flags.append("source_roi_native_t1_coverage_below_80pct_requires_review")
        if t.fraction_target_in_existing_cortical_atlas.min() <= .1:
            flags.append("source_roi_cortical_atlas_coverage_at_or_below_original_10pct_rule")
    metrics = {"subject": subject, "created_utc": stamp(), "derivative_root": str(DERIVATIVES),
               "queue_success": queue["finished"][subject], "automatic_screens": screens,
               "automatic_flags": flags, "native_support_retained": float(a.fraction_native_support_retained.iloc[0]),
               "outer_support_retained": float(a.fraction_outer_8mm_support_retained.iloc[0]),
               "bold_t1_mask_dice_min": float(f.bold_to_masked_t1_mni_dice.min()),
               "bold_t1_mask_dice_max": float(f.bold_to_masked_t1_mni_dice.max()),
               "bold_transform_r_min": float(f.bold_transform_reproduction_r.min()),
               "target_rows": len(t), "volume_counts": shapes, **image_metrics,
               "review_status": "pending_visual_review", "functional_outcomes_read": False,
               "screen_note": "Screens flag problems for review; they never exclude participants or move ROIs."}
    write_json(out / "metrics.json", metrics)
    write_json(out / "queue_completion_evidence.json", {"captured_utc": stamp(), "queue_updated_utc": queue.get("updated_utc"),
                                                       "subject": subject, "completion": queue["finished"][subject]})
    print(json.dumps({"audited": subject, "flags": flags, "native_support": metrics["native_support_retained"],
                      "minimum_bold_t1_dice": metrics["bold_t1_mask_dice_min"]}), flush=True)
    return True


def summary():
    rows = []
    for subject in SUBJECTS:
        folder = COHORT / f"sub-{subject}"
        if not (folder / "metrics.json").exists():
            rows.append({"subject": subject, "status": "not_audited"})
            continue
        m = json.loads((folder / "metrics.json").read_text())
        review = json.loads((folder / "visual_review.json").read_text()) if (folder / "visual_review.json").exists() else {}
        rows.append({"subject": subject, "status": review.get("decision", "pending_visual_review"),
                     "automatic_flags": ";".join(m["automatic_flags"]), "native_support_retained": m["native_support_retained"],
                     "outer_support_retained": m["outer_support_retained"], "bold_t1_dice_min": m["bold_t1_mask_dice_min"],
                     "template_t1_r": m["template_t1_intensity_r"], "target_rows": m["target_rows"]})
    pd.DataFrame(rows).to_csv(COHORT / "subject_status.csv", index=False)
    return rows


def record_review(subject, decision, note):
    out = COHORT / f"sub-{subject}"
    metrics = json.loads((out / "metrics.json").read_text())
    if not note:
        raise ValueError("A concrete visual finding is required")
    if decision == "pass" and metrics["automatic_flags"]:
        raise ValueError("Automatic flags remain: use flag decision, investigate and document before release")
    reference_folder = out / "reference_intensity"
    if decision == "pass" and reference_folder.exists():
        resolution_path = reference_folder / "resolution_receipt.json"
        if not resolution_path.exists() or json.loads(resolution_path.read_text()).get("passed") is not True:
            raise ValueError("Reference-intensity investigation must be resolved before visual pass")
    files = [out / "native_anatomy_and_atlas.png", out / "standard_anatomy_all_sessions.png", out / "metrics.json"]
    files += sorted(out.glob("*target_chain.png"))
    bilateral = sorted(out.glob("*bilateral_target_planes.png"))
    if metrics["target_rows"] and not bilateral:
        raise ValueError("Bilateral three-plane target views must be generated and reviewed first")
    files += bilateral
    coverage_path = out / "atlas_functional_coverage_summary.json"
    if not coverage_path.exists():
        raise ValueError("Atlas functional coverage must be checked before visual review is recorded")
    files.append(coverage_path)
    files += sorted((out / "reference_intensity").glob("*.png"))
    if (out / "reference_intensity/resolution_receipt.json").exists():
        files += [out / "reference_intensity/resolution_receipt.json", out / "reference_intensity/RESOLUTION.md"]
    write_json(out / "visual_review.json", {"subject": subject, "reviewed_utc": stamp(), "decision": decision,
                                            "note": note, "evidence": [evidence(p) for p in files],
                                            "scope": "Anatomical alignment; not coil/peak-stimulation verification or functional effects."})


def combine_audits():
    out = COHORT / "full_audit"
    out.mkdir(exist_ok=True)
    hashes, headers = {}, []
    for table in TABLES:
        frames = [read_table(COHORT / f"sub-{s}/{table}.csv") for s in SUBJECTS]
        pd.concat([f for f in frames if not f.empty], ignore_index=True).to_csv(out / f"{table}.csv", index=False)
    for subject in SUBJECTS:
        folder = COHORT / f"sub-{subject}"
        p = json.loads((folder / "provenance.json").read_text())
        for path, digest in p["input_sha256"].items():
            if path in hashes and hashes[path] != digest:
                raise ValueError(f"Inconsistent audited input: {path}")
            hashes[path] = digest
        headers.extend(json.loads((folder / "headers.json").read_text()))
    write_json(out / "headers.json", headers)
    write_json(out / "provenance.json", {"created_utc": stamp(), "derivatives_root": str(DERIVATIVES),
                                         "subjects": SUBJECTS, "input_sha256": hashes,
                                         "individual_audits": [str(COHORT / f"sub-{s}") for s in SUBJECTS],
                                         "functional_outcomes_read": False, "production_inputs_modified": False,
                                         "support_definition": "Nonzero native T1 image support, not independent tissue segmentation.",
                                         "purpose": "Combined completed-subject audit with explicit visual QA."})
    return out


def release_inputs():
    """Hash the exact 3D anatomical evidence and transforms, never 4D signals."""
    paths = {ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz"}
    template = ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym"
    paths.update(template / name for name in ["tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz",
                                            "tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz"])
    paths.update(r["path"] for r in audit.selected_masks())
    for subject in SUBJECTS:
        base = DERIVATIVES / f"sub-{subject}"
        stem = f"sub-{subject}_ses-anat"
        paths.add(ROOT / f"data/aachen_itbs/sub-{subject}/ses-anat/anat/{stem}_T1w.nii.gz")
        paths.update((base / "ses-anat/anat").glob("*.h5"))
        for suffix in ["_desc-preproc_T1w.nii.gz", "_desc-brain_mask.nii.gz",
                       "_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz",
                       "_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz"]:
            paths.add(base / f"ses-anat/anat/{stem}{suffix}")
        paths.update((base / "log").glob("*/fmriprep.toml"))
        for session in SESSIONS:
            stem = f"sub-{subject}_ses-{session}_task-rest"
            func = base / f"ses-{session}/func"
            for suffix in ["_desc-coreg_boldref.nii.gz", "_from-boldref_to-T1w_mode-image_desc-coreg_xfm.txt",
                           "_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz",
                           "_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz"]:
                paths.add(func / f"{stem}{suffix}")
    result = [evidence(p) for p in sorted(paths)]
    write_json(COHORT / "release_input_manifest.json", {"created_utc": stamp(), "files": result,
        "scope": "Main cohort audit: 3D anatomy, references, masks, transforms, source ROIs and configurations. Bounded reference-intensity investigations sampled volumes 1/150 of raw and processed pre90/post90 BOLD; their input hashes are in the corresponding sub-XX/reference_intensity folders."})


def cohort_findings(combined):
    corrected = read_table(combined / "anatomical_support.csv")
    original = read_table(audit.REVIEW / "original/anatomical_support.csv")
    original = original[original.subject.isin(SUBJECTS)]
    paired = original.merge(corrected, on="subject", suffixes=("_original", "_corrected"), validate="one_to_one")
    fields = ["fraction_native_support_retained", "fraction_outer_8mm_support_retained", "fmriprep_native_mask_mL"]
    paired[["subject"]+[f"{f}_{version}" for f in fields for version in ["original", "corrected"]]].to_csv(
        COHORT / "anatomy_original_corrected_comparison.csv", index=False)
    original_target = pd.read_csv(ROOT / "output/aachen_target_anatomy/parcel_coverage.csv")
    corrected_target = pd.read_csv(COHORT / "target_mapping/parcel_coverage.csv")
    target = original_target.merge(corrected_target, on=["subject_id", "side"],
                                   suffixes=("_original", "_corrected"), validate="one_to_one")
    if len(target) != 16:
        raise ValueError("Expected eight participants with unchanged bilateral source definitions")
    target.to_csv(COHORT / "target_coverage_original_corrected_comparison.csv", index=False)
    function = read_table(combined / "functional_anatomical_agreement.csv")
    coverage = read_table(COHORT / "atlas_functional_coverage_summary.csv")
    if set(coverage.subject) != set(SUBJECTS) or not (coverage.sessions == 4).all():
        raise ValueError("Full-cohort atlas functional coverage evidence is incomplete")
    weights = pd.read_csv(COHORT / "target_mapping/subject_region_weights.csv")
    groups = weights.groupby(["subject_id", "side"])
    expected_groups = int(corrected_target.eligible.sum())
    if len(groups) != expected_groups or not (groups.size() == 454).all():
        raise ValueError("Target weight rows are incomplete")
    np.testing.assert_allclose(groups.weight.sum(), 1, atol=1e-10)
    if not np.isfinite(weights.weight).all() or (weights.weight < 0).any() or (weights.loc[weights.roi > 400, "weight"] != 0).any():
        raise ValueError("Target weights violate the cortical-only definition")
    quantitative = {
        "native_support_retained_original_range": [float(original.fraction_native_support_retained.min()), float(original.fraction_native_support_retained.max())],
        "native_support_retained_corrected_range": [float(corrected.fraction_native_support_retained.min()), float(corrected.fraction_native_support_retained.max())],
        "outer_8mm_support_retained_corrected_range": [float(corrected.fraction_outer_8mm_support_retained.min()), float(corrected.fraction_outer_8mm_support_retained.max())],
        "bold_t1_mask_dice_corrected_range": [float(function.bold_to_masked_t1_mni_dice.min()), float(function.bold_to_masked_t1_mni_dice.max())],
        "bold_transform_reproduction_r_minimum": float(function.bold_transform_reproduction_r.min()),
        "target_cortical_fraction_corrected_range": [float(corrected_target.cortical_fraction.min()), float(corrected_target.cortical_fraction.max())],
        "eligible_target_sides_original": int(original_target.eligible.sum()),
        "eligible_target_sides_corrected": int(corrected_target.eligible.sum()),
        "zero_reference_parcel_sessions": int(coverage.zero_reference_parcel_sessions.sum()),
        "zero_brainmask_parcel_sessions": int(coverage.zero_brainmask_all_parcel_sessions.sum()),
        "minimum_parcel_reference_support_fraction": float(coverage.minimum_nonzero_reference_coverage.min()),
        "minimum_parcel_brainmask_support_fraction": float(coverage.minimum_brainmask_coverage.min()),
    }
    write_json(COHORT / "cohort_quantitative_summary.json", quantitative)
    text = ["# Corrected Aachen anatomical release", "",
        "All 15 included participants passed native and standard-space anatomical review after successful corrected processing. All four sessions were reviewed per participant (60 runs). The separate sub-02 slice-timing exclusion is unchanged.", "",
        "The original forced skull stripping truncated already stripped T1 inputs. Skipping this second extraction preserves the supplied cortical and cerebellar support and yields coherent native-atlas, template/T1 and BOLD/T1 registration. The target voxel masks were never relocated or replaced.", "",
        f"Native T1 support retained: original {quantitative['native_support_retained_original_range'][0]:.1%}–{quantitative['native_support_retained_original_range'][1]:.1%}; corrected {quantitative['native_support_retained_corrected_range'][0]:.2%}–{quantitative['native_support_retained_corrected_range'][1]:.2%}. Corrected outer 8-mm support retention: {quantitative['outer_8mm_support_retained_corrected_range'][0]:.2%}–{quantitative['outer_8mm_support_retained_corrected_range'][1]:.2%}. These are image-support measurements, not independent gray-matter segmentation.", "",
        f"Across 60 runs, corrected BOLD/T1 mask Dice is {quantitative['bold_t1_mask_dice_corrected_range'][0]:.3f}–{quantitative['bold_t1_mask_dice_corrected_range'][1]:.3f}; transform reconstruction correlation is at least {quantitative['bold_transform_reproduction_r_minimum']:.4f}.", "",
        f"Eight participants have deposited usable bilateral analysis ROIs. Corrected cortical overlap ranges from {quantitative['target_cortical_fraction_corrected_range'][0]:.1%} to {quantitative['target_cortical_fraction_corrected_range'][1]:.1%}. The unchanged >10% rule retains {quantitative['eligible_target_sides_corrected']} of 16 sides (originally {quantitative['eligible_target_sides_original']}). Weights sum to one per eligible subject/side and have nonzero entries only for cortical atlas parcels.", "",
        f"There are {quantitative['zero_reference_parcel_sessions']} parcel/session combinations with zero reference support and {quantitative['zero_brainmask_parcel_sessions']} with zero brainmask coverage. Partial support remains: the minimum parcel reference and brainmask coverage fractions are {quantitative['minimum_parcel_reference_support_fraction']:.1%} and {quantitative['minimum_parcel_brainmask_support_fraction']:.1%}. This is coverage evidence, not a temporal signal-quality assessment.", "",
        "Pre90 reference-contrast differences in sub-10 and sub-14 were investigated before release. Raw/HMC input references were unchanged, while enhanced native references differed. Saved workflow metadata and four actual processed BOLD samples per investigated participant distinguish the N4 reference enhancement from the 4D intensity input; sampled-volume/T1 alignment remained coherent. The detailed resolutions, input hashes and limitations are in the corresponding `sub-XX/reference_intensity/` folders. Spatial intensity nonuniformity remains and no full temporal assessment is implied.", "",
        "Limits: supplied inferior tissue remnants, EPI susceptibility distortion, incomplete inferior field of view, and mean-functional/reference motion remain. Source descriptions identify intended stimulation-site analysis ROIs, but exact coil positions and peak electric fields were not independently available. Deposited native ROI voxels are unchanged; their MNI centroids follow the corrected transforms, without shifts toward nominal DLPFC coordinates. The right ROI is contralateral, not an independently stimulated target. These corrected weights require downstream time series generated through the same corrected anatomical chain. No functional effects were examined.", "",
        "The release receipt links hashed inputs, per-subject quantitative evidence and visual decisions, combined audits and mapping outputs. `subject_status.csv` records every participant; no participant was discarded based on the anatomical findings or treatment effects.", ""]
    (COHORT / "COHORT_FINDINGS.md").write_text("\n".join(text), encoding="utf-8")
    return quantitative


def finalize(queue):
    rows = summary()
    if not all(queue_ready(queue, s) for s in SUBJECTS):
        raise ValueError("Not all 15 processing jobs completed successfully")
    if not all(r["status"] == "pass" for r in rows):
        raise ValueError("Not all 15 participants have a visual pass")
    reference_checks = []
    for subject in SUBJECTS:
        folder = COHORT / f"sub-{subject}/reference_intensity"
        if not folder.exists():
            continue
        resolution = json.loads((folder / "resolution_receipt.json").read_text())
        if resolution.get("passed") is not True:
            raise ValueError(f"Reference-versus-4D intensity distinction remains unresolved: {subject}")
        for item in resolution["evidence"]:
            if evidence(Path(item["path"]))["sha256"] != item["sha256"]:
                raise ValueError(f"Reference investigation evidence changed: {subject}")
        reference_checks.append({"subject": subject, "sessions": resolution["sampled_sessions"],
            "raw_and_processed_volume_numbers": resolution["raw_and_processed_volume_numbers"],
            "reference_enhancement_distinction_resolved": True})
    if {item["subject"] for item in reference_checks} != {"10", "14"}:
        raise ValueError("Expected documented reference investigations for sub-10 and sub-14")
    for subject in SUBJECTS:
        folder = COHORT / f"sub-{subject}"
        review = json.loads((folder / "visual_review.json").read_text())
        metrics = json.loads((folder / "metrics.json").read_text())
        if metrics["automatic_flags"]:
            raise ValueError(f"Unresolved automatic anatomical flags: {subject}")
        required = set(folder.glob("*bilateral_target_planes.png"))
        if 2*len(required) != metrics["target_rows"]:
            raise ValueError(f"Missing bilateral target views: {subject}")
        if not (folder / "atlas_functional_coverage_summary.json").exists():
            raise ValueError(f"Missing atlas functional coverage: {subject}")
        if not required <= {Path(item["path"]) for item in review["evidence"]}:
            raise ValueError(f"Bilateral three-plane views were not included in visual review: {subject}")
        for item in review["evidence"]:
            if evidence(Path(item["path"]))["sha256"] != item["sha256"]:
                raise ValueError("Reviewed evidence changed")
    combined = combine_audits()
    argv, previous = sys.argv, mapper.REVIEW
    try:
        mapper.REVIEW = COHORT
        sys.argv = ["aachen_alignment_map_targets.py", "--derivatives-root", str(DERIVATIVES),
                    "--audit-dir", str(combined), "--label", "target_mapping"]
        mapper.main()
    finally:
        sys.argv, mapper.REVIEW = argv, previous
    quantitative = cohort_findings(combined)
    release_inputs()
    write_json(COHORT / "queue_final_snapshot.json", queue)
    paths = [COHORT / "subject_status.csv", COHORT / "queue_final_snapshot.json", combined / "provenance.json",
             COHORT / "target_mapping/anatomy_provenance.json", COHORT / "target_mapping/subject_region_weights.csv",
             COHORT / "target_mapping/parcel_coverage.csv", COHORT / "release_input_manifest.json",
             COHORT / "COHORT_FINDINGS.md", COHORT / "COHORT_QA.md", COHORT / "release_guard_check.json", Path(__file__),
             ROOT / "setup/aachen_alignment_atlas_coverage.py", ROOT / "setup/aachen_alignment_target_planes.py",
             ROOT / "setup/aachen_alignment_reference_intensity.py", ROOT / "setup/aachen_alignment_sampled_bold.py",
             ROOT / "setup/aachen_alignment_verify_release.py",
             ROOT / "setup/aachen_alignment_reference_spotcheck.py",
             ROOT / "setup/aachen_alignment_anatomical_review.py", ROOT / "setup/aachen_alignment_map_targets.py"]
    paths += [COHORT / f"sub-{s}/visual_review.json" for s in SUBJECTS]
    paths += sorted(COHORT.glob("*.csv")) + [COHORT / "cohort_quantitative_summary.json"]
    paths += sorted(combined.glob("*.csv"))
    paths += sorted((COHORT / "target_mapping").glob("*.nii.gz"))
    paths += sorted((COHORT / "target_mapping").glob("*.csv"))
    for subject in SUBJECTS:
        folder = COHORT / f"sub-{subject}"
        paths += sorted(folder.glob("*.csv")) + sorted(folder.glob("*.json"))
        paths += sorted((folder / "reference_intensity").glob("*"))
    paths = sorted(set(paths))
    receipt = {"schema_version": 1, "created_utc": stamp(), "passed": True,
               "scope": "Full corrected Aachen anatomical QA and deposited-mask mapping",
               "derivative_root": str(DERIVATIVES), "subjects": SUBJECTS, "required_sessions": SESSIONS,
               "participants_passed": len(rows), "functional_runs_audited": len(SUBJECTS)*len(SESSIONS),
               "mapping_subjects": sorted({r['subject'] for r in audit.selected_masks()}),
               "source_mask_definitions_unchanged": True, "functional_outcomes_read": False,
               "temporal_signal_quality_assessed": False,
               "bounded_spatial_sample_checks": reference_checks,
               "release_for": "Anatomically consistent downstream extraction and unchanged deposited-ROI mapping; functional signal quality and treatment effects require separate assessment.",
               "quantitative_summary": quantitative,
               "evidence": [evidence(p) for p in paths],
               "limitations": ["Image support is not an independent gray-matter segmentation.",
                               "Source ROI identity is supported by study methods/README, not independently verified coil locations or electric fields.",
                               "Deposited native ROI voxels are unchanged; MNI centroids follow the corrected transforms without shifts toward nominal DLPFC coordinates.",
                               "Residual EPI distortion and mean-functional/reference motion remain.",
                               "Inferior field-of-view and uneven parcel coverage remain; reference support is not a temporal signal-quality test.",
                               "Sub-10 and sub-14 pre90 have changed reference-enhancement contrast. Saved dataflow and sampled actual BOLD volumes distinguish reference-only N4 enhancement from the 4D intensity input; full temporal quality remains separate.",
                               "Downstream time series and effects must use the corrected anatomical chain; old effects are not validated.",
                               "Sub-02 remains excluded for the separate documented slice-timing issue."],
               "production_inputs_modified_by_review": False}
    write_json(COHORT / "cohort_release_receipt.json", receipt)
    print("Full cohort release passed", COHORT / "cohort_release_receipt.json", flush=True)


def main():
    parser = ArgumentParser()
    parser.add_argument("--subjects", nargs="*")
    parser.add_argument("--record-review")
    parser.add_argument("--decision", choices=["pass", "flag"])
    parser.add_argument("--note")
    parser.add_argument("--finalize", action="store_true")
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    COHORT.mkdir(parents=True, exist_ok=True)
    if args.watch:
        print("Watching corrected queue for completed participants; 45-second checks, no processing jobs started.", flush=True)
        reported_failures = set()
        while True:
            try:
                queue = json.loads(QUEUE.read_text())
            except (OSError, json.JSONDecodeError):
                time.sleep(45)
                continue
            changed = False
            for subject in SUBJECTS:
                item = queue.get("finished", {}).get(subject, {})
                if item and item.get("success") is False and subject not in reported_failures:
                    print(json.dumps({"queue_failure": subject, "completion": item}), flush=True)
                    reported_failures.add(subject)
                if queue_ready(queue, subject):
                    changed = audit_subject(subject, queue) or changed
            if changed:
                summary()
            if all((COHORT / f"sub-{s}/metrics.json").exists() for s in SUBJECTS):
                print("All 15 quantitative audits and review figures are ready; explicit visual passes are still required.", flush=True)
                return
            time.sleep(45)
    queue = json.loads(QUEUE.read_text())
    if args.record_review:
        record_review(args.record_review.zfill(2), args.decision, args.note)
    elif args.finalize:
        finalize(queue)
    else:
        for subject in args.subjects or SUBJECTS:
            subject = subject.zfill(2)
            if queue_ready(queue, subject):
                audit_subject(subject, queue)
    rows = summary()
    print(json.dumps({"queue_status": queue.get("status"), "queue_running": queue.get("running", []),
                      "qa_passed": [r['subject'] for r in rows if r['status'] == 'pass'],
                      "qa_pending_visual": [r['subject'] for r in rows if r['status'] == 'pending_visual_review'],
                      "qa_flagged": [r['subject'] for r in rows if r['status'] == 'flag']}), flush=True)


if __name__ == "__main__":
    main()
