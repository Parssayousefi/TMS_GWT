"""Summarize completed numerical audit and the recorded 25 Nov 2026 visual review.

The visual observations below apply to the reviewed contact sheets, not to any
future regeneration. Existing review hashes must match before rewriting report.
"""
import hashlib
import json
from pathlib import Path
import pandas as pd
from audit_sudmex_registration import ROOT, OUT


def md_table(frame):
    header = "| " + " | ".join(frame.columns) + " |"
    rule = "| " + " | ".join(["---"] * len(frame.columns)) + " |"
    return "\n".join([header, rule] + ["| " + " | ".join(str(v) for v in row) + " |"
                                        for row in frame.itertuples(index=False, name=None)])


def main():
    n = pd.read_csv(OUT / "native_anatomy_metrics.csv", dtype={"subject": str})
    a = pd.read_csv(OUT / "standard_anatomy_metrics.csv", dtype={"subject": str})
    r = pd.read_csv(OUT / "run_qc_with_coverage_and_motion.csv", dtype={"subject": str})
    ns = json.loads((OUT / "native_anatomy_summary.json").read_text())
    rs = json.loads((OUT / "alignment_summary.json").read_text())
    figures = ns["native_contact_sheets"] + rs["standard_contact_sheets"]
    reviewed = []
    for figure in figures:
        path = OUT / figure["path"]
        reviewed.append({**figure, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "reviewed": True})
    review_path = OUT / "visual_review_record.json"
    if review_path.exists():
        previous = json.loads(review_path.read_text())
        if previous["figures"] != reviewed:
            raise RuntimeError("Contact sheet content changed: repeat visual review before updating its record")
    review = {"review_date": "2026-11-25", "reviewer": "Codex visual image review",
              "native_scope": "All 44 subjects; raw and native preprocessed T1 at mid-sagittal, anterior-coronal and upper-axial cuts (11 pages).",
              "standard_scope": "All 142 runs; BOLD and T1 at MNI x=-32 and z=36 with both masks (24 pages).",
              "raw_full_head_tissue_visible_subjects": n.subject.tolist(),
              "gross_Aachen_like_native_cortical_overstripping_seen": [],
              "native_observation": "Scalp/skull and neck tissue visible in every raw T1; orange masks follow brain envelope and retain peripheral cortical tissue in the sampled views.",
              "standard_observation": "No gross orientation/normalization failure seen in sampled cuts. Inferior BOLD coverage restriction recurs; focal mask/image boundary differences remain.",
              "figures": reviewed,
              "limits": "Qualitative, slice-sampled review; not a manual whole-volume segmentation or certification of local registration accuracy."}
    review_path.write_text(json.dumps(review, indent=2), encoding="utf-8")
    visual_subjects = n[["subject"]].copy()
    visual_subjects["raw_head_tissue_visible"] = True
    visual_subjects["gross_cortical_overstripping_seen_in_reviewed_planes"] = False
    visual_subjects["native_contact_sheet"] = [f"native_contact_sheets/page_{i // 4 + 1:02d}.png" for i in range(len(n))]
    visual_subjects.to_csv(OUT / "native_visual_review.csv", index=False)
    notes = {
        ("025", "t0"): "Pronounced inferior restriction; already excluded for motion.",
        ("025", "t1"): "Pronounced inferior restriction despite motion INCLUDE; native T1 extraction retained cortex.",
        ("025", "t2"): "Pronounced inferior restriction despite motion INCLUDE; native T1 extraction retained cortex.",
        ("025", "t3"): "Lowest whole-brain mask Dice; pronounced inferior truncation; already excluded for motion.",
        ("033", "t1"): "Pronounced inferior restriction; superior brain alignment retained in sampled views.",
        ("034", "t0"): "Pronounced inferior restriction; superior brain alignment retained in sampled views.",
        ("024", "t3"): "Inferior restriction with lower whole-brain Dice; no gross native extraction failure.",
        ("024", "t1"): "Left HIP-head-m1 parcel has 36.3% of anatomical support in BOLD mask; individual regional coverage caveat.",
        ("001", "t1"): "Relatively conservative cortical mask; many low-coverage parcels are orbitofrontal/temporal-pole; positive BOLD reference extends beyond mask.",
        ("001", "t4"): "Relatively conservative cortical mask and local frontal/posterior differences; positive reference extends beyond mask.",
    }
    for session in ["t0", "t1", "t14", "t2", "t3"]:
        notes[("036", session)] = "Frontal BOLD/mask extends beyond T1 boundary at x=-32; inferior restriction. Local mismatch, not native cortical extraction loss."
    records = []
    for (subject, session), note in notes.items():
        item = r[(r.subject == subject) & (r.session == session)].iloc[0]
        records.append({"subject": subject, "session": session, "existing_motion_status": item.existing_motion_status,
                        "dice": item.dice, "anatomical_cortical_atlas_in_bold_mask_fraction": item.anatomical_cortical_atlas_in_bold_mask_fraction,
                        "note": note, "action": "Preserve existing motion status; consider regional coverage for the planned analysis."})
    pd.DataFrame(records).sort_values(["subject", "session"]).to_csv(OUT / "spatial_review_notes.csv", index=False)
    worst = r.nsmallest(8, "dice")[["subject", "session", "existing_motion_status", "dice", "anatomical_mask_below_z_minus30_in_bold_mask_fraction"]].copy()
    worst.columns = ["Subject", "Session", "Existing motion QC", "BOLD/T1 Dice", "T1 support below z=-30 in BOLD mask"]
    worst["BOLD/T1 Dice"] = worst["BOLD/T1 Dice"].map(lambda x: f"{x:.3f}")
    worst[worst.columns[-1]] = worst.iloc[:, -1].map(lambda x: f"{100*x:.1f}%")
    report = f"""# SUDMEX full-cohort anatomical and registration QC

Audit date: 25 November 2026. Scope: all **44 configured processed participants and 142 saved resting-state runs**, including scans excluded by the existing motion QC. This is a read-only audit of existing preprocessing, not a rerun or a treatment-effect analysis.

**Conclusion: the reviewed SUDMEX inputs and outputs do not show the Aachen failure mechanism. There is no evidence here for repeating SUDMEX T1 extraction solely because `skull_strip_t1w=force` was used.** Every SUDMEX raw T1 visibly retains head tissue, and every reviewed native mask retains the outer cortical envelope. The functional images do have material inferior coverage limitations and some local mismatch; this report does not label every run spatially clean.

## Inventory and effective processing

- All 44 baseline raw T1s, native preprocessed T1s, native masks, standard T1s/masks and participant reports exist.
- All 142 expected runs have a standard BOLD, BOLD reference, brain mask and confounds. Raw BOLD, saved BOLD and confounds each have 300 volumes/rows. Counts: t0=44, t1=44, t14=14, t2=20, t3=14, t4=6. Duplicate t5 sessions remain outside the configured scope.
- The latest saved configuration for every participant is fMRIPrep 23.2.1, `skull_strip_t1w=force`, OASIS30ANTs, no FreeSurfer reconstruction, and baseline-only (`ses-t0`) anatomy. All 142 functional summary reports describe FSL FLIRT boundary-based registration with 6 degrees of freedom and **no susceptibility distortion correction**.
- The inventory contains 45 saved configuration records. The sole retained crash belongs to the superseded 5 November sub-052 run. All 44 latest run directories have zero crash files; current sub-052 t0/t1 outputs and all six rerun sub-002 outputs are present. Historical failures are not counted as current failures.
- Existing motion decisions remain 124 INCLUDE, 9 TRUNCATE, 9 EXCLUDE. This audit did not change them.

## Native extraction: all 44 participants reviewed

All 11 native contact sheets were opened and inspected. Each participant has raw and preprocessed T1 images with the saved brain-mask outline at three cuts: mid-sagittal, anterior-coronal and upper-axial. Every raw T1 shows scalp/skull and neck tissue; several have facial tissue removed, which is not whole-head brain extraction. The raw inputs are therefore unlike the already brain-stripped Aachen inputs.

No broad loss of peripheral cortex, cerebellum or brainstem from the native mask was visible in these cuts. Raw and native images have matching geometry to the checked tolerance. All native masks have a single connected component, and no raw/native T1 nonfinite voxels were found. Native mask volumes range from {n.native_mask_volume_ml.min():.0f} to {n.native_mask_volume_ml.max():.0f} mL (median {n.native_mask_volume_ml.median():.0f} mL); this is supporting descriptive evidence, not a normative pass threshold.

The raw positive-signal volume includes nonbrain head/neck tissue and background signal. A mask/raw-positive-support ratio must **not** be interpreted as the fraction of brain removed. The full-head native preprocessed image itself is also correctly marked `SkullStripped: false`; the separate mask specifies the brain extraction.

## Standard-space correspondence and regional limits

All 24 standard contact sheets were opened and inspected: all 142 runs at MNI x=-32 and z=36, side-by-side with T1, with BOLD and anatomical masks overlaid. No gross orientation or normalization failure was apparent in the reviewed cuts. Independent comparison of anatomical masks to the matching MNI template gives Dice {a.MNI_anatomy_template_mask_dice.min():.3f}–{a.MNI_anatomy_template_mask_dice.max():.3f} (median {a.MNI_anatomy_template_mask_dice.median():.3f}); template-brain containment in anatomical masks is {100*a.MNI_template_brain_in_anatomy_mask_fraction.min():.1f}–{100*a.MNI_template_brain_in_anatomy_mask_fraction.max():.1f}%.

Across all 142 runs, BOLD/T1 mask Dice is {r.dice.min():.3f}–{r.dice.max():.3f} (median {r.dice.median():.3f}); BOLD/T1 mask volume ratio is {r.bold_to_T1_mask_volume_ratio.min():.3f}–{r.bold_to_T1_mask_volume_ratio.max():.3f}. These whole-brain summaries conceal local limitations:

- Above or at MNI z=0, median anatomical-mask support retained by the BOLD mask is {100*r.anatomical_mask_at_or_above_z0_in_bold_mask_fraction.median():.1f}% (range {100*r.anatomical_mask_at_or_above_z0_in_bold_mask_fraction.min():.1f}–{100*r.anatomical_mask_at_or_above_z0_in_bold_mask_fraction.max():.1f}%). Below z=-30 it is {100*r.anatomical_mask_below_z_minus30_in_bold_mask_fraction.median():.1f}% (range {100*r.anatomical_mask_below_z_minus30_in_bold_mask_fraction.min():.1f}–{100*r.anatomical_mask_below_z_minus30_in_bold_mask_fraction.max():.1f}%). These are geometric slabs, not tissue labels.
- Of cortical atlas support lying inside each participant's T1 mask, {100*r.anatomical_cortical_atlas_in_bold_mask_fraction.min():.1f}–{100*r.anatomical_cortical_atlas_in_bold_mask_fraction.max():.1f}% lies in the BOLD mask (median {100*r.anatomical_cortical_atlas_in_bold_mask_fraction.median():.1f}%). Corresponding subcortical support is {100*r.anatomical_subcortical_atlas_in_bold_mask_fraction.min():.1f}–{100*r.anatomical_subcortical_atlas_in_bold_mask_fraction.max():.1f}%.
- Per run, 2–29 cortical parcels have less than half of their anatomical support inside the BOLD mask. Many are orbitofrontal or temporal-pole parcels. This threshold is a descriptive flag, not a validated exclusion rule. Sub-024 t1 also has 36.3% coverage of the left HIP-head-m1 parcel. Mask coverage does not measure temporal signal quality.
- Sub-025 has pronounced inferior restriction across all saved sessions; its t3 run has the lowest Dice and was already excluded for motion. Sub-033 t1 and sub-034 t0 also show pronounced inferior restriction. Sub-036 shows frontal BOLD/mask extension beyond the T1 boundary in the reviewed sagittal cut across sessions. Sub-001 t1/t4 have relatively conservative cortical masks despite much broader positive BOLD-reference support.

{md_table(worst)}

The inferior mismatch is visible in the functional references and masks while the native T1 cortical boundary is preserved. This supports a functional coverage/local alignment problem rather than the Aachen T1 over-extraction mechanism. Its precise contribution from acquisition field of view, susceptibility distortion, registration and mask estimation was not separated by an intervention. A blanket T1 rerun is not justified by this evidence; inferior/regional analyses still need explicit coverage handling, and the noted focal cases should not be described as unqualified passes.

## Deliverables and limits

- `native_anatomy_metrics.csv`, `fmriprep_config_inventory.csv`, `anatomy_input_manifest.csv`, `native_geometry.json`: native measures, effective settings, input hashes and geometry.
- `alignment_metrics.csv`, `standard_anatomy_metrics.csv`, `run_qc_with_coverage_and_motion.csv`: all-run and all-subject quantitative audit.
- `parcel_mask_coverage.csv`, `poorly_covered_parcels_with_labels.csv`, `spatial_review_notes.csv`: local coverage and specific reviewed cases, separate from existing motion exclusions.
- `native_contact_sheets/` (11 pages), `standard_contact_sheets/` (24 pages), `visual_review_record.json`: reviewed images, exact scope and image hashes. The six older single-subject PNGs are from the earlier acute-only audit; the contact sheets are the full-cohort review.

Visual review sampled three native and two standard planes; it was not manual whole-volume segmentation. No ground-truth extraction mask, quantitative distortion map, or pilot rerun was used. The volume-count checks and finite BOLD-reference check do not constitute full inspection of every 4-D time series, denoising, motion-salvage implementation, or downstream parcellated data. This audit addresses the requested preprocessing anatomy/registration question without claiming that those other steps have been validated.

Reproduce numerical outputs using `C:/Miniforge3/envs/gwt/python.exe -X utf8` with `setup/sudmex_anatomy_audit.py`, `setup/audit_sudmex_registration.py`, and `setup/sudmex_anatomy_coverage.py`. `setup/sudmex_anatomy_report.py` summarizes them and refuses to reuse the visual review if the recorded contact-sheet hashes change. All writes are confined to this audit output and its setup scripts.
"""
    (OUT / "SUDMEX_FULL_PREPROCESSING_QC.md").write_text(report, encoding="utf-8")
    print(OUT / "SUDMEX_FULL_PREPROCESSING_QC.md")
    print(f"Reviewed {len(ns['native_contact_sheets'])} native and {len(rs['standard_contact_sheets'])} standard sheets")


if __name__ == "__main__":
    main()
