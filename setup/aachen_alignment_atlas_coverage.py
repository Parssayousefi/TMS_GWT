"""Read-only atlas/FOV coverage check for completed corrected Aachen runs."""
from pathlib import Path
import json

import nibabel as nib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/aachen_alignment_review/corrected_cohort"
DERIV = ROOT / "data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125"
SESSIONS = ["pre90", "post90", "pre120", "post120"]


def main():
    queue = json.loads((ROOT / "logs/aachen_corrected_20261125/queue_status.json").read_text())
    atlas = nib.load(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    labels = np.asarray(atlas.dataobj).astype(int)
    counts = np.bincount(labels.ravel(), minlength=455)[1:455]
    assert counts.min() > 0
    summaries = []
    for subject, completed in sorted(queue.get("finished", {}).items()):
        if not (completed.get("success") is True and completed.get("state", {}).get("ExitCode") == 0):
            continue
        folder = OUT / f"sub-{subject}"
        if not (folder / "metrics.json").exists():
            continue
        table = folder / "atlas_functional_coverage.csv"
        summary_path = folder / "atlas_functional_coverage_summary.json"
        if summary_path.exists():
            summaries.append(json.loads(summary_path.read_text()))
            continue
        rows = []
        for session in SESSIONS:
            stem = f"sub-{subject}_ses-{session}_task-rest"
            func = DERIV / f"sub-{subject}/ses-{session}/func"
            mask_image = nib.load(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz")
            ref_image = nib.load(func / f"{stem}_space-MNI152NLin2009cAsym_res-2_boldref.nii.gz")
            np.testing.assert_allclose(mask_image.affine, atlas.affine, atol=1e-4)
            np.testing.assert_allclose(ref_image.affine, atlas.affine, atol=1e-4)
            mask = np.asarray(mask_image.dataobj) > 0
            ref = np.asarray(ref_image.dataobj)
            mask_counts = np.bincount(labels[mask], minlength=455)[1:455]
            nonzero_counts = np.bincount(labels[ref > 0], minlength=455)[1:455]
            for i in range(454):
                rows.append({"subject": subject, "session": session, "roi": i+1,
                             "cortical": i < 400, "atlas_voxels": int(counts[i]),
                             "bold_brainmask_voxels": int(mask_counts[i]),
                             "nonzero_bold_reference_voxels": int(nonzero_counts[i]),
                             "brainmask_coverage": float(mask_counts[i]/counts[i]),
                             "nonzero_reference_coverage": float(nonzero_counts[i]/counts[i])})
        df = pd.DataFrame(rows)
        df.to_csv(table, index=False)
        summary = {"subject": subject, "sessions": len(SESSIONS),
                   "parcel_session_rows": len(df),
                   "zero_brainmask_cortical_parcel_sessions": int(((df.bold_brainmask_voxels == 0) & df.cortical).sum()),
                   "zero_brainmask_all_parcel_sessions": int((df.bold_brainmask_voxels == 0).sum()),
                   "zero_reference_parcel_sessions": int((df.nonzero_bold_reference_voxels == 0).sum()),
                   "minimum_brainmask_coverage": float(df.brainmask_coverage.min()),
                   "minimum_nonzero_reference_coverage": float(df.nonzero_reference_coverage.min()),
                   "mean_cortical_brainmask_coverage": float(df.loc[df.cortical, "brainmask_coverage"].mean()),
                   "limitation": "Reference support/brainmask coverage, not a temporal signal-quality or effect test."}
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    if summaries:
        pd.DataFrame(summaries).to_csv(OUT / "atlas_functional_coverage_summary.csv", index=False)


if __name__ == "__main__":
    main()
