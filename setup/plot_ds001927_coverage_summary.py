"""Compact scientific QC plot from the already-audited pilot diagnostic means."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/ds001927_coverage_trace_20261126"
atlas = nib.load(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
template = nib.load(ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_desc-brain_mask.nii.gz").get_fdata() > 0
index = int(round(nib.affines.apply_affine(np.linalg.inv(atlas.affine), [0, 0, -20])[2]))
fig, axes = plt.subplots(2, 3, figsize=(10, 6.7))
for row, subject in enumerate(("s001", "s022")):
    stem = f"sub-{subject}_ses-fro_task-restpre"
    paths = [ROOT / "data/ds001927/derivatives/fmriprep_23_2_1" / f"sub-{subject}/anat/sub-{subject}_space-MNI152NLin2009cAsym_res-2_desc-preproc_T1w.nii.gz",
             OUT / f"{stem}_raw_mean_MNI.nii.gz", OUT / f"{stem}_preprocessed_mean.nii.gz"]
    for col, path in enumerate(paths):
        im = nib.load(path)
        np.testing.assert_allclose(im.affine, atlas.affine, atol=1e-5)
        values = im.get_fdata(dtype=np.float32)
        ax = axes[row, col]
        ax.imshow(values[:, :, index].T, origin="lower", cmap="gray", vmin=0,
                  vmax=np.percentile(values[values > 0], 98))
        ax.contour(template[:, :, index].T, levels=[.5], colors=["#20afc2"], linewidths=.7)
        ax.set_title(["Anatomy", "Raw BOLD: display warp only", "Preprocessed BOLD"][col], fontsize=11)
        if col == 0:
            ax.text(.03,.96,subject,transform=ax.transAxes,va="top",color="white",fontsize=11)
        ax.axis("off")
fig.suptitle("Inferior signal loss is already visible in the raw scans", fontsize=15)
fig.text(.5,.015,"Frontal-session baseline examples; all 12 pilot runs were checked. MNI z ≈ −20 mm.\nCyan: template outline. The raw mean has not undergone motion correction.",ha="center",fontsize=10)
fig.tight_layout(rect=(0,.07,1,.95))
fig.savefig(OUT / "coverage_summary.png",dpi=150)
plt.close(fig)
