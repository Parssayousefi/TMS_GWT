"""Derive SUDMEX sham workspaces with the unchanged GWT core functions.

This is descriptive QC, not a treatment-effect test or an atlas-fitting step.
Only t0 and the randomized t1 phase enter; later open-label scans do not.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import nibabel as nib
from nilearn import plotting
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "pipeline/scripts/analysis/05_group_average_workspace.py"
REFERENCE = ROOT / "pipeline/results/group_analysis/hcp/rest"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def subset_workspace(gwt, original, subjects):
    """Subset participants, then apply the production averaging and ranking."""
    assert len(original["subjects"]) == len(set(original["subjects"]))
    index = [original["subjects"].index(s) for s in subjects]
    result = dict(original)
    result["subjects"] = list(subjects)
    result["all_matrices"] = {k: v[index] for k, v in original["all_matrices"].items()}
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        for key, matrices in result["all_matrices"].items():
            result["group_" + key] = np.nanmean(matrices, axis=0)
    for key in ["synergy", "redundancy"]:
        result[key + "_strength"] = gwt.compute_node_strength(result["group_" + key])
    (result["workspace_mask"], result["gradient"], result["synergy_ranks"],
     result["redundancy_ranks"]) = gwt.identify_workspace(
        result["synergy_strength"], result["redundancy_strength"])
    independent = (pd.Series(result["synergy_strength"]).rank()
                   - pd.Series(result["redundancy_strength"]).rank()).to_numpy()
    np.testing.assert_array_equal(independent, result["gradient"])
    np.testing.assert_array_equal(independent > 0, result["workspace_mask"])
    return result


def summarize(gwt, result, hcp, hcp_gradient, out):
    n = len(result["subjects"])
    mask = result["workspace_mask"]
    strengths = {}
    offdiag = ~np.eye(454, dtype=bool)
    for key in ["synergy", "redundancy"]:
        matrices = result["all_matrices"][key]
        assert np.isfinite(matrices[:, offdiag]).all()
        strengths[key] = np.nansum(matrices, axis=2)
        np.testing.assert_allclose(strengths[key].mean(0), result[key + "_strength"], atol=1e-10)
    rng = np.random.default_rng(2026112541)
    bootstrap = []
    split = []
    membership = np.zeros(454)
    for draw in range(2000):
        ix = rng.integers(n, size=n)
        m, g, _, _ = gwt.identify_workspace(
            strengths["synergy"][ix].mean(0), strengths["redundancy"][ix].mean(0))
        membership += m
        bootstrap.append({"draw": draw, "dice_hcp": gwt.dice_coefficient(m, hcp),
                          "rho_hcp": stats.spearmanr(g, hcp_gradient).statistic})
    for draw in range(1000):
        ix = rng.permutation(n)
        masks = [gwt.identify_workspace(strengths["synergy"][j].mean(0),
                                        strengths["redundancy"][j].mean(0))[0]
                 for j in [ix[:n // 2], ix[n // 2:]]]
        split.append({"draw": draw, "dice_halves": gwt.dice_coefficient(*masks)})
    boot = pd.DataFrame(bootstrap)
    halves = pd.DataFrame(split)
    boot.to_csv(out / "participant_bootstrap.csv", index=False)
    halves.to_csv(out / "split_half_stability.csv", index=False)
    pd.DataFrame({"roi": np.arange(1, 455), "workspace": mask, "hcp_workspace": hcp,
                  "gradient": result["gradient"], "hcp_gradient": hcp_gradient,
                  "bootstrap_membership_probability": membership / 2000
                  }).to_csv(out / "regional_comparison.csv", index=False)
    shared = int((mask & hcp).sum())
    return {"n": n, "n_workspace": int(mask.sum()), "hcp_n_workspace": int(hcp.sum()),
            "shared_regions": shared, "dice_hcp": gwt.dice_coefficient(mask, hcp),
            "jaccard_hcp": shared / int((mask | hcp).sum()),
            "gradient_rho_hcp": stats.spearmanr(result["gradient"], hcp_gradient).statistic,
            "bootstrap_dice_95_percentiles": boot.dice_hcp.quantile([.025, .975]).tolist(),
            "split_half_dice_median": halves.dice_halves.median(),
            "subjects": result["subjects"]}


def plot_comparison(results, hcp, summaries, out, dataset_label="SUDMEX", title_labels=None):
    atlas = nib.load(ROOT / "pipeline/atlases/schaefer_400_tian_54_combined.nii.gz")
    labels = np.asarray(atlas.dataobj).astype(int)
    background = ROOT / "tmp/templateflow/tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz"
    colors = ["#24799e", "#dd8745", "#9a67ae"]
    fig = plt.figure(figsize=(13, 8.5))
    for row, (name, result) in enumerate(results.items()):
        mask = result["workspace_mask"]
        lookup = np.zeros(455, dtype=np.int16)
        lookup[1:][mask & hcp] = 1
        lookup[1:][mask & ~hcp] = 2
        lookup[1:][~mask & hcp] = 3
        img = nib.Nifti1Image(lookup[labels], atlas.affine)
        nib.save(img, out / (name + "_hcp_overlap.nii.gz"))
        bottom = .68 - row * .29
        plotting.plot_roi(img, bg_img=background, cmap=ListedColormap(colors), vmin=1, vmax=3,
                          display_mode="z", cut_coords=[-20, 0, 20, 40, 60], colorbar=False,
                          alpha=.8, figure=fig, axes=[.02, bottom, .96, .22], annotate=True,
                          draw_cross=False, black_bg=False)
        s = summaries[name]
        labels_for_title = title_labels or {"all_baseline": "All pretreatment scans", "sham_t0": "Sham arm: pretreatment",
                                           "sham_t1": "Sham arm: after two weeks"}
        title = labels_for_title[name]
        fig.text(.03, bottom + .235,
                 f"{title} (n={s['n']}): {s['shared_regions']} shared regions; Dice={s['dice_hcp']:.3f}",
                 fontsize=12)
    fig.suptitle(f"{dataset_label} control workspaces versus the frozen HCP/GWT mask", fontsize=15, y=.98)
    fig.legend(handles=[Patch(facecolor=c, label=t) for c, t in zip(colors,
               ["Shared workspace", f"{dataset_label} only", "HCP only"])],
               loc="lower center", bbox_to_anchor=(.5, .015), ncol=3, frameon=False)
    fig.savefig(out / "control_workspace_comparison.png", dpi=170)
    fig.savefig(out / "control_workspace_comparison.pdf")
    plt.close(fig)


def main():
    spec = importlib.util.spec_from_file_location("gwt_workspace_sham_qc", SCRIPT)
    gwt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gwt)
    clinical_path = ROOT / "data/master_clinical_effects.csv"
    participants_path = ROOT / "data/sudmex_tms/participants.tsv"
    clinical = pd.read_csv(clinical_path)
    clinical = clinical.loc[clinical.dataset.eq("sudmex_tms")]
    assert clinical.bids_subject_id.is_unique
    participants = pd.read_csv(participants_path, sep="\t").set_index("participant_id")
    candidates = clinical.loc[clinical.arm.eq("sham") & clinical.mri_pair_usable.eq(1)]
    subjects = sorted(candidates.bids_subject_id.str.removeprefix("sub-").tolist())
    assert len(subjects) == 18, "Unexpected control sample; inspect cohort before interpretation"
    assert participants.loc[["sub-" + s for s in subjects], "group"].eq(1).all()
    assert participants.loc[["sub-" + s for s in subjects], "exclusion"].eq(0).all()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        baseline = gwt.compute_workspace_for_dataset("sudmex_tms", ROOT / "pipeline", task="t0")
        post = gwt.compute_workspace_for_dataset("sudmex_tms", ROOT / "pipeline", task="t1")
    # Verify this wrapper exactly reproduces the production result before subsetting.
    same = subset_workspace(gwt, baseline, baseline["subjects"])
    for key in ["workspace_mask", "gradient", "synergy_strength", "redundancy_strength"]:
        np.testing.assert_array_equal(same[key], baseline[key])
    assert set(subjects) <= set(baseline["subjects"]) & set(post["subjects"])
    results = {"all_baseline": baseline,
               "sham_t0": subset_workspace(gwt, baseline, subjects),
               "sham_t1": subset_workspace(gwt, post, subjects)}
    hcp_path = REFERENCE / "hcp_rest_workspace_mask.npy"
    grad_path = REFERENCE / "hcp_rest_synergy_redundancy_gradient.npy"
    original = Path("D:/projects/GWT/pipeline/results/group_analysis/hcp/rest")
    for path in [hcp_path, grad_path]:
        assert sha256(path) == sha256(original / path.name)
    hcp = np.load(hcp_path).astype(bool)
    gradient = np.load(grad_path)
    out = ROOT / "output/control_workspace_qc" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True)
    summaries = {}
    for name, result in results.items():
        folder = out / name
        folder.mkdir()
        gwt.save_single_dataset_results(result, folder)
        summaries[name] = summarize(gwt, result, hcp, gradient, folder)
    direct = {"n_paired_sham": len(subjects),
              "dice_pre_post": gwt.dice_coefficient(results["sham_t0"]["workspace_mask"],
                                                     results["sham_t1"]["workspace_mask"]),
              "gradient_rho_pre_post": stats.spearmanr(results["sham_t0"]["gradient"],
                                                         results["sham_t1"]["gradient"]).statistic}
    sources = [SCRIPT, clinical_path, participants_path, hcp_path, grad_path]
    for task, result in [("t0", baseline), ("t1", results["sham_t1"])]:
        for subject in result["subjects"]:
            for measure in ["synergy", "redundancy", "phi", "phi_r"]:
                sources.append(ROOT / "pipeline/results/phiid_matrices" /
                               f"dataset-sudmex_tms_sub-{subject}_ses-{task}_task-rest_run-1_{measure}.npy")
    provenance = {str(p.relative_to(ROOT)): sha256(p) for p in sources}
    summary = {"conditions": summaries, "sham_pre_post": direct,
               "reference_identical_to_original_gwt": True,
               "method": "Unchanged GWT core: arithmetic group matrix means, node-strength ranks, synergy rank > redundancy rank.",
               "scope": "All baseline scans plus the same 18 sham participants at t0 and randomized-phase t1; no later open-label visits.",
               "interpretation": "Descriptive QC; no spatial p-values, no tuned overlap threshold, no proof that registration or every parcel is valid.",
               "source_sha256": provenance}
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plot_comparison(results, hcp, summaries, out)
    lines = ["# Control-condition workspace QC", "",
             "The GWT workspace algorithm was applied independently to the control data; the HCP reference was not used to construct these workspaces. The production functions were imported unchanged. A subset wrapper was verified against the full production result, and independent pandas ranks reproduce each mask.", "",
             "| Condition | Participants | Workspace regions | Shared with HCP | Dice | Gradient Spearman rho |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, s in summaries.items():
        lines.append(f"| {name} | {s['n']} | {s['n_workspace']} | {s['shared_regions']} | {s['dice_hcp']:.3f} | {s['gradient_rho_hcp']:.3f} |")
    lines += ["", f"The same {len(subjects)} sham participants contribute to both visits. Sham pre/post workspace Dice is {direct['dice_pre_post']:.3f}; the rank-gradient correlation is {direct['gradient_rho_pre_post']:.3f}. Later open-label active treatment is excluded.", "",
              "The 235-region reference mask and gradient are byte-identical to the files in the original GWT project. Those are therefore the same comparator. No HCP or psychedelic computation was rerun.", "",
              "Each condition includes 2,000 participant bootstrap draws and 1,000 random split halves. These estimate sample stability with the reference held fixed. No spatial significance claim or predefined pass threshold is implied.", "",
              "The independent full SUDMEX anatomical audit found no Aachen-like over-extraction, but inferior functional coverage and some focal alignment limitations remain. Workspace resemblance cannot validate every parcel or establish a treatment effect. See output/sudmex_registration_audit/SUDMEX_FULL_PREPROCESSING_QC.md.", "",
              "Aachen control workspaces remain pending corrected anatomy-dependent preprocessing. Old Aachen results must not be accepted based on apparent resemblance to HCP.", "",
              f"Output directory: `{out.relative_to(ROOT).as_posix()}`."]
    report = "\n".join(lines) + "\n"
    (out / "README.md").write_text(report, encoding="utf-8")
    (ROOT / "docs/TMS_CONTROL_WORKSPACE_QC_20261125.md").write_text(report, encoding="utf-8")
    (out.parent / "LATEST.txt").write_text(out.name, encoding="utf-8")
    print(json.dumps({"conditions": summaries, "sham_pre_post": direct, "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
