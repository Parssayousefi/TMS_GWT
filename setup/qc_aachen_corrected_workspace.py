"""Control-workspace QC from the fully corrected Aachen branch only.

Two pretreatment visits are checked separately and averaged within participant.
HCP/GWT reference outputs remain frozen read-only inputs.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import importlib.util
import json
import warnings

import numpy as np
import pandas as pd
from scipy import stats

from qc_sham_workspace import ROOT, REFERENCE, sha256, subset_workspace, summarize, plot_comparison

STAGE = ROOT / "reruns/aachen_skullstrip_skip_20261125/pipeline"
SUBJECTS = ["01", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--anatomy-release", type=Path, required=True)
    args = parser.parse_args()
    release = json.loads(args.anatomy_release.read_text(encoding="utf-8-sig"))
    assert release.get("passed") is True, "Full corrected-cohort anatomical QA must pass first"
    downstream = json.loads((ROOT / "logs/aachen_corrected_20261125/downstream/status.json").read_text())
    assert downstream["phase"] == "complete" and downstream["hrf_scans"] == downstream["phiid_scans"] == 60
    assert Path(downstream["anatomy_release"]).resolve() == args.anatomy_release.resolve()
    assert downstream["anatomy_release_sha256"] == sha256(args.anatomy_release)
    queue = json.loads((ROOT / "logs/aachen_corrected_20261125/queue_status.json").read_text())
    assert queue["status"] == "complete" and set(queue["finished"]) == set(SUBJECTS)
    assert all(r["success"] for r in queue["finished"].values())
    script = STAGE / "scripts/analysis/05_group_average_workspace.py"
    assert sha256(script) == sha256(ROOT / "pipeline/scripts/analysis/05_group_average_workspace.py")
    spec = importlib.util.spec_from_file_location("gwt_corrected_aachen_workspace", script)
    gwt = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gwt)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        pre90 = gwt.compute_workspace_for_dataset("aachen_itbs", STAGE, task="pre90")
        pre120 = gwt.compute_workspace_for_dataset("aachen_itbs", STAGE, task="pre120")
    for result in [pre90, pre120]:
        assert set(result["subjects"]) == set(SUBJECTS), "Unexpected baseline cohort; inspect exclusions/completeness"
    pre90 = subset_workspace(gwt, pre90, SUBJECTS)
    pre120 = subset_workspace(gwt, pre120, SUBJECTS)
    combined = dict(pre90)
    combined["control_task"] = "pre90_pre120_mean"
    combined["all_matrices"] = {
        key: (pre90["all_matrices"][key] + pre120["all_matrices"][key]) / 2
        for key in pre90["all_matrices"]}
    combined = subset_workspace(gwt, combined, SUBJECTS)
    results = {"baseline_mean": combined, "pre90": pre90, "pre120": pre120}
    hcp_path = REFERENCE / "hcp_rest_workspace_mask.npy"
    gradient_path = REFERENCE / "hcp_rest_synergy_redundancy_gradient.npy"
    original = Path("D:/projects/GWT/pipeline/results/group_analysis/hcp/rest")
    for path in [hcp_path, gradient_path]:
        assert sha256(path) == sha256(original / path.name)
    hcp = np.load(hcp_path).astype(bool)
    gradient = np.load(gradient_path)
    out = ROOT / "output/control_workspace_qc" / ("aachen_corrected_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ"))
    out.mkdir(parents=True)
    summaries = {}
    for name, result in results.items():
        folder = out / name
        folder.mkdir()
        gwt.save_single_dataset_results(result, folder)
        summaries[name] = summarize(gwt, result, hcp, gradient, folder)
        # Uniform influence QC, motivated by a baseline intensity-gradient flag.
        # All participants are omitted in turn; no omission selects the main map.
        strengths = {key: np.nansum(result["all_matrices"][key], axis=2)
                     for key in ["synergy", "redundancy"]}
        influence = []
        for index, subject in enumerate(SUBJECTS):
            keep = np.arange(len(SUBJECTS)) != index
            mask, rank_gradient, _, _ = gwt.identify_workspace(
                strengths["synergy"][keep].mean(0), strengths["redundancy"][keep].mean(0))
            influence.append({"omitted_subject": subject, "n": int(keep.sum()),
                              "workspace_regions": int(mask.sum()),
                              "dice_hcp": gwt.dice_coefficient(mask, hcp),
                              "gradient_rho_hcp": stats.spearmanr(rank_gradient, gradient).statistic,
                              "dice_full_cohort": gwt.dice_coefficient(mask, result["workspace_mask"])})
        influence = pd.DataFrame(influence)
        influence.to_csv(folder / "leave_one_participant_out.csv", index=False)
        summaries[name]["leave_one_out_dice_hcp_range"] = [float(influence.dice_hcp.min()), float(influence.dice_hcp.max())]
        summaries[name]["leave_one_out_dice_full_cohort_range"] = [float(influence.dice_full_cohort.min()), float(influence.dice_full_cohort.max())]
    visit_similarity = {
        "dice_pre90_pre120": gwt.dice_coefficient(pre90["workspace_mask"], pre120["workspace_mask"]),
        "gradient_rho_pre90_pre120": stats.spearmanr(pre90["gradient"], pre120["gradient"]).statistic}
    sources = [Path(__file__), ROOT / "setup/qc_sham_workspace.py", script, args.anatomy_release, hcp_path, gradient_path,
               STAGE / "configs/datasets/aachen_itbs.yaml", STAGE / "results/qc/aachen_itbs/exclusions.csv"]
    for subject in SUBJECTS:
        for session in ["pre90", "pre120"]:
            for measure in ["synergy", "redundancy", "phi", "phi_r"]:
                sources.append(STAGE / "results/phiid_matrices" /
                               f"dataset-aachen_itbs_sub-{subject}_ses-{session}_task-rest_run-1_{measure}.npy")
    summary = {"conditions": summaries, "baseline_visit_similarity": visit_similarity,
               "source_pipeline": str(STAGE), "anatomy_release": str(args.anatomy_release),
               "method": "Unchanged GWT core; persistent synergy/redundancy, arithmetic means, node-strength ranks, gradient > 0.",
               "averaging": "Average the two pretreatment PhiID matrices within each participant, then average participants; no time-series concatenation or post scans.",
               "interpretation": "Descriptive QC; no spatial p-values or tuned overlap threshold. Aachen has no sham condition.",
               "source_sha256": {str(p.resolve()): sha256(p) for p in sources}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plot_comparison(results, hcp, summaries, out, dataset_label="Aachen",
                    title_labels={"baseline_mean": "Two pretreatment visits averaged within participant",
                                  "pre90": "Before the 90% intensity visit", "pre120": "Before the 120% intensity visit"})
    lines = ["# Corrected Aachen control-workspace QC", "",
             "Both pretreatment scans are included with equal weight per participant. Separate visit maps provide a stability check. These are prestimulation controls; Aachen has no sham group. All source matrices are from the isolated corrected preprocessing branch.", "",
             "| Condition | Participants | Workspace regions | Shared with HCP | Dice | Gradient rho |",
             "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for name, s in summaries.items():
        lines.append(f"| {name} | {s['n']} | {s['n_workspace']} | {s['shared_regions']} | {s['dice_hcp']:.3f} | {s['gradient_rho_hcp']:.3f} |")
    lines += ["", f"Between-baseline Dice: {visit_similarity['dice_pre90_pre120']:.3f}; gradient Spearman rho: {visit_similarity['gradient_rho_pre90_pre120']:.3f}.", "",
              "The fixed HCP mask and gradient are byte-identical to the GWT project reference. The independently derived Aachen workspaces are QC outputs only; they do not replace the fixed reference for treatment-effect inference. No HCP or psychedelic calculations were rerun.", "",
              "Participant bootstraps, split-half outputs and uniform leave-one-participant-out maps quantify sampling stability and individual influence conditional on the pipeline. The full cohort remains the main map; no participant omission was selected to improve overlap. Spatial resemblance is descriptive and cannot itself certify registration or demonstrate a treatment effect.", "",
              "Independent review identified uneven inferior frontal/temporal functional coverage. Reference contrast in subjects 10 and 14 was investigated: saved dataflow and sampled actual volumes distinguish reference enhancement from the unchanged 4D intensity input. Both anatomical releases passed with this evidence; spatial bias and temporal-quality limits remain documented in the cohort audit. Nonzero reference support alone does not establish usable temporal signal.", "",
              f"Outputs: `{out.relative_to(ROOT).as_posix()}`."]
    report = "\n".join(lines) + "\n"
    (out / "README.md").write_text(report, encoding="utf-8")
    (ROOT / "docs/AACHEN_CORRECTED_WORKSPACE_QC_20261125.md").write_text(report, encoding="utf-8")
    (out.parent / "LATEST_aachen_corrected.txt").write_text(out.name, encoding="utf-8")
    print(json.dumps({"conditions": summaries, "baseline_visit_similarity": visit_similarity,
                      "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()
