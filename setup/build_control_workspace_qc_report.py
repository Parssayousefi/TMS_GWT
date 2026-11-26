"""Combine the completed control QC runs without fitting or changing any mask."""
from datetime import datetime, timezone
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy import stats

from qc_sham_workspace import ROOT, REFERENCE, plot_comparison, sha256


def main():
    base = ROOT / "output/control_workspace_qc"
    runs = {"SUDMEX": base / (base / "LATEST.txt").read_text().strip(),
            "Aachen": base / (base / "LATEST_aachen_corrected.txt").read_text().strip()}
    summaries = {name: json.loads((path / "summary.json").read_text()) for name, path in runs.items()}
    downstream = json.loads((ROOT / "logs/aachen_corrected_20261125/downstream/status.json").read_text())
    assert downstream["phase"] == "complete" and downstream["phiid_scans"] == 60
    hcp_path = REFERENCE / "hcp_rest_workspace_mask.npy"
    gradient_path = REFERENCE / "hcp_rest_synergy_redundancy_gradient.npy"
    original = Path("D:/projects/GWT/pipeline/results/group_analysis/hcp/rest")
    for path in [hcp_path, gradient_path]:
        assert sha256(path) == sha256(original / path.name)
    hcp = np.load(hcp_path).astype(bool)
    gradient = np.load(gradient_path)
    out = base / ("combined_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ"))
    out.mkdir()
    selections = [("SUDMEX", "sham_t0", "SUDMEX: sham arm before treatment"),
                  ("SUDMEX", "sham_t1", "SUDMEX: the same sham participants after two weeks"),
                  ("Aachen", "baseline_mean", "Aachen: corrected pretreatment visits averaged within participant")]
    figures, specs, titles = {}, {}, {}
    sources = [Path(__file__), ROOT / "setup/qc_sham_workspace.py", hcp_path, gradient_path,
               *[p / "summary.json" for p in runs.values()]]
    for dataset, condition, title in selections:
        path = runs[dataset] / condition / "regional_comparison.csv"
        sources.append(path)
        nodes = pd.read_csv(path).sort_values("roi")
        np.testing.assert_array_equal(nodes.roi, np.arange(1, 455))
        np.testing.assert_array_equal(nodes.hcp_workspace, hcp)
        mask = nodes.workspace.to_numpy(bool)
        spec = summaries[dataset]["conditions"][condition]
        np.testing.assert_allclose(2 * (mask & hcp).sum() / (mask.sum() + hcp.sum()), spec["dice_hcp"])
        np.testing.assert_allclose(stats.spearmanr(nodes.gradient, gradient).statistic, spec["gradient_rho_hcp"])
        key = dataset + "_" + condition
        figures[key], specs[key], titles[key] = {"workspace_mask": mask}, spec, title
    plot_comparison(figures, hcp, specs, out, dataset_label="TMS", title_labels=titles)
    labels = {"all_baseline": "All pretreatment scans", "sham_t0": "Matched sham: pretreatment",
              "sham_t1": "Matched sham: two weeks", "baseline_mean": "Two pretreatment visits averaged",
              "pre90": "Before 90% session", "pre120": "Before 120% session"}
    lines = ["# Control workspace quality check", "",
             "Workspaces were independently derived using the unchanged GWT core: arithmetic means of persistent synergy/redundancy matrices, node-strength ranks, and synergy rank greater than redundancy rank. The HCP mask was not used to construct or tune these maps.", "",
             "| Dataset | Condition | Participants | Workspace regions | Shared with HCP | Dice | Gradient rho |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for dataset, summary in summaries.items():
        for condition, s in summary["conditions"].items():
            lines.append(f"| {dataset} | {labels[condition]} | {s['n']} | {s['n_workspace']} | {s['shared_regions']} | {s['dice_hcp']:.3f} | {s['gradient_rho_hcp']:.3f} |")
    sham = summaries["SUDMEX"]["sham_pre_post"]
    aachen = summaries["Aachen"]["baseline_visit_similarity"]
    aachen_mean = summaries["Aachen"]["conditions"]["baseline_mean"]
    boot = aachen_mean["bootstrap_dice_95_percentiles"]
    loo = aachen_mean["leave_one_out_dice_hcp_range"]
    lines += ["", "Dice measures mask overlap (1 means identical); rho is the Spearman correlation of the unthresholded synergy-minus-redundancy rank gradients.", "",
              f"The same 18 sham participants contribute at both SUDMEX visits: between-visit Dice {sham['dice_pre_post']:.3f}, gradient rho {sham['gradient_rho_pre_post']:.3f}. Later open-label treatment is excluded.", "",
              f"Aachen has no sham arm. Its two prestimulation visits are paired within the same 15 participants; between-visit Dice {aachen['dice_pre90_pre120']:.3f}, gradient rho {aachen['gradient_rho_pre90_pre120']:.3f}. The combined map averages matrices within participant before averaging participants, with no time-series concatenation.", "",
              f"For the combined Aachen map, the central 95% participant-bootstrap range of Dice against HCP is {boot[0]:.3f}–{boot[1]:.3f}; omitting each participant in turn gives {loo[0]:.3f}–{loo[1]:.3f}. These are stability diagnostics with the reference held fixed.", "",
              "Our 235-region HCP mask and its gradient are byte-identical to the original GWT project files. They are one frozen reference, not two independent comparators. HCP and psychedelic calculations were not rerun.", "",
              "Each condition includes 2,000 participant bootstraps and 1,000 split halves. Aachen also includes every leave-one-participant-out map. These describe sampling stability conditional on the pipeline; no spatial significance test, tuned overlap threshold, or selected participant exclusion was used.", "",
              "Aachen uses the corrected branch only: all 15 participants and 60 runs passed anatomical review and output validation through deconvolution and PhiID. Native anatomical support retention improved from 44.2–52.9% to 99.63–100%. Regional functional coverage remains uneven, particularly in inferior frontal and temporal parcels. The SUDMEX audit found no comparable anatomical over-extraction, but also documented inferior functional-coverage limitations. Workspace resemblance does not certify the temporal signal quality of every parcel.", "",
              "The previous Aachen treatment-effect estimates remain on hold pending statistical reanalysis of these corrected inputs. No treatment-effect inference is made in this control QC.", "",
              f"![Control workspaces compared with the frozen HCP mask]({(out / 'control_workspace_comparison.png').as_posix()})", "",
              f"SUDMEX source run: `{runs['SUDMEX'].relative_to(ROOT).as_posix()}`.",
              f"Corrected Aachen source run: `{runs['Aachen'].relative_to(ROOT).as_posix()}`.",
              "Full anatomical review: `output/aachen_alignment_review/corrected_cohort/COHORT_FINDINGS.md`."]
    report = "\n".join(lines) + "\n"
    (out / "README.md").write_text(report, encoding="utf-8")
    (ROOT / "docs/TMS_CONTROL_WORKSPACE_QC_20261125.md").write_text(report, encoding="utf-8")
    (out / "provenance.json").write_text(json.dumps({
        "source_runs": {k: str(v) for k, v in runs.items()},
        "source_sha256": {str(p.resolve()): sha256(p) for p in sources},
        "new_numerical_analysis": False,
        "figure_checks": "All source masks, overlap values and rank correlations verified against saved summaries."}, indent=2), encoding="utf-8")
    (base / "LATEST_combined.txt").write_text(out.name, encoding="utf-8")
    print(str(out))


if __name__ == "__main__":
    main()
