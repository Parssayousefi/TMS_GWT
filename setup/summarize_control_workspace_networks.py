"""Descriptive network/coverage context for independently derived control maps."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NAMES = ["Visual", "Somatomotor", "Dorsal attention", "Ventral attention", "Limbic", "Frontoparietal", "Default mode", "Subcortical"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=["sudmex_tms", "aachen_itbs"], default="sudmex_tms")
    parser.add_argument("--run", type=Path)
    args = parser.parse_args()
    base = ROOT / "output/control_workspace_qc"
    latest = "LATEST_aachen_corrected.txt" if args.dataset == "aachen_itbs" else "LATEST.txt"
    run = args.run or base / (base / latest).read_text().strip()
    summary = json.loads((run / "summary.json").read_text())
    network_path = ROOT / "pipeline/results/group_analysis/hcp/rest/hcp_rest_rsn_labels.npy"
    networks = np.load(network_path, allow_pickle=True).astype(int)
    if args.dataset == "aachen_itbs":
        assert set(summary["conditions"]) == {"baseline_mean", "pre90", "pre120"}
        coverage_paths = [ROOT / "output/aachen_alignment_review/corrected_cohort" /
                          f"sub-{subject}/atlas_functional_coverage.csv"
                          for subject in summary["conditions"]["baseline_mean"]["subjects"]]
        coverage = pd.concat([pd.read_csv(p, dtype={"subject": str}) for p in coverage_paths], ignore_index=True)
        coverage = coverage.rename(columns={"roi": "parcel_id", "brainmask_coverage": "fraction_of_anatomical_parcel_in_bold_mask"})
    else:
        assert set(summary["conditions"]) == {"all_baseline", "sham_t0", "sham_t1"}
        coverage_paths = [ROOT / "output/sudmex_registration_audit/parcel_mask_coverage.csv"]
        coverage = pd.read_csv(coverage_paths[0], dtype={"subject": str})
    assert not coverage.duplicated(["subject", "session", "parcel_id"]).any()
    records = []
    for condition, spec in summary["conditions"].items():
        nodes = pd.read_csv(run / condition / "regional_comparison.csv")
        if args.dataset == "aachen_itbs":
            sessions = ["pre90", "pre120"] if condition == "baseline_mean" else [condition]
        else:
            sessions = ["t1" if condition == "sham_t1" else "t0"]
        subset = coverage.loc[coverage.subject.isin(spec["subjects"]) & coverage.session.isin(sessions)].copy()
        assert len(subset) == spec["n"] * 454 * len(sessions)
        assert subset.groupby(["subject", "session"]).size().eq(454).all()
        field = "fraction_of_anatomical_parcel_in_bold_mask"
        subset["below_half_mask_support"] = subset[field] < .5
        means = subset.groupby("parcel_id").agg(mean_bold_mask_coverage=(field, "mean"),
                    n_scans_below_half_mask_support=("below_half_mask_support", "sum"))
        nodes = nodes.merge(means, left_on="roi", right_index=True, validate="one_to_one")
        nodes["network"] = [NAMES[i] for i in networks]
        nodes.to_csv(run / condition / "regional_comparison_with_coverage.csv", index=False)
        for network, group in nodes.groupby("network", sort=False):
            mask = group.workspace.to_numpy(bool)
            hcp = group.hcp_workspace.to_numpy(bool)
            shared = int((mask & hcp).sum())
            records.append({"condition": condition, "network": network, "n_regions": len(group),
                            "control_workspace": int(mask.sum()), "hcp_workspace": int(hcp.sum()),
                            "shared": shared, "dice": 2 * shared / (mask.sum() + hcp.sum()) if mask.any() or hcp.any() else np.nan,
                            "mean_bold_mask_coverage": group.mean_bold_mask_coverage.mean(),
                            "regions_with_any_scan_below_half_mask_support": int((group.n_scans_below_half_mask_support > 0).sum())})
    results = pd.DataFrame(records)
    results.to_csv(run / "network_overlap_and_coverage.csv", index=False)
    paths = [Path(__file__), run / "summary.json", network_path, *coverage_paths]
    provenance = {"dataset": args.dataset,
                  "note": "Descriptive coverage context; no exclusions or inferential tests. Two Aachen baselines contribute equally per participant.",
                  "source_sha256": {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (run / "network_coverage_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    primary = "baseline_mean" if args.dataset == "aachen_itbs" else "all_baseline"
    print(results.loc[results.condition.eq(primary)].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
