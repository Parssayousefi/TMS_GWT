"""Publish existing QC-refitted analyses without continuous FD covariates.

The user's current model specification uses motion only to select eligible
scans. Existing permitted fits are reused byte-for-value; no inference is
selected by significance and no correction family is narrowed. Superseded
continuous-FD models remain in their historical output directories.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "output/tms_fd_reanalysis_20261126"
HASHES = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    path = path.resolve()
    HASHES[path.relative_to(ROOT).as_posix()] = sha(path)
    return pd.read_csv(path)


def json_input(path):
    HASHES[path.relative_to(ROOT).as_posix()] = sha(path)
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verified_source(folder):
    receipt = json_input(folder / "completion_receipt.json")
    if receipt.get("completed") is not True:
        raise ValueError(f"Unfinished source: {folder}")
    for relative, expected in receipt["output_sha256"].items():
        if sha(folder / relative) != expected:
            raise ValueError(f"Source output changed since refit: {relative}")
    return receipt


def targets():
    source, out = BASE / "targets", BASE / "targets_qc_only"
    verified_source(source)
    extra = json_input(source / "independent_refit_verification.json")
    for relative, expected in extra["source_sha256"].items():
        if sha(source / relative) != expected:
            raise ValueError(f"Independent target audit input changed: {relative}")
    out.mkdir(exist_ok=True)
    summaries = []
    for branch in ("registered", "strict"):
        original, destination = source / branch, out / branch
        destination.mkdir(exist_ok=True)
        inference = read(original / "inference.csv")
        keep = ((inference.dataset.eq("aachen_itbs") & inference.model.eq("paired_mean")) |
                (inference.dataset.eq("sudmex_tms") & inference.model.eq("baseline")))
        inference = inference.loc[keep].copy()
        if len(inference) != 90 or inference.model.str.contains("motion", case=False).any():
            raise ValueError("Unexpected permitted target inference inventory")
        inference.to_csv(destination / "inference.csv", index=False)
        for filename in ("aachen_intensity_descriptive.csv", "aachen_intensity_exact_sensitivity.csv",
                         "participant_visit_changes.csv", "input_visit_pairs.csv", "excluded_participant_visits.csv"):
            frame = read(original / filename)
            # The group-analysis views do not carry unused FD covariate fields.
            frame.drop(columns=[c for c in ("delta_fd", "mean_fd", "legacy_mean_fd", "retained_mean_fd") if c in frame]).to_csv(destination / filename, index=False)
        loo = read(original / "leave_one_out.csv")
        # analyse_targeted_information_flow.sensitivities uses ancova without fd
        # for SUDMEX and a paired mean for Aachen; these existing fits are valid.
        loo["model"] = loo.dataset.map({"aachen_itbs": "paired_mean_leave_one_out", "sudmex_tms": "baseline_adjusted_leave_one_out"})
        if loo.model.isna().any():
            raise ValueError("Unknown target leave-one-out model")
        loo.to_csv(destination / "leave_one_out.csv", index=False)
        diagnostics = read(original / "model_diagnostics.csv")
        diagnostics = diagnostics.loc[diagnostics.model.isin(["paired_mean", "baseline"])]
        diagnostics.to_csv(destination / "model_diagnostics.csv", index=False)
        for (dataset, definition), block in inference.groupby(["dataset", "definition"]):
            primary = block.loc[block.variant.eq("lag1_full")]
            correction = "p_maxT5" if dataset == "aachen_itbs" else "p_maxT10"
            summaries.append({"qc_policy": branch, "dataset": dataset, "definition": definition,
                              "model": primary.model.iloc[0], "n": int(primary.n.iloc[0]),
                              "n_active": int(primary.n_active.iloc[0]) if dataset == "sudmex_tms" else None,
                              "n_sham": int(primary.n_sham.iloc[0]) if dataset == "sudmex_tms" else None,
                              "correction": correction, "minimum_adjusted_p": float(primary[correction].min()),
                              "primary_survivors": int(primary[correction].lt(.05).sum())})
    pd.DataFrame(summaries).to_csv(out / "cohort_summary.csv", index=False)
    return out, {"analysis": "Existing target-flow inference with motion used for scan exclusion only",
                 "sudmex_model": "Post outcome ~ randomized arm + baseline outcome; HC3; no FD term",
                 "aachen_model": "Paired within-person average of both complete intensity changes; no FD term",
                 "tests_per_policy": 90, "family_changes": "None: original 5-measure/2-scope families and all three variants retained",
                 "omitted_models": "All continuous-FD target ANCOVAs and Aachen baseline+FD sensitivity fits",
                 "leave_one_out": "Original SUDMEX baseline-only and Aachen paired fits; no FD covariate",
                 "summaries": summaries}


def sensitivities():
    source, out = BASE / "sensitivities", BASE / "sensitivities_qc_only"
    verified_source(source)
    out.mkdir(exist_ok=True)
    filenames = ["estimator_sensitivity_results.csv", "temporal_component_results.csv", "network_results.csv",
                 "visit_consistency_results.csv", "participant_influence.csv", "nonparametric_checks.csv",
                 "method_interactions.csv", "signal_participant_changes.csv", "network_participant_changes.csv"]
    for branch in ("registered", "strict"):
        destination = out / branch
        destination.mkdir(exist_ok=True)
        for filename in filenames:
            frame = read(source / branch / filename)
            frame.drop(columns=[c for c in ("delta_fd", "mean_fd") if c in frame]).to_csv(destination / filename, index=False)
        for filename in ("cohort.json", "independent_signflip_verification.json"):
            write_json(destination / filename, json_input(source / branch / filename))
    family = read(source / "family_summary.csv")
    family = family.loc[family.table.ne("baseline_motion_sensitivity.csv")].copy()
    if family.analysis_family.str.contains("baseline_motion").any():
        raise ValueError("A continuous-FD sensitivity entered the current view")
    family.to_csv(out / "family_summary.csv", index=False)
    return out, {"analysis": "Cached corrected-Aachen representation/network checks after scan exclusion",
                 "continuous_motion_covariates": False,
                 "omitted_models": "The entire historical baseline-motion6 sensitivity block is omitted",
                 "family_changes": "None for retained analyses; exact maxT estimator12, temporal9, network35/8, visit12 and method8 remain intact",
                 "normalization": "Cached original corrected-baseline Phi reference; no claim of fresh calibration here",
                 "interpretation": "Method/robustness diagnostics, not cross-dataset causal replication"}


def main():
    completed = [targets(), sensitivities()]
    helper = ROOT / "setup/analyse_targeted_information_flow.py"
    HASHES[helper.relative_to(ROOT).as_posix()] = sha(helper)
    for relative, expected in HASHES.items():
        if sha(ROOT / relative) != expected:
            raise ValueError(f"Source changed while publishing views: {relative}")
    for out, details in completed:
        write_json(out / "current_model_receipt.json", {
            "created_utc": datetime.now(timezone.utc).isoformat(), "complete": True,
            "model_specification": "User requested motion QC for exclusions only, with no continuous FD covariates in group analyses",
            "operation": "Publish already computed permitted fits; no new outcome selection, p-value recomputation, or family narrowing",
            "historical_outputs_preserved": True, **details, "source_sha256": HASHES,
            "script_sha256": sha(Path(__file__)),
            "output_sha256": {path.relative_to(out).as_posix(): sha(path) for path in out.rglob("*")
                              if path.is_file() and path.name != "current_model_receipt.json"}})
        print("COMPLETE", out)


if __name__ == "__main__":
    main()
