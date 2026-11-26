"""Refit existing target-flow tests after exact-window motion exclusions.

No time series or atoms are recomputed. The corrected Aachen scalar release and
the SUDMEX-only part of the earlier target release are immutable inputs. Aachen
retains the original average of two complete intensity visits per participant.
The existing five-measure/two-scope families, variants, and random seed are used.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

sys.dont_write_bytecode = True
import analyse_targeted_information_flow as target

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "output/existing_tms_fd_audit_20261126"
AACHEN = ROOT / "output/corrected_aachen_target_flow/20261126T113156_656159Z"
SUDMEX = ROOT / "output/targeted_information_flow/20261125T193623_492441Z"
CORE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
OUT = ROOT / "output/tms_fd_reanalysis_20261126/targets"
SCAN_KEY = ["dataset", "subject_id", "session"]
PAIR_KEY = ["dataset", "definition", "variant", "subject_id", "comparison", "scope"]
TEST_KEY = ["dataset", "definition", "variant", "model", "scope", "metric"]
HASHES: dict[str, str] = {}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def track(path: Path, expected: str | None = None) -> Path:
    path = path.resolve()
    relative = path.relative_to(ROOT).as_posix()
    digest = sha(path)
    if expected is not None and digest != expected:
        raise ValueError(f"Released input hash mismatch: {relative}")
    if relative in HASHES and HASHES[relative] != digest:
        raise ValueError(f"Input changed during rerun: {relative}")
    HASHES[relative] = digest
    return path


def save_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def compare_frames(actual: pd.DataFrame, expected: pd.DataFrame, keys, label: str):
    """Compare all shared scalar columns, including original corrected p values."""
    if actual.duplicated(keys).any() or expected.duplicated(keys).any():
        raise ValueError(f"Duplicate parity keys: {label}")
    a = actual.set_index(keys).sort_index()
    b = expected.set_index(keys).sort_index()
    if not a.index.equals(b.index):
        raise ValueError(f"Parity keys differ: {label}")
    errors = {}
    for column in a.columns.intersection(b.columns):
        if pd.api.types.is_numeric_dtype(a[column]) and pd.api.types.is_numeric_dtype(b[column]):
            x, y = a[column].to_numpy(float), b[column].to_numpy(float)
            np.testing.assert_allclose(x, y, rtol=2e-7, atol=2e-10, equal_nan=True,
                                       err_msg=f"{label}: {column}")
            finite = np.isfinite(x) & np.isfinite(y)
            errors[column] = float(np.max(np.abs(x[finite] - y[finite]))) if finite.any() else 0.0
        else:
            if not a[column].fillna("<NA>").equals(b[column].fillna("<NA>")):
                raise ValueError(f"Parity column differs: {label}: {column}")
    return {"label": label, "rows": len(a), "max_absolute_errors": errors}


def load_sources():
    release = json.loads(track(AACHEN / "target_flow_release_receipt.json").read_text())
    if not release["analysis_complete"] or not release["required_numerical_checks_passed"]:
        raise ValueError("Corrected Aachen target release is incomplete")
    needed = ["scan_metrics.csv", "input_visit_pairs.csv", "participant_visit_changes.csv",
              "inference.csv", "aachen_baseline_motion_sensitivity.csv", "aachen_intensity_descriptive.csv",
              "leave_one_out.csv", "independent_verification.json", "analysis_specification.json"]
    for filename in needed:
        track(AACHEN / filename, release["files"][filename])
    sudmex_verified = json.loads(track(SUDMEX / "independent_verification.json").read_text())
    if sudmex_verified["status"] != "passed":
        raise ValueError("Prior SUDMEX target verification was not passed")
    scans_a = pd.read_csv(track(AACHEN / "scan_metrics.csv"))
    # Historical Aachen rows in this legacy table must never enter the analysis.
    scans_s = pd.read_csv(track(SUDMEX / "scan_metrics.csv")).query("dataset == 'sudmex_tms'").copy()
    scans = pd.concat([scans_a, scans_s], ignore_index=True)
    plan_a = pd.read_csv(track(AACHEN / "input_visit_pairs.csv"))
    plan_s = pd.read_csv(track(CORE / "raw_paired_changes.csv"))
    plan_s = plan_s.loc[plan_s.dataset.eq("sudmex_tms"), plan_a.columns].copy()
    if not plan_s.comparison.eq("t1_minus_t0").all():
        raise ValueError("SUDMEX target comparison must be the randomized t0 to t1 visit")
    plan = pd.concat([plan_a, plan_s], ignore_index=True)
    old_pairs = pd.concat([
        pd.read_csv(track(AACHEN / "participant_visit_changes.csv")),
        pd.read_csv(track(SUDMEX / "participant_visit_changes.csv")).query("dataset == 'sudmex_tms'")],
        ignore_index=True)
    old_inference = pd.concat([
        pd.read_csv(track(AACHEN / "inference.csv")),
        pd.read_csv(track(SUDMEX / "inference.csv")).query("dataset == 'sudmex_tms'")], ignore_index=True)
    if scans_a.subject_id.nunique() != 8 or scans_s.query("definition == 'common_baseline_target'").subject_id.nunique() != 40:
        raise ValueError("Unexpected original target cohort")
    if scans.duplicated(SCAN_KEY + ["definition", "variant", "side"]).any():
        raise ValueError("Duplicate target scan rows")
    for _, group in scans.groupby(SCAN_KEY + ["definition", "variant"]):
        if set(group.side) != {"l", "r"}:
            raise ValueError("Missing bilateral target support")
    if not set(scans.definition).issubset({"deposited_individual_corrected", "common_baseline_target", "individual_baseline_target"}):
        raise ValueError("Historical anatomical definition entered the rerun")
    for path in (Path(target.__file__), ROOT / "setup/analyse_tms_familywise.py"):
        track(path)
    return scans, plan, old_pairs, old_inference


def validate_pairs(pairs: pd.DataFrame) -> None:
    if pairs.duplicated(PAIR_KEY).any():
        raise ValueError("Duplicate participant target contrasts")
    for (dataset, _, _, _), block in pairs.groupby(["dataset", "definition", "variant", "subject_id"]):
        comparisons = {"post90_minus_pre90", "post120_minus_pre120"} if dataset == "aachen_itbs" else {"t1_minus_t0"}
        if set(block.comparison) != comparisons or len(block) != len(comparisons) * 2:
            raise ValueError("Incomplete intensity pairs or bilateral scopes")
        if set(block.scope) != set(target.SCOPES):
            raise ValueError("Target scopes are incomplete")
    numerical = pairs.select_dtypes(include=np.number)
    if not np.isfinite(numerical).all().all():
        raise ValueError("Non-finite target scalar")


def design_diagnostics(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, definition, variant, scope), block in pairs.groupby(["dataset", "definition", "variant", "scope"]):
        d = block.groupby("subject_id").mean(numeric_only=True) if dataset == "aachen_itbs" else block.set_index("subject_id")
        models = ["baseline", "baseline_motion"] if dataset == "sudmex_tms" and variant == "lag1_full" else ["baseline"]
        if dataset == "aachen_itbs":
            models = ["paired_mean"] + (["baseline_motion_sensitivity"] if variant == "lag1_full" else [])
        for metric in target.METRICS:
            for model in models:
                x = np.ones((len(d), 1))
                if model != "paired_mean":
                    x = np.column_stack([x, d[metric + "_baseline"] - d[metric + "_baseline"].mean()])
                    if dataset == "sudmex_tms":
                        x = np.column_stack([x, d.arm.eq("active").astype(float)])
                    if "motion" in model:
                        fd = d.delta_fd.to_numpy()
                        # Aachen helper's intercept estimates change at delta FD=0.
                        x = np.column_stack([x, fd if dataset == "aachen_itbs" else fd - fd.mean()])
                rank = np.linalg.matrix_rank(x)
                if rank != x.shape[1] or len(d) - rank < 2:
                    raise ValueError(f"Insufficient model rank or residual df: {dataset}/{definition}/{metric}/{model}")
                leverage = np.diag(x @ np.linalg.pinv(x))
                if leverage.max() >= 1 - 1e-9:
                    raise ValueError("Unit leverage in target model")
                rows.append({"dataset": dataset, "definition": definition, "variant": variant,
                             "scope": scope, "metric": metric, "model": model, "n": len(d),
                             "parameters": x.shape[1], "rank": int(rank), "residual_df": len(d) - int(rank),
                             "max_leverage": float(leverage.max()), "design_condition_number": float(np.linalg.cond(x))})
    return pd.DataFrame(rows)


def infer(pairs: pd.DataFrame, permutations: int) -> pd.DataFrame:
    validate_pairs(pairs)
    results = []
    for (dataset, definition, variant), block in pairs.groupby(["dataset", "definition", "variant"]):
        models = ["baseline", "baseline_motion"] if dataset == "sudmex_tms" and variant == "lag1_full" else ["baseline"]
        if dataset == "aachen_itbs":
            models = ["paired_mean"]
        for model in models:
            result = (target.aachen_stats(block) if dataset == "aachen_itbs" else
                      target.trial_stats(block, b=permutations, motion=model == "baseline_motion"))
            result["model"] = model
            result["dataset"] = dataset
            result["definition"] = definition
            result["variant"] = variant
            results.append(result)
        print(f"Tested {dataset} / {definition} / {variant}", flush=True)
    result = pd.concat(results, ignore_index=True)
    if len(result) != 110 or result.duplicated(TEST_KEY).any():
        raise ValueError("The original 110-test target inventory changed")
    return result


def select_qc(scans, plan, audit, cohorts, branch):
    selected = audit.loc[audit[f"{branch}_pass"], SCAN_KEY + ["retained_mean_fd"]]
    selected_scans = scans.merge(selected, on=SCAN_KEY, how="inner", validate="many_to_one")
    selected_scans["legacy_mean_fd"] = selected_scans.mean_fd
    selected_scans["mean_fd"] = selected_scans.retained_mean_fd
    allowed = set().union(*(set(ids) for ids in cohorts[branch].values()))
    selected_plan = plan.loc[plan.subject_id.isin(allowed)].copy()
    fd = audit.set_index(SCAN_KEY).retained_mean_fd
    scan_keys = set(map(tuple, selected[SCAN_KEY].to_numpy()))
    for index, pair in selected_plan.iterrows():
        baseline = (pair.dataset, pair.subject_id, pair.baseline_session)
        followup = (pair.dataset, pair.subject_id, pair.followup_session)
        if baseline not in scan_keys or followup not in scan_keys:
            raise ValueError("Primary cohort lists a scan failing QC")
        selected_plan.loc[index, "delta_fd"] = fd.loc[followup] - fd.loc[baseline]
    # Only complete planned participants contribute; excluded visits cannot enter make_pairs.
    selected_scans = selected_scans.loc[selected_scans.subject_id.isin(set(selected_plan.subject_id))].copy()
    pairs = target.make_pairs(selected_scans, selected_plan)
    validate_pairs(pairs)
    exclusions = plan.loc[~plan.subject_id.isin(allowed)].copy()
    exclusions["reason"] = "At least one required participant scan fails exact-window QC; both Aachen intensity pairs required"
    return selected_scans, selected_plan, pairs, exclusions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    args = parser.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(OUT):
        raise ValueError("Output must remain within the assigned target rerun directory")
    out.mkdir(parents=True, exist_ok=True)
    if (out / "completion_receipt.json").exists():
        raise FileExistsError("Completed rerun exists; use a child output directory for another run")
    permutations = 32768
    with threadpool_limits(limits=2):
        scans, plan, old_pairs, old_inference = load_sources()
        full = out / "full_sample_reproduction"
        full.mkdir(exist_ok=True)
        pairs = target.make_pairs(scans, plan)
        parity = [compare_frames(pairs, old_pairs, PAIR_KEY, "all target participant scalars")]
        design_diagnostics(pairs).to_csv(full / "model_diagnostics.csv", index=False)
        inference = infer(pairs, permutations)
        inference.to_csv(full / "inference.csv", index=False)
        parity.append(compare_frames(inference, old_inference, TEST_KEY, "all original target inference"))
        target.sensitivities(pairs, full)
        for filename, keys in [
            ("aachen_baseline_motion_sensitivity.csv", ["scope", "metric"]),
            ("aachen_intensity_descriptive.csv", ["scope", "comparison", "metric"]),
        ]:
            parity.append(compare_frames(pd.read_csv(full / filename), pd.read_csv(track(AACHEN / filename)), keys, filename))
        loo = pd.read_csv(full / "leave_one_out.csv")
        original_loo = pd.concat([pd.read_csv(track(AACHEN / "leave_one_out.csv")),
            pd.read_csv(track(SUDMEX / "leave_one_out.csv")).query("dataset == 'sudmex_tms'")], ignore_index=True)
        parity.append(compare_frames(loo, original_loo, ["dataset", "definition", "scope", "metric", "omitted"], "leave-one-out fits"))
        save_json(full / "parity.json", {"passed": True, "checks": parity})

        audit_dir = args.audit.resolve()
        receipt = json.loads(track(audit_dir / "audit_receipt.json").read_text())
        if receipt["state"] != "complete" or not receipt["all_193_windows_verified"]:
            raise ValueError("Exact-window audit is not complete")
        for filename in ("scan_audit.csv", "cohort_subjects.json", "pair_audit.csv"):
            track(audit_dir / filename, receipt["output_sha256"][filename])
        audit = pd.read_csv(audit_dir / "scan_audit.csv")
        cohorts = json.loads((audit_dir / "cohort_subjects.json").read_text())
        if len(audit) != 193 or audit.duplicated(SCAN_KEY).any() or not audit.bounds_verified.all():
            raise ValueError("Invalid exact-window whitelist")
        summaries = []
        for branch in ("registered", "strict"):
            dest = out / branch
            dest.mkdir(exist_ok=True)
            branch_scans, branch_plan, branch_pairs, exclusions = select_qc(scans, plan, audit, cohorts, branch)
            branch_scans.to_csv(dest / "scan_metrics.csv", index=False)
            branch_plan.to_csv(dest / "input_visit_pairs.csv", index=False)
            branch_pairs.to_csv(dest / "participant_visit_changes.csv", index=False)
            exclusions.to_csv(dest / "excluded_participant_visits.csv", index=False)
            diagnostics = design_diagnostics(branch_pairs)
            diagnostics.to_csv(dest / "model_diagnostics.csv", index=False)
            result = infer(branch_pairs, permutations)
            result.to_csv(dest / "inference.csv", index=False)
            target.sensitivities(branch_pairs, dest)
            for (dataset, definition), block in result.groupby(["dataset", "definition"]):
                primary = block.loc[block.variant.eq("lag1_full") & block.model.eq("paired_mean" if dataset == "aachen_itbs" else "baseline_motion")]
                p_column = "p_maxT5" if dataset == "aachen_itbs" else "p_maxT10"
                summaries.append({"branch": branch, "dataset": dataset, "definition": definition,
                                  "n": int(primary.n.iloc[0]),
                                  "n_active": int(primary.n_active.iloc[0]) if dataset == "sudmex_tms" else None,
                                  "n_sham": int(primary.n_sham.iloc[0]) if dataset == "sudmex_tms" else None,
                                  "primary_model": str(primary.model.iloc[0]), "primary_adjustment": p_column,
                                  "minimum_primary_adjusted_p": float(primary[p_column].min()),
                                  "primary_survivors": int(primary[p_column].lt(.05).sum())})
            print(f"Completed {branch}: " + str(summaries[-3:]), flush=True)
        pd.DataFrame(summaries).to_csv(out / "cohort_summary.csv", index=False)
        changed = [path for path, digest in HASHES.items() if sha(ROOT / path) != digest]
        if changed:
            raise ValueError(f"Source files changed during rerun: {changed}")
        outputs = {str(path.relative_to(out)).replace("\\", "/"): sha(path) for path in out.rglob("*") if path.is_file()}
        save_json(out / "completion_receipt.json", {
            "created_utc": datetime.now(timezone.utc).isoformat(), "completed": True,
            "original_result_parity_passed": True, "all_original_110_tests_refitted_per_branch": True,
            "new_preprocessing_or_atom_estimation": False, "source_covariate": "Exact retained-window mean FD",
            "scan_exclusion_policy": "Exclude complete participant from primary if any required visit fails; no salvage/censoring",
            "aachen_estimand": "Mean of the two matched intensity changes in the corrected anatomical target subset",
            "sudmex_estimand": "Randomized t0 to t1 active versus sham; common target primary, individual target sensitivity",
            "sudmex_target_definition": "Fixed original common anatomical target from all 21 published localization coordinates",
            "families": {"aachen": "Five measures per scope exact maxT5; Holm10 wider sensitivity",
                         "sudmex": "Five measures by two scopes shared-label permutation maxT10; maxT5 and Holm10 also retained",
                         "aachen_motion_sensitivity": "Holm10", "aachen_intensity_sensitivity": "Holm20"},
            "seed": target.SEED, "permutations": permutations, "summaries": summaries,
            "limitations": ["Aachen target samples are very small and lack sham control; nominal and covariate sensitivities are not a replicated causal effect.",
                            "SUDMEX label permutation assumes exchangeable treatment labels in retained imaging completers; original randomization blocks not reconstructed.",
                            "All reported confidence intervals are marginal, not simultaneous.",
                            "The original fresh direct MATLAB target-atom check was unavailable; this scalar rerun uses previously verified corrected production matrices and independent Python identities."],
            "source_sha256": HASHES, "script_sha256": sha(Path(__file__)), "output_sha256": outputs})
    print("COMPLETE", out, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
