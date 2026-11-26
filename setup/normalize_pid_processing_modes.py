"""Conditional-null sensitivity for the six disjoint PhiID processing modes.

The production estimator and frozen cohort outputs are read-only. Each mode
is formed in information units BEFORE its own conditional CDF is calculated.
The two established baseline-only VAR pools are reused after hash verification.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import sys
import types
import shutil

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "setup"), str(ROOT / "pipeline/scripts/analysis")]
import explore_pid_processing_modes as raw
import tms_numit as tn
import run_tms_numit_sensitivity as legacy

RAW = ROOT / "output/pid_processing_modes/20261126T130022_906563Z"
SOURCES = {
    "aachen_itbs": raw.CORE / "normalization",
    "sudmex_tms": ROOT / "output/tms_numit/20261125T155252_384543Z",
}
SCOPES = list(json.loads((RAW / "scope_pair_counts.json").read_text()))
ZERO_TOL = 1e-12
HASHES = {}


def track(path):
    HASHES[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def canonical_modes(atoms):
    # Resolve numerical cancellation at exact MMI ties. Signed values larger
    # than this absolute tolerance are retained, never rectified or clipped.
    return {m: np.where(abs(v) <= ZERO_TOL, 0., v) for m, v in raw.modes(atoms).items()}


def expanded_atoms(sigma):
    atoms = tn.pg.phiid_atoms(sigma)
    return {**atoms, **canonical_modes(atoms), "_atom_tdmi": atoms["tdmi"]}


def independent_modes(sigma):
    shape = np.asarray(sigma).shape[:-2]
    atoms = raw.independent_atoms(np.asarray(sigma).reshape(-1, 4, 4))
    values = canonical_modes(atoms)
    return {**{m: v.reshape(shape) for m, v in values.items()},
            "tdmi": sum(atoms.values()).reshape(shape)}


def runtime():
    """Isolated namespaces; no modifications to imported production modules."""
    null = types.ModuleType("mode_null_runtime")
    null.__file__ = tn.nn.__file__
    exec(compile(Path(tn.nn.__file__).read_text(), tn.nn.__file__, "exec"), null.__dict__)
    null.pg = types.SimpleNamespace(normalize_sigma=tn.pg.normalize_sigma,
                                    phiid_atoms=expanded_atoms)
    null.MEASURES = [*raw.MODES, *legacy.MEASURES, "_atom_tdmi"]

    def valid(pool, target, count, seed):
        rng = np.random.default_rng(seed)
        collected = {m: [] for m in raw.MODES}
        needed, rejected, generated, maximum = count, 0, 0, 0.
        for attempt in range(20):
            size = needed + 64
            result = null.null_atoms_for_targets(pool, np.array([target]), size, rng)
            error = np.maximum(result["_tdmi_err"][0], abs(result["_atom_tdmi"][0] - target))
            good = np.isfinite(error) & (error <= 1e-6)
            for m in null.MEASURES:
                good &= np.isfinite(result[m][0])
            ix = np.flatnonzero(good)[:needed]
            generated += size
            rejected += int((~good).sum())
            if len(ix):
                maximum = max(maximum, float(error[ix].max()))
                for m in raw.MODES:
                    collected[m].append(result[m][0, ix])
                needed -= len(ix)
            if needed == 0:
                return {m: np.concatenate(v) for m, v in collected.items()}, {
                    "generated": generated, "rejected": rejected,
                    "max_kept_tdmi_error": maximum, "replenishment_batches": attempt}
        raise RuntimeError("Could not obtain the specified number of valid common draws")

    helper = types.ModuleType("mode_calibration_runtime")
    helper.__file__ = legacy.__file__
    exec(compile(Path(legacy.__file__).read_text(), legacy.__file__, "exec"), helper.__dict__)
    helper.MEASURES = list(raw.MODES)
    helper.tms = types.SimpleNamespace(DATASETS=tuple(SOURCES))
    helper.tn = types.SimpleNamespace(valid_reference=valid, empirical_cdf=tn.empirical_cdf,
                                     conditional_cdf=tn.conditional_cdf,
                                     stable_phiid_atoms=independent_modes)
    helper.pg = types.SimpleNamespace(**{k: getattr(tn.pg, k) for k in dir(tn.pg) if not k.startswith("__")})
    helper.pg.phiid_atoms = expanded_atoms
    return helper, null


def preparation(out, helper, null):
    cache = out / "cache"
    cache.mkdir(exist_ok=True)
    prepared = {"tdmi_ranges": {}, "calibration_targets": {}, "baseline_pool_sizes": {}}
    checks = []
    for dataset, source in SOURCES.items():
        record = json.loads(track(source / "cache/prepared.json").read_text())
        files = {p: h for p, h in record["input_sha256"].items()
                 if f"dataset-{dataset}_" in p and p.endswith("_timeseries_deconv.npy")}
        assert len(files) == (60 if dataset == "aachen_itbs" else 80)
        for rel, digest in files.items():
            path = track(ROOT / rel)
            assert HASHES[str(path.relative_to(ROOT))] == digest, path
            if dataset == "aachen_itbs":
                assert "reruns" in path.parts and "aachen_skullstrip_skip_20261125" in path.parts
        source_pool = track(source / "cache" / f"{dataset}_baseline_A_pool.npy")
        pool = np.load(source_pool)
        selected = np.load(track(source / "cache/pool_edge_indices.npy"))
        edges = np.array(np.triu_indices(454, k=1)).T[selected]
        # Reconstruct the baseline-only pool, including the original order.
        pieces = []
        scans = pd.read_csv(raw.CORE / "raw_scan_metrics.csv").sort_values(["dataset", "subject_id", "session"])
        plan = pd.read_csv(raw.CORE / "raw_paired_changes.csv")
        baseline_keys = set(zip(plan.subject_id, plan.baseline_session))
        scans = scans.loc[[(s.subject_id, s.session) in baseline_keys for s in scans.itertuples()]]
        for scan in scans.loc[scans.dataset.eq(dataset) & scans.session.isin(["t0", "pre90", "pre120"])].itertuples():
            base = raw.CORRECTED if dataset == "aachen_itbs" else ROOT
            path = base / "pipeline/results/timeseries" / f"{scan.prefix}_timeseries_deconv.npy"
            assert str(path.relative_to(ROOT)) in files
            sig = tn.pg.lagged_sigma_from_timeseries(np.load(path), 1, edges)
            a, _, _ = tn.pg.var1_from_sigma(sig)
            pieces.append(a[np.max(abs(np.linalg.eigvals(a)), axis=1) < 1])
        np.testing.assert_allclose(np.concatenate(pieces), pool, rtol=0, atol=1e-12)
        shutil.copyfile(source_pool, cache / source_pool.name)
        for key in prepared:
            prepared[key][dataset] = record[key][dataset]
        # Identical draws in the unchanged null engine reproduce the five
        # existing composites exactly; only requested summaries are extended.
        old = tn.nn.null_atoms_for_targets(pool, np.array([.8]), 2048, np.random.default_rng(9231))
        new = null.null_atoms_for_targets(pool, np.array([.8]), 2048, np.random.default_rng(9231))
        for m in [*legacy.MEASURES, "_tdmi_err"]:
            np.testing.assert_array_equal(old[m], new[m])
        np.testing.assert_allclose(sum(new[m] for m in raw.MODES), new["_atom_tdmi"], atol=3e-7, rtol=1e-6)
        checks.append({"dataset": dataset, "baseline_pool": len(pool), "source_scans": len(files),
                       "reconstructed_baseline_pool": True, "identical_legacy_draws": True})
    raw.save_json(cache / "prepared.json", prepared)
    raw.save_json(out / "input_verification.json", checks)
    return prepared


def scopes_and_scans():
    plan = pd.read_csv(track(raw.CORE / "raw_paired_changes.csv"))
    scans = pd.read_csv(track(raw.CORE / "raw_scan_metrics.csv"))
    wanted = set(zip(plan.subject_id, plan.baseline_session)) | set(zip(plan.subject_id, plan.followup_session))
    scans = scans.loc[[(r.subject_id, r.session) in wanted for r in scans.itertuples()]].sort_values(["dataset", "subject_id", "session"])
    assert len(scans) == 140
    ref = raw.tms.load_reference(ROOT)
    for path in ref["input_files"]:
        track(path)
    edges = np.array(np.triu_indices(454, 1)).T
    i, j = edges.T
    ws, lab = ref["workspace"], ref["modules"]
    striatum = ref["labels"].Structure_Name.eq("Striatum").to_numpy()
    def between(a, b):
        return (a[i] & b[j]) | (b[i] & a[j])
    scopes = {"whole_brain": np.ones(len(i), bool), "within_workspace": ws[i] & ws[j],
              "outside_workspace": ~ws[i] & ~ws[j], "workspace_boundary": ws[i] != ws[j],
              "FPN_striatum": between(lab == 5, striatum), "FPN_DAN": between(lab == 5, lab == 2),
              "FPN_DMN": between(lab == 5, lab == 6)}
    assert {k: int(v.sum()) for k, v in scopes.items()} == json.loads((RAW / "scope_pair_counts.json").read_text())
    return scans, plan, edges, scopes


def observed(scan, edges):
    base = raw.CORRECTED if scan.dataset == "aachen_itbs" else ROOT
    ts = np.load(base / "pipeline/results/timeseries" / f"{scan.prefix}_timeseries_deconv.npy")
    sigma = tn.pg.lagged_sigma_from_timeseries(ts, 1, edges)
    atoms = expanded_atoms(sigma)
    return atoms


def interpolation(out, scans, edges, helper, refs):
    rng = np.random.default_rng(2026112619)
    rows = []
    for dataset in SOURCES:
        pool = np.load(out / "cache" / f"{dataset}_baseline_A_pool.npy")
        block = scans.loc[scans.dataset.eq(dataset)]
        queries = []
        for scan in block.iloc[np.linspace(0, len(block)-1, 8, dtype=int)].itertuples():
            atoms = observed(scan, edges)
            tdmi = atoms["tdmi"]
            for k in [tdmi.argmin(), tdmi.argmax(), *rng.choice(len(tdmi), 2, replace=False)]:
                queries.append((tdmi[k], {m: atoms[m][k] for m in raw.MODES}))
        for seed in helper.SEEDS:
            grid, tables = refs[(dataset, seed)]
            for k, (target, values) in enumerate(queries):
                direct, _ = helper.tn.valid_reference(pool, target, helper.COUNT, seed)
                for m in raw.MODES:
                    expected = float(tn.empirical_cdf(np.sort(direct[m]), values[m]))
                    actual = float(tn.conditional_cdf(grid, tables[m], np.array([target]), np.array([values[m]]))[0])
                    rows.append(dict(dataset=dataset, seed=seed, query=k, mode=m, tdmi=target,
                                     direct=expected, interpolated=actual, absolute_error=abs(actual-expected)))
            print(f"Validated interpolation: {dataset}, seed {seed}", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "interpolation_validation.csv", index=False)
    assert frame.absolute_error.max() < .002, frame.sort_values("absolute_error").tail()
    return float(frame.absolute_error.max())


def extract(out, scans, edges, scopes, refs, helper):
    rows = []
    for k, scan in enumerate(scans.itertuples()):
        atoms = observed(scan, edges)
        results = {str(seed): {} for seed in helper.SEEDS}
        for m in raw.MODES:
            values = [tn.conditional_cdf(*[refs[(scan.dataset, seed)][0], refs[(scan.dataset, seed)][1][m]],
                                         atoms["tdmi"], atoms[m]) for seed in helper.SEEDS]
            for s, mask in scopes.items():
                for seed, v in zip(helper.SEEDS, values):
                    results[str(seed)][f"{s}__{m}"] = float(v[mask].mean())
        results["combined"] = {m: np.mean([v[m] for v in results.values()]) for m in next(iter(results.values()))}
        metadata = {n: getattr(scan, n) for n in ["dataset", "subject_id", "session", "arm", "mean_fd"]}
        rows.extend({**metadata, "reference_seed": seed, **v} for seed, v in results.items())
        if (k+1) % 10 == 0:
            print(f"Normalized six modes: {k+1}/140 scans", flush=True)
    frame = pd.DataFrame(rows)
    assert np.isfinite(frame[[f"{s}__{m}" for s in scopes for m in raw.MODES]]).all().all()
    frame.to_csv(out / "scan_modes_all_seeds.csv", index=False)
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    out = args.resume or ROOT / "output/pid_modes_normalized" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True, exist_ok=True)
    print("OUTPUT:", out, flush=True)
    for p in [Path(__file__), Path(raw.__file__), Path(legacy.__file__), Path(tn.__file__), Path(tn.nn.__file__), Path(tn.pg.__file__),
              RAW / "analysis_specification.json", RAW / "verification.json"]:
        track(p)
    spec = {"status": "Exploratory normalization sensitivity; specifications fixed before normalized treatment results",
            "created_utc": datetime.now(timezone.utc).isoformat(), "modes": raw.MODES, "scopes": SCOPES,
            "normalization": "Mode sum in nats then its own conditional CDF at each pair's TDMI, then mean across pairs. Never sum percentile atoms.",
            "pools": {k: str(v.relative_to(ROOT)) for k, v in SOURCES.items()},
            "nulls": "32768 common valid draws per .005-nat TDMI grid point and per seed; two seeds; baseline-only cohort-specific VAR pools",
            "numerics": "Common rejection on finite values and TDMI error <=1e-6; values with absolute magnitude <=1e-12 nats canonicalized to zero to handle MMI ties",
            "calibration": "Independent scipy Lyapunov/root draws and LAPACK atom inversion; 128 population/finite-record replicates at three TDMI levels and two scan lengths per cohort; interpolation check tolerance .002",
            "inference": "Same matched Aachen person averages and SUDMEX baseline+motion HC3 active-sham models as raw six-mode analysis. Holm6 per scope, Holm42 across scopes; Aachen exact maxT6/maxT42; replication conjunction Holm6/Holm42.",
            "seed_sensitivity": "Both numerical seeds reported separately, not independent replications. Raw and normalized branches remain visible.",
            "limits": "Conditional percentiles are not nats and the normalized modes are not an additive decomposition. Normalization does not establish a shared intervention effect or neural causation."}
    if args.resume:
        previous = json.loads((out / "run_source_sha256.json").read_text())
        for rel, digest in previous.items():
            assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest, rel
    else:
        raw.save_json(out / "analysis_specification.json", spec)
    helper, null = runtime()
    prepared = preparation(out, helper, null)
    scans, plan, edges, scopes = scopes_and_scans()
    raw.save_json(out / "run_source_sha256.json", HASHES)
    refs = helper.build_references(out, prepared)
    error = interpolation(out, scans, edges, helper, refs)
    calibration = helper.calibrate_finite_records(out, prepared, refs)
    frame = extract(out, scans, edges, scopes, refs, helper)
    result_frames, conjunction_frames, change_frames = [], [], []
    for seed in [str(s) for s in helper.SEEDS] + ["combined"]:
        destination = out if seed == "combined" else out / seed
        destination.mkdir(exist_ok=True)
        results, cc = raw.inference(frame.loc[frame.reference_seed.eq(seed)], plan, SCOPES, destination)
        result_frames.append(results.assign(reference_seed=seed))
        conjunction_frames.append(cc.assign(reference_seed=seed))
        change_frames.append(pd.read_csv(destination / "participant_changes.csv").assign(reference_seed=seed))
    pd.concat(result_frames).to_csv(out / "mode_treatment_results_all_seeds.csv", index=False)
    pd.concat(conjunction_frames).to_csv(out / "cross_dataset_conjunction_all_seeds.csv", index=False)
    seed_rows = []
    for dataset in SOURCES:
        for kind, data, keys, cols in [
            ("scans", frame, ["subject_id", "session"], [f"{s}__{m}" for s in scopes for m in raw.MODES]),
            ("changes", pd.concat(change_frames), ["subject_id", "comparison"], [f"{s}__{m}_delta" for s in scopes for m in raw.MODES])]:
            selected = data.loc[data.dataset.eq(dataset)]
            a, b = [selected.loc[selected.reference_seed.eq(str(seed))].set_index(keys).sort_index() for seed in helper.SEEDS]
            assert a.index.equals(b.index)
            for m in cols:
                difference = a[m] - b[m]
                seed_rows.append(dict(dataset=dataset, kind=kind, metric=m, rms_seed_difference=np.sqrt(np.mean(difference**2)),
                                      max_absolute_seed_difference=abs(difference).max(), data_sd=((a[m]+b[m])/2).std()))
    pd.DataFrame(seed_rows).to_csv(out / "monte_carlo_stability.csv", index=False)
    for rel, digest in HASHES.items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest, rel
    raw.save_json(out / "verification.json", {"passed": True, "scans": len(scans), "source_files_unchanged": len(HASHES),
                    "max_interpolation_error": error, "calibration": calibration,
                    "shared_holm6_results": int((conjunction_frames[-1].p_holm6 < .05).sum()),
                    "shared_holm42_results": int((conjunction_frames[-1].p_holm42 < .05).sum()),
                    "note": "Numerical verification under the specified null prior; not proof that BOLD follows a VAR(1) model."})
    (out.parent / "LATEST.txt").write_text(out.name + "\n", encoding="utf-8")
    print("COMPLETE", out, flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
