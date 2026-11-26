"""Rebuild baseline-only normalization after verified retained-frame FD QC.

The unchanged five-composite and six-mode acceptance rules are preserved. Each
CDF is applied pairwise before averaging; role metrics use the original frozen
HCP mask. No HCP or psychedelic data are reprocessed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
import sys
import time
import types

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "setup"), str(ROOT / "pipeline/scripts/analysis")]
import normalize_pid_processing_modes as mode_norm
import rerun_pid_modes_qc_only as rerun
import tms_numit as tn
import tms_publication as tms
import analyse_tms_workspace_scope as ws

OUT = ROOT / "output/tms_fd_reanalysis_20261126/recalibration"
AUDIT = ROOT / "output/existing_tms_fd_audit_20261126"
CORE = rerun.CORE
SOURCES = mode_norm.SOURCES
MODE_SOURCE = rerun.NORMALIZED
MODES = list(mode_norm.raw.MODES)
COMPOSITES = list(tms.MATRIX_NAMES)
MEASURES = COMPOSITES + MODES
SEEDS = [2026112511, 2026112512]
HASHES = {}
KEY = ["dataset", "subject_id", "session"]


def track(path):
    path = Path(path).resolve()
    key = str(path.relative_to(ROOT))
    if key not in HASHES:
        HASHES[key] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read_json(path):
    return json.loads(track(path).read_text(encoding="utf-8-sig"))


def read_csv(path, **kwargs):
    return pd.read_csv(track(path), **kwargs)


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + "\n", encoding="utf-8")


def source_time_series(scan):
    base = mode_norm.raw.CORRECTED if scan.dataset == "aachen_itbs" else ROOT
    return track(base / "pipeline/results/timeseries" / f"{scan.prefix}_timeseries_deconv.npy")


def runtime():
    helper, null = mode_norm.runtime()
    old_mode_valid = helper.tn.valid_reference

    def dual_valid(pool, target, count, seed):
        # Identical initial random draws, distinct unchanged family acceptance rules.
        values = null.null_atoms_for_targets(pool, np.array([target]), count + 64, np.random.default_rng(seed))
        basic_error = values["_tdmi_err"][0]
        composite_good = np.isfinite(basic_error) & (basic_error <= 1e-6)
        for measure in COMPOSITES:
            composite_good &= np.isfinite(values[measure][0])
        mode_error = np.maximum(basic_error, abs(values["_atom_tdmi"][0] - target))
        mode_good = np.isfinite(mode_error) & (mode_error <= 1e-6)
        for measure in null.MEASURES:
            mode_good &= np.isfinite(values[measure][0])
        answer, qc = {}, {}
        for family, names, good, error, original_valid in [
            ("composite", COMPOSITES, composite_good, basic_error, tn.valid_reference),
            ("mode", MODES, mode_good, mode_error, old_mode_valid),
        ]:
            indices = np.flatnonzero(good)[:count]
            if len(indices) != count:
                draws, checks = original_valid(pool, target, count, seed)
                answer.update(draws)
                qc[family + "_fallback"] = True
                qc.update({family + "_" + key: value for key, value in checks.items()})
            else:
                answer.update({m: values[m][0, indices] for m in names})
                qc.update({family + "_fallback": False, family + "_generated": count + 64,
                           family + "_rejected": int((~good).sum()),
                           family + "_max_kept_tdmi_error": float(error[indices].max())})
        return answer, qc

    def independent_all(sigma):
        return {**tn.stable_phiid_atoms(sigma), **mode_norm.independent_modes(sigma)}

    helper.MEASURES = MEASURES
    helper.tn = types.SimpleNamespace(valid_reference=dual_valid, empirical_cdf=tn.empirical_cdf,
                                     conditional_cdf=tn.conditional_cdf, stable_phiid_atoms=independent_all)
    assert helper.SEEDS == SEEDS and helper.COUNT == 32768 and helper.STEP == .005
    return helper, old_mode_valid


def load_previous_refs(dataset):
    refs = {}
    for seed in SEEDS:
        old_folder = SOURCES[dataset] / "references" / dataset / str(seed)
        mode_folder = MODE_SOURCE / "references" / dataset / str(seed)
        grid = np.load(track(old_folder / "tdmi_grid.npy"))
        np.testing.assert_array_equal(grid, np.load(track(mode_folder / "tdmi_grid.npy")))
        values = {m: np.load(track((old_folder if m in COMPOSITES else mode_folder) / f"{m}.npy"), mmap_mode="r")
                  for m in MEASURES}
        assert all(v.shape == (len(grid), 32768) for v in values.values())
        refs[(dataset, seed)] = grid, values
    return refs


def bind_reference_cache(destination, prepared, datasets, helper, old_refs):
    """Refuse reuse when a pool, mathematical helper, grid, or draw rule changed."""
    binding = dict(measures=MEASURES, count=helper.COUNT, seeds=SEEDS, step=helper.STEP,
                   dual_acceptance_sha256=hashlib.sha256(inspect.getsource(helper.tn.valid_reference).encode()).hexdigest(),
                   mathematical_source_sha256={str(Path(p).relative_to(ROOT)): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                       for p in [mode_norm.__file__, mode_norm.raw.__file__, mode_norm.legacy.__file__,
                                 tn.__file__, tn.nn.__file__, tn.pg.__file__]}, datasets={})
    for dataset in datasets:
        pool = destination / "cache" / f"{dataset}_baseline_A_pool.npy"
        binding["datasets"][dataset] = dict(pool_sha256=hashlib.sha256(pool.read_bytes()).hexdigest(),
            grid_sha256=hashlib.sha256(old_refs[(dataset, SEEDS[0])][0].tobytes()).hexdigest(),
            baseline_pool_size=prepared["baseline_pool_sizes"][dataset])
    target = destination / "reference_cache_binding.json"
    if target.exists():
        prior = json.loads(target.read_text())
        if prior != binding and not list((destination / "references").glob("*/*/complete.json")):
            # An interrupted first build has no reusable completed reference.
            save_json(target, binding)
        else:
            assert prior == binding, "Reference cache binding changed; refusing stale reuse"
    elif (destination / "references").exists():
        raise AssertionError("Existing references have no binding; refusing unverifiable reuse")
    else:
        save_json(target, binding)


def prepare(out, scans, audit, old_refs, helper, old_mode_valid):
    old_records = {d: read_json(folder / "cache/prepared.json") for d, folder in SOURCES.items()}
    original_plan = read_csv(CORE / "raw_paired_changes.csv")
    baseline_keys = set(zip(original_plan.subject_id, original_plan.baseline_session))
    baseline = scans.loc[[(r.subject_id, r.session) in baseline_keys for r in scans.itertuples()]].sort_values(KEY)
    audit_index = audit.set_index(KEY)
    edges = np.column_stack(np.triu_indices(454, 1))
    preparations, membership, parity = {}, [], []
    for policy in ["registered", "strict"]:
        (out / policy / "cache").mkdir(parents=True, exist_ok=True)
        preparations[policy] = {"tdmi_ranges": {}, "calibration_targets": {}, "baseline_pool_sizes": {}}
    for dataset, folder in SOURCES.items():
        record = old_records[dataset]
        pool_path = track(folder / "cache" / f"{dataset}_baseline_A_pool.npy")
        pool = np.load(pool_path)
        selected = np.load(track(folder / "cache/pool_edge_indices.npy"))
        pieces = {"registered": [], "strict": []}
        offset = 0
        for scan in baseline.loc[baseline.dataset.eq(dataset)].itertuples():
            ts = np.load(source_time_series(scan))
            sigma = tn.pg.lagged_sigma_from_timeseries(ts, 1, edges[selected])
            a, _, _ = tn.pg.var1_from_sigma(sigma)
            reconstructed = a[np.max(abs(np.linalg.eigvals(a)), axis=1) < 1]
            # Slice the original pool after verifying ordering and source values,
            # so filtering does not introduce new floating-point reconstruction drift.
            source_piece = pool[offset:offset + len(reconstructed)]
            np.testing.assert_allclose(reconstructed, source_piece, rtol=0, atol=1e-12)
            error = float(abs(reconstructed - source_piece).max())
            audited = audit_index.loc[(scan.dataset, scan.subject_id, scan.session)]
            assert audited.normalization_baseline_pool_member
            for policy in pieces:
                keep = bool(audited[policy + "_pass"])
                if keep:
                    pieces[policy].append(source_piece)
                membership.append(dict(qc_policy=policy, dataset=dataset, subject_id=scan.subject_id,
                                       session=scan.session, prefix=scan.prefix, retained=keep,
                                       first_pool_row=offset, end_pool_row=offset + len(reconstructed),
                                       pool_rows=len(reconstructed), reconstruction_max_error=error,
                                       retained_mean_fd=audited.retained_mean_fd,
                                       retained_peak_fd=audited.retained_peak_fd,
                                       fail_reasons=audited[policy + "_fail_reasons"]))
            offset += len(reconstructed)
        assert offset == len(pool)
        for policy, selected_pieces in pieces.items():
            fresh = np.concatenate(selected_pieces)
            np.save(out / policy / "cache" / f"{dataset}_baseline_A_pool.npy", fresh)
            for name in ["tdmi_ranges", "calibration_targets"]:
                preparations[policy][name][dataset] = record[name][dataset]
            preparations[policy]["baseline_pool_sizes"][dataset] = len(fresh)
            if policy == "registered" and dataset == "aachen_itbs":
                np.testing.assert_array_equal(fresh, pool)
        # Exact old-pool parity at low/middle/high grid points for both seeds and
        # both numerical definitions, before any new treatment inference.
        for seed in SEEDS:
            grid, tables = old_refs[(dataset, seed)]
            for index in [0, len(grid) // 2, len(grid) - 1]:
                draws, quality = helper.tn.valid_reference(pool, grid[index], helper.COUNT, seed)
                for measure in MEASURES:
                    current = np.sort(draws[measure])
                    expected = tables[measure][index]
                    np.testing.assert_allclose(current, expected, rtol=1e-7, atol=3e-12,
                                               err_msg=f"Old reference parity: {dataset}/{seed}/{index}/{measure}")
                    parity.append(dict(dataset=dataset, seed=seed, grid_index=index, tdmi=grid[index],
                                       measure=measure, max_absolute_difference=float(abs(current - expected).max()),
                                       byte_identical=bool(np.array_equal(current, expected)), passed=True))
    for policy, prepared in preparations.items():
        save_json(out / policy / "cache/prepared.json", prepared)
    pd.DataFrame(membership).to_csv(out / "baseline_membership.csv", index=False)
    pd.DataFrame(parity).to_csv(out / "old_reference_parity.csv", index=False)
    return preparations


def masks(reference):
    edges, core_masks = ws.edge_scopes(reference["workspace"])
    i, j = edges.T
    lab = reference["modules"]
    striatum = reference["labels"].Structure_Name.eq("Striatum").to_numpy()
    def between(a, b):
        return (a[i] & b[j]) | (b[i] & a[j])
    mode_masks = {"whole_brain": core_masks["whole_brain"], "within_workspace": core_masks["within_workspace"],
                  "outside_workspace": core_masks["outside_workspace"], "workspace_boundary": core_masks["between"],
                  "FPN_striatum": between(lab == 5, striatum), "FPN_DAN": between(lab == 5, lab == 2),
                  "FPN_DMN": between(lab == 5, lab == 6)}
    assert {s: int(v.sum()) for s, v in mode_masks.items()} == read_json(rerun.RAW / "scope_pair_counts.json")
    return edges, core_masks, mode_masks


def observed(scan, edges):
    sigma = tn.pg.lagged_sigma_from_timeseries(np.load(source_time_series(scan)), 1, edges)
    atoms = mode_norm.expanded_atoms(sigma)
    core_tdmi = tn.pg.tdmi_from_sigma4(sigma)
    base = mode_norm.raw.CORRECTED if scan.dataset == "aachen_itbs" else ROOT
    for measure in COMPOSITES:
        stored = np.load(track(base / "pipeline/results/phiid_matrices" / f"{scan.prefix}_{measure}.npy"))[edges[:, 0], edges[:, 1]]
        np.testing.assert_allclose(atoms[measure], stored, rtol=1e-9, atol=2e-9)
        atoms[measure] = stored
    return atoms, core_tdmi


def interpolation(out, scans, refs, helper, edges):
    records = []
    rng = np.random.default_rng(2026112619)
    for dataset in sorted(scans.dataset.unique()):
        pool = np.load(out / "cache" / f"{dataset}_baseline_A_pool.npy")
        block = scans.loc[scans.dataset.eq(dataset)]
        queries = []
        for scan in block.iloc[np.linspace(0, len(block) - 1, 8, dtype=int)].itertuples():
            atoms, core_tdmi = observed(scan, edges)
            for index in [atoms["tdmi"].argmin(), atoms["tdmi"].argmax(), *rng.choice(len(edges), 2, replace=False)]:
                queries.append((float(atoms["tdmi"][index]), {m: float(atoms[m][index]) for m in MEASURES}))
        for seed in SEEDS:
            grid, tables = refs[(dataset, seed)]
            for query, (target, values) in enumerate(queries):
                direct, _ = helper.tn.valid_reference(pool, target, helper.COUNT, seed)
                for measure in MEASURES:
                    expected = float(tn.empirical_cdf(np.sort(direct[measure]), values[measure]))
                    actual = float(tn.conditional_cdf(grid, tables[measure], np.array([target]), np.array([values[measure]]))[0])
                    records.append(dict(dataset=dataset, seed=seed, query=query, measure=measure, tdmi=target,
                                        direct=expected, interpolated=actual, absolute_error=abs(expected - actual)))
            print(f"Interpolation validated {out.name}/{dataset}/{seed}", flush=True)
    frame = pd.DataFrame(records)
    frame.to_csv(out / "interpolation_validation.csv", index=False)
    if frame.absolute_error.max() > .002:
        raise AssertionError(f"Interpolation maximum exceeds .002: {frame.absolute_error.max()}")
    return float(frame.absolute_error.max())


def extraction(out, scans, refs, reference, edges, core_masks, mode_masks, acute_keys):
    core_rows, mode_rows, saturation = [], [], []
    for index, scan in enumerate(scans.itertuples()):
        atoms, core_tdmi = observed(scan, edges)
        q = {str(seed): {} for seed in SEEDS}
        for seed in SEEDS:
            grid, tables = refs[(scan.dataset, seed)]
            for measure in MEASURES:
                tdmi = core_tdmi if measure in COMPOSITES else atoms["tdmi"]
                q[str(seed)][measure] = tn.conditional_cdf(grid, tables[measure], tdmi, atoms[measure])
        q["combined"] = {m: (q[str(SEEDS[0])][m] + q[str(SEEDS[1])][m]) / 2 for m in MEASURES}
        metadata = {key: getattr(scan, key) for key in [*KEY, "arm", "mean_fd", "prefix", "n_timepoints", "n_valid_rois", "valid_roi_signature"]}
        for seed, values in q.items():
            matrices = {}
            for measure in COMPOSITES:
                matrix = np.full((454, 454), np.nan)
                matrix[edges[:, 0], edges[:, 1]] = matrix[edges[:, 1], edges[:, 0]] = values[measure]
                matrices[measure] = matrix
            metrics, _ = tms.metric_values(matrices, reference)
            scoped = {scope + "__" + m: float(values[m][mask].mean()) for scope, mask in core_masks.items() for m in COMPOSITES}
            core_rows.append({**metadata, "reference_seed": seed, **metrics, **scoped})
            if (scan.subject_id, scan.session) in acute_keys:
                mode_rows.append({**metadata, "reference_seed": seed,
                                  **{scope + "__" + m: float(values[m][mask].mean()) for scope, mask in mode_masks.items() for m in MODES}})
        saturation.extend({**{key: getattr(scan, key) for key in KEY}, "measure": m,
                           "endpoint_fraction": float(((v == 0) | (v == 1)).mean())} for m, v in q["combined"].items())
        if (index + 1) % 20 == 0 or index + 1 == len(scans):
            print(f"Extracted {out.name}: {index + 1}/{len(scans)} QC-passing scans", flush=True)
    core, modes = pd.DataFrame(core_rows), pd.DataFrame(mode_rows)
    core.to_csv(out / "normalized_scan_metrics_all_seeds.csv", index=False)
    core.loc[core.reference_seed.eq("combined")].to_csv(out / "normalized_scan_metrics.csv", index=False)
    modes.to_csv(out / "normalized_mode_scan_metrics_all_seeds.csv", index=False)
    pd.DataFrame(saturation).to_csv(out / "endpoint_saturation.csv", index=False)
    records = []
    for kind, frame, metrics in [("core", core, [*tms.METRICS, *[s + "__" + m for s in core_masks for m in COMPOSITES]]),
                                  ("modes", modes, rerun.METRICS)]:
        for dataset in sorted(frame.dataset.unique()):
            blocks = [frame.loc[frame.dataset.eq(dataset) & frame.reference_seed.eq(str(seed))].set_index(KEY).sort_index() for seed in SEEDS]
            a, b = blocks
            assert a.index.equals(b.index)
            for metric in metrics:
                difference = a[metric] - b[metric]
                records.append(dict(kind=kind, dataset=dataset, metric=metric, n=len(a),
                                    rms_seed_difference=np.sqrt(np.mean(difference ** 2)),
                                    max_absolute_seed_difference=abs(difference).max(), data_sd=((a[metric] + b[metric]) / 2).std()))
    pd.DataFrame(records).to_csv(out / "monte_carlo_stability.csv", index=False)
    return core, modes


def unchanged_aachen_check(out, core, modes):
    before_core = read_csv(SOURCES["aachen_itbs"] / "scan_metrics_all_seeds_scoped.csv", dtype={"reference_seed": str}).set_index(KEY + ["reference_seed"])
    before_modes = read_csv(MODE_SOURCE / "scan_modes_all_seeds.csv", dtype={"reference_seed": str}).set_index(KEY + ["reference_seed"])
    records = []
    for label, current, previous in [("core", core, before_core), ("modes", modes, before_modes)]:
        selected = current.loc[current.dataset.eq("aachen_itbs")].set_index(KEY + ["reference_seed"])
        columns = [c for c in previous if c in selected and pd.api.types.is_numeric_dtype(previous[c]) and c != "mean_fd"]
        a, b = selected[columns].to_numpy(float), previous.loc[selected.index, columns].to_numpy(float)
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=2e-12)
        records.append(dict(family=label, scans_seeds=len(selected), columns=len(columns), max_absolute_difference=float(abs(a - b).max()), passed=True))
    pd.DataFrame(records).to_csv(out / "unchanged_registered_aachen_check.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    for path in [Path(__file__), Path(mode_norm.__file__), Path(mode_norm.raw.__file__), Path(rerun.__file__),
                 Path(tn.__file__), Path(tn.nn.__file__), Path(tn.pg.__file__), Path(tms.__file__), Path(ws.__file__),
                 Path(mode_norm.legacy.__file__)]:
        track(path)
    receipt = read_json(AUDIT / "audit_receipt.json")
    assert receipt["state"] == "complete" and receipt["all_193_windows_verified"] and not receipt["changed_source_files"]
    for name, digest in receipt["output_sha256"].items():
        path = track(AUDIT / name)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    audit = read_csv(AUDIT / "scan_audit.csv")
    for col in ["registered_pass", "strict_pass", "bounds_verified", "normalization_baseline_pool_member"]:
        audit[col] = rerun.bool_column(audit[col])
    assert len(audit) == 193 and audit.bounds_verified.all() and not audit.duplicated(KEY).any()
    scans = read_csv(CORE / "raw_scan_metrics.csv").sort_values(KEY).reset_index(drop=True)
    audit_index = audit.set_index(KEY)
    assert set(map(tuple, scans[KEY].to_numpy())) == set(audit_index.index)
    scans["mean_fd"] = [audit_index.loc[(r.dataset, r.subject_id, r.session), "retained_mean_fd"] for r in scans.itertuples()]
    helper, old_mode_valid = runtime()
    old_refs = {}
    for dataset in SOURCES:
        old_refs.update(load_previous_refs(dataset))
    save_json(OUT / "analysis_specification.json", dict(created_utc=datetime.now(timezone.utc).isoformat(),
        baseline_eligibility="QC-passing records in original 30 Aachen / 40 SUDMEX baseline pool, independent of post availability.",
        calibration="Same baseline VAR slices/order, seeds 2026112511/12, 32768 draws per .005-nat target, unchanged TDMI grids.",
        numerical_rules="Five composites retain original null-TDMI rejection. Six modes retain original extra atom-TDMI rejection and 1e-12 tie canonicalization. Shared initial draws, separate accepted subsets, original-family replenishment fallback.",
        registered_aachen="All 30 baseline scans pass: reuse original validated five- and six-measure references exactly.",
        other_pools="Rebuild registered SUDMEX, strict SUDMEX and strict Aachen from QC-clean baseline pool slices.",
        output="Combined and separate seeds for 40 core measures on all passing scans; 42 mode measures on passing scans from original acute analysis.",
        interpretation="Numerical seed repeats are not independent replications; conditional percentiles are not additive nats. Baseline prior uncertainty remains unpropagated."))
    prepared = prepare(OUT, scans, audit, old_refs, helper, old_mode_valid)
    save_json(OUT / "source_sha256.json", HASHES)
    print("Old-pool reconstruction and all sampled reference parity checks passed", flush=True)
    if args.prepare_only:
        save_json(OUT / "status.json", {"state": "prepared", "elapsed_seconds": time.perf_counter() - started})
        return
    reference = tms.load_reference(ROOT)
    for path in reference["input_files"]:
        track(path)
    edges, core_masks, mode_masks = masks(reference)
    plan = read_csv(CORE / "raw_paired_changes.csv")
    acute_keys = set(zip(plan.subject_id, plan.baseline_session)) | set(zip(plan.subject_id, plan.followup_session))
    verification = []
    for policy in ["registered", "strict"]:
        destination = OUT / policy
        pass_flag = audit.set_index(KEY)[policy + "_pass"]
        selected = scans.loc[[bool(pass_flag.loc[(r.dataset, r.subject_id, r.session)]) for r in scans.itertuples()]].copy()
        helper.tms = types.SimpleNamespace(DATASETS=("sudmex_tms",) if policy == "registered" else tuple(SOURCES))
        bind_reference_cache(destination, prepared[policy], helper.tms.DATASETS, helper, old_refs)
        refs = helper.build_references(destination, prepared[policy])
        if policy == "registered":
            refs.update({key: value for key, value in old_refs.items() if key[0] == "aachen_itbs"})
        for key, (grid, _) in refs.items():
            np.testing.assert_array_equal(grid, old_refs[key][0])
        helper.tms = types.SimpleNamespace(DATASETS=tuple(SOURCES))
        max_error = interpolation(destination, selected, refs, helper, edges)
        calibration = helper.calibrate_finite_records(destination, prepared[policy], refs)
        save_json(destination / "validation_before_effects.json", {"max_interpolation_error": max_error, **calibration})
        core, modes = extraction(destination, selected, refs, reference, edges, core_masks, mode_masks, acute_keys)
        if policy == "registered":
            unchanged_aachen_check(destination, core, modes)
        whitelists = read_json(AUDIT / "cohort_subjects.json")
        participant_ids = {subject for ids in whitelists[policy].values() for subject in ids}
        analysis_plan = plan.loc[plan.subject_id.isin(participant_ids)].copy()
        eligible_modes = modes.loc[modes.subject_id.isin(participant_ids)].copy()
        mode_out = rerun.DEFAULT_OUT / policy / "normalized_recalibrated"
        treatment, conjunction, _ = rerun.normalized_inference(eligible_modes, analysis_plan, mode_out)
        save_json(mode_out / "verification.json", dict(passed=True, state="complete", qc_policy=policy,
            normalized_scan_source=str((destination / "normalized_mode_scan_metrics_all_seeds.csv").relative_to(ROOT)),
            normalized_scan_sha256=hashlib.sha256((destination / "normalized_mode_scan_metrics_all_seeds.csv").read_bytes()).hexdigest(),
            normalized_reference="FD-clean recalibrated baseline prior",
            shared_holm6_results=int(((conjunction.reference_seed == "combined") & (conjunction.p_holm6 < .05)).sum()),
            shared_holm42_results=int(((conjunction.reference_seed == "combined") & (conjunction.p_holm42 < .05)).sum())))
        record = dict(qc_policy=policy, passed=True, state="complete", passing_scans=len(selected), core_seed_rows=len(core),
                      mode_seed_rows=len(modes), max_interpolation_error=max_error, calibration=calibration,
                      baseline_pool_sizes=prepared[policy]["baseline_pool_sizes"],
                      normalized_reference="FD-clean baseline prior; registered Aachen reference unchanged")
        save_json(destination / "verification.json", record)
        verification.append(record)
        print(f"Complete fresh calibration and metrics for {policy}", flush=True)
    save_json(OUT / "source_sha256.json", HASHES)
    for relative, digest in HASHES.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest, relative
    output_hashes = {str(p.relative_to(OUT)): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in OUT.rglob("*") if p.is_file() and p.name not in ["verification.json", "output_sha256.json"]}
    save_json(OUT / "output_sha256.json", output_hashes)
    save_json(OUT / "verification.json", {"passed": True, "state": "complete", "policies": verification,
              "sources_unchanged": True, "source_files_hashed": len(HASHES), "elapsed_seconds": time.perf_counter() - started})
    print("COMPLETE", OUT, flush=True)


if __name__ == "__main__":
    with threadpool_limits(limits=2):
        main()
