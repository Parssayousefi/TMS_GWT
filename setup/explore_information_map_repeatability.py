"""Bounded exploratory alternative to the haemodynamic publication narrative.

Read corrected TMS matrices only. Derive the production synergy-minus-redundancy
rank gradient independently per scan, without fitting a new reference atlas.
This is a new exploratory analysis, not a preregistered treatment endpoint.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util
import json

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "output/corrected_tms_analysis/20261126T112054_284547Z"
CORRECTED = ROOT / "reruns/aachen_skullstrip_skip_20261125"
REFERENCE = ROOT / "pipeline/results/group_analysis/hcp/rest"
SEED = 20261126141
N_PERM = 19999


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, default=str) + "\n", encoding="utf-8")


def corr_matrix(a, b):
    ar = stats.rankdata(a, axis=1)
    br = stats.rankdata(b, axis=1)
    ar -= ar.mean(1, keepdims=True)
    br -= br.mean(1, keepdims=True)
    ar /= np.linalg.norm(ar, axis=1, keepdims=True)
    br /= np.linalg.norm(br, axis=1, keepdims=True)
    return ar @ br.T


def identify_stats(matrix, rng, test=True):
    n = len(matrix)
    diagonal = np.diag(matrix)
    unmatched = (matrix.sum() - diagonal.sum()) / (n * (n - 1))
    advantage = diagonal.mean() - unmatched
    # Strict winners: ties do not count as successful identification.
    forward = diagonal > np.where(np.eye(n, dtype=bool), -np.inf, matrix).max(0)
    reverse = diagonal > np.where(np.eye(n, dtype=bool), -np.inf, matrix).max(1)
    p = np.nan
    if test:
        # Shuffle identities, preserving each complete map and all spatial dependence.
        indices = np.arange(n)
        perm = np.stack([rng.permutation(n) for _ in range(N_PERM)])
        perm_mean = matrix[indices[None, :], perm].mean(1)
        p = (1 + np.count_nonzero(perm_mean >= diagonal.mean() - 1e-12)) / (N_PERM + 1)
    boot = diagonal[rng.integers(n, size=(4000, n))].mean(1)
    return dict(n=n, matched_rho_mean=diagonal.mean(), unmatched_rho_mean=unmatched,
                identification_advantage=advantage,
                matched_mean_ci_low=np.quantile(boot, .025),
                matched_mean_ci_high=np.quantile(boot, .975),
                forward_correct=int(forward.sum()), reverse_correct=int(reverse.sum()),
                forward_accuracy=forward.mean(), reverse_accuracy=reverse.mean(),
                chance_accuracy=1/n, p=p)


def main():
    out = ROOT / "output/information_map_repeatability" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out.mkdir(parents=True)
    print("OUTPUT:", out, flush=True)
    specification = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "New exploratory narrative viability analysis after review of earlier results; not preregistered.",
        "feature": "454-region production node-strength synergy rank minus redundancy rank, separately estimated per scan.",
        "similarity": "Spearman correlation of the complete continuous regional gradients; no feature selection or training.",
        "primary_groups": {"aachen": "15 people, same-visit pre/post similarity matrices averaged over 90% and 120% visits before identity permutation", "sudmex_active": "22 people, t0 to t1", "sudmex_sham": "18 people, t0 to t1"},
        "family": "Holm4 across three one-sided identity-permutation tests and one two-sided SUDMEX active-minus-sham retention contrast.",
        "permutation_statistic": "Mean matched correlation minus mean unmatched correlation. Shuffle participant identities, jointly across both Aachen visits.",
        "permutations": N_PERM, "seed": SEED,
        "arm_contrast": "Fisher-z(pre/post gradient Spearman rho) ~ active + pre FD + post FD; OLS HC3, t reference.",
        "secondary_descriptive": "Individual Aachen visits; Aachen baseline-to-baseline similarity; cortical-only profiles; HCP map correspondence. No spatial-parcel p-values.",
        "limits": "Identity retention does not establish absence of stimulation effects, neural specificity, equivalence, a biomarker, or superiority to conventional FC. Repeatable anatomy, motion and preprocessing may contribute. Fixed HCP mask does not construct participant maps.",
        "routing": {"aachen": str(CORRECTED), "sudmex": str(ROOT / 'pipeline/results')},
    }
    # Specification is written before reading participant matrices.
    save_json(out / "analysis_specification.json", specification)
    paths = [Path(__file__), CORE / "raw_scan_metrics.csv", CORE / "raw_paired_changes.csv",
             REFERENCE / "hcp_rest_workspace_mask.npy", REFERENCE / "hcp_rest_synergy_redundancy_gradient.npy"]
    script = ROOT / "pipeline/scripts/analysis/05_group_average_workspace.py"
    paths.append(script)
    loader = importlib.util.spec_from_file_location("production_workspace", script)
    production = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(production)
    scans = pd.read_csv(CORE / "raw_scan_metrics.csv")
    pairs = pd.read_csv(CORE / "raw_paired_changes.csv")
    assert pairs.groupby("dataset").size().to_dict() == {"aachen_itbs": 30, "sudmex_tms": 40}
    keys = {(r.subject_id, r.baseline_session) for r in pairs.itertuples()} | {(r.subject_id, r.followup_session) for r in pairs.itertuples()}
    selected = scans.loc[[(r.subject_id, r.session) in keys for r in scans.itertuples()]].copy()
    assert len(selected) == len(keys) == 140
    assert not selected.duplicated(["subject_id", "session"]).any()
    profiles, strengths = {}, {}
    source_hashes = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    hcp_mask = np.load(REFERENCE / "hcp_rest_workspace_mask.npy").astype(bool)
    hcp_gradient = np.load(REFERENCE / "hcp_rest_synergy_redundancy_gradient.npy")
    assert hcp_mask.sum() == 235
    obs = []
    offdiag = ~np.eye(454, dtype=bool)
    for r in selected.itertuples():
        base = CORRECTED if r.dataset == "aachen_itbs" else ROOT
        ss = []
        for measure in ["synergy", "redundancy"]:
            path = base / "pipeline/results/phiid_matrices" / f"{r.prefix}_{measure}.npy"
            source_hashes[str(path.relative_to(ROOT))] = sha(path)
            matrix = np.load(path)
            assert matrix.shape == (454, 454) and np.isfinite(matrix[offdiag]).all()
            strength = production.compute_node_strength(matrix)
            np.testing.assert_allclose(strength, np.nansum(matrix, axis=1))
            ss.append(strength)
        mask, gradient, _, _ = production.identify_workspace(*ss)
        independent = pd.Series(ss[0]).rank().to_numpy() - pd.Series(ss[1]).rank().to_numpy()
        np.testing.assert_array_equal(gradient, independent)
        assert np.isfinite(gradient).all() and np.std(gradient) > 0
        key = (r.subject_id, r.session)
        profiles[key], strengths[key] = gradient, np.asarray(ss)
        obs.append(dict(dataset=r.dataset, subject_id=r.subject_id, session=r.session, arm=r.arm,
                        mean_fd=r.mean_fd, hcp_rho=stats.spearmanr(gradient, hcp_gradient).statistic,
                        hcp_dice=production.dice_coefficient(mask, hcp_mask), workspace_count=mask.sum()))
    observations = pd.DataFrame(obs)
    observations.to_csv(out / "scan_profiles_summary.csv", index=False)
    selected = selected.set_index(["subject_id", "session"])

    rng = np.random.default_rng(SEED)
    groups = {
        "aachen": (pairs.query("dataset == 'aachen_itbs'"), [("pre90", "post90"), ("pre120", "post120")]),
        "sudmex_active": (pairs.query("dataset == 'sudmex_tms' and arm == 'active'"), [("t0", "t1")]),
        "sudmex_sham": (pairs.query("dataset == 'sudmex_tms' and arm == 'sham'"), [("t0", "t1")]),
    }
    primary, descriptive, participants, group_maps, similarities = [], [], [], [], {}
    for name, (group, sessions) in groups.items():
        ids = sorted(group.subject_id.unique())
        matrices, cortical = [], []
        for pre, post in sessions:
            a = np.stack([profiles[(s, pre)] for s in ids])
            b = np.stack([profiles[(s, post)] for s in ids])
            m = corr_matrix(a, b)
            # Independently check every diagonal and one off-diagonal per row.
            for i in range(len(ids)):
                for j in [i, (i+1) % len(ids)]:
                    np.testing.assert_allclose(m[i,j], stats.spearmanr(a[i], b[j]).statistic, atol=1e-12)
            matrices.append(m)
            cortical.append(corr_matrix(a[:, :400], b[:, :400]))
            descriptive.append(dict(group=name, comparison=f"{post}_vs_{pre}", **identify_stats(m, rng, test=False)))
        m = np.mean(matrices, axis=0)
        similarities[name] = m
        np.save(out / f"{name}_similarity.npy", m)
        primary.append(dict(group=name, endpoint="identity_retention", **identify_stats(m, rng)))
        descriptive.append(dict(group=name, comparison="cortical_only", **identify_stats(np.mean(cortical, axis=0), rng, test=False)))
        for i, s in enumerate(ids):
            participants.append(dict(group=name, subject_id=s, matched_rho=m[i,i],
                                     mean_other_rho=np.delete(m[i], i).mean(),
                                     pre_fd=np.mean([selected.loc[(s, pre), "mean_fd"] for pre, _ in sessions]),
                                     post_fd=np.mean([selected.loc[(s, post), "mean_fd"] for _, post in sessions])))
        map_values = {}
        for phase, index in [("pre", 0), ("post", 1)]:
            st = np.stack([np.mean([strengths[(s, visit[index])] for visit in sessions], axis=0) for s in ids]).mean(0)
            mask, grad, _, _ = production.identify_workspace(*st)
            map_values[phase] = (mask, grad)
            group_maps.append(dict(group=name, phase=phase, n=len(ids), workspace_count=int(mask.sum()),
                                   hcp_dice=production.dice_coefficient(mask, hcp_mask),
                                   hcp_rho=stats.spearmanr(grad, hcp_gradient).statistic))
        group_maps.append(dict(group=name, phase="pre_post_similarity", n=len(ids),
                               hcp_dice=np.nan, hcp_rho=np.nan,
                               pre_post_dice=production.dice_coefficient(map_values["pre"][0], map_values["post"][0]),
                               pre_post_rho=stats.spearmanr(map_values["pre"][1], map_values["post"][1]).statistic))
    ids = sorted(groups["aachen"][0].subject_id.unique())
    aa_control = corr_matrix(np.stack([profiles[(s, "pre90")] for s in ids]),
                             np.stack([profiles[(s, "pre120")] for s in ids]))
    descriptive.append(dict(group="aachen", comparison="baseline_between_visits", **identify_stats(aa_control, rng, test=False)))
    pp = pd.DataFrame(participants)
    pp.to_csv(out / "participant_similarity.csv", index=False)
    sud = pp.loc[pp.group.str.startswith("sudmex")].copy()
    sud["active"] = sud.group.eq("sudmex_active").astype(float)
    x = sm.add_constant(sud[["active", "pre_fd", "post_fd"]])
    y = np.arctanh(sud.matched_rho)
    model = sm.OLS(y, x).fit(cov_type="HC3", use_t=True)
    # Independent HC3 sandwich and coefficient check.
    xx, yy = x.to_numpy(), y.to_numpy()
    pinv = np.linalg.pinv(xx)
    beta = pinv @ yy
    residual = yy - xx @ beta
    leverage = np.sum(xx * pinv.T, axis=1)
    covariance = (pinv * (residual / (1-leverage))**2) @ pinv.T
    np.testing.assert_allclose(beta, model.params, atol=1e-12)
    np.testing.assert_allclose(covariance, model.cov_params(), atol=1e-12)
    ci = model.conf_int().loc["active"]
    primary.append(dict(group="sudmex_active_minus_sham", endpoint="retention_arm_contrast",
                        n=len(sud), estimate=model.params["active"], ci_low=ci.iloc[0], ci_high=ci.iloc[1],
                        p=model.pvalues["active"], units="Fisher-z of matched gradient correlation"))
    primary = pd.DataFrame(primary)
    primary["p_holm4"] = multipletests(primary.p, method="holm")[1]
    primary.to_csv(out / "primary_results.csv", index=False)
    pd.DataFrame(descriptive).to_csv(out / "descriptive_sensitivities.csv", index=False)
    pd.DataFrame(group_maps).to_csv(out / "group_map_correspondence.csv", index=False)

    # Reproduce archived group QC before accepting the independently derived maps.
    for name, archive, condition, phase in [
        ("aachen", "aachen_corrected_20261125T234156_788094Z", "baseline_mean", "pre"),
        ("sudmex_sham", "20261125T201530_636948Z", "sham_t0", "pre"),
        ("sudmex_sham", "20261125T201530_636948Z", "sham_t1", "post"),
    ]:
        archived = json.loads((ROOT / "output/control_workspace_qc" / archive / "summary.json").read_text())["conditions"][condition]
        current = next(row for row in group_maps if row["group"] == name and row["phase"] == phase)
        np.testing.assert_allclose(current["hcp_dice"], archived["dice_hcp"], atol=1e-12)
        np.testing.assert_allclose(current["hcp_rho"], archived["gradient_rho_hcp"], atol=1e-12)
    for path, digest in source_hashes.items():
        assert sha(ROOT / path) == digest, path
    save_json(out / "source_sha256.json", source_hashes)
    save_json(out / "verification.json", {"passed": True, "n_scans": 140, "n_matrices": 280,
              "all_gradients_match_independent_ranks": True, "similarities_match_scipy": True,
              "three_archived_group_maps_reproduced": True, "hc3_matches_independent_sandwich": True,
              "inputs_unchanged": True})
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), constrained_layout=True)
    titles = ["Aachen · two visits", "SUDMEX · active", "SUDMEX · sham"]
    for ax, (name, m), title in zip(axes, similarities.items(), titles):
        im = ax.imshow(m, vmin=-.2, vmax=1, cmap="viridis")
        ax.set(title=title, xlabel="Post-stimulation participant", ylabel="Baseline participant")
        ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=axes, label="Regional gradient Spearman correlation", shrink=.75)
    fig.savefig(out / "information_map_repeatability.png", dpi=160)
    plt.close(fig)
    (out.parent / "LATEST.txt").write_text(out.name + "\n", encoding="utf-8")
    print(primary.to_string(index=False), flush=True)
    print("DESCRIPTIVE", pd.DataFrame(descriptive).to_string(index=False), flush=True)
    print("MAPS", pd.DataFrame(group_maps).to_string(index=False), flush=True)
    print("Verification passed", flush=True)


if __name__ == "__main__":
    main()
