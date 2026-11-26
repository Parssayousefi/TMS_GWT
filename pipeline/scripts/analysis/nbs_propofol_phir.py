"""
Network-Based Statistic (NBS) on per-scan Phi_R matrices for the Cambridge
propofol dataset, replicating Luppi et al. 2024.

Three analyses:
  - Figure 4B: deep vs awake (paired)             (decreases / increases)
  - Figure 4C: deep vs recovery (paired)          (decreases / increases)
  - Figure 4D: conjunction of decreases-during-anaesthesia AND
               restored-on-recovery: edges with t > 3 in BOTH
               (awake - deep) AND (recovery - deep) within subjects who
               have all three conditions; permutation null = sign-flip both
               contrasts independently per subject.

Vectorised paired-NBS implementation: ~25 perms/sec on 454-ROI matrices.

Method (Zalesky et al. 2010; Luppi 2024 Methods):
  - primary threshold t = 3   (~= F = 9 in Luppi convention; one-tail)
  - 5000 permutations
  - alpha = 0.05 cluster-corrected (max-component-size null)
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from analysis_loader import load_group_matrices

DATASET = 'cambridge_propofol'
T_THRESH = 3.0
N_PERMUTATIONS = 5000
SEED = 42


def symmetrise(mat: np.ndarray) -> np.ndarray:
    out = np.where(np.isnan(mat), 0.0, mat)
    out = 0.5 * (out + out.T)
    np.fill_diagonal(out, 0.0)
    return out


def edge_stack(matrices: np.ndarray, iu: tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    """(n_subj, N, N) -> (n_subj, n_edges) of upper-triangle edge values."""
    out = np.stack([symmetrise(m) for m in matrices], axis=0)
    return out[:, iu[0], iu[1]]


def paired_t(diffs: np.ndarray) -> np.ndarray:
    n = diffs.shape[0]
    mean = diffs.mean(axis=0)
    sd = diffs.std(axis=0, ddof=1)
    sd_safe = np.where(sd > 0, sd, np.inf)
    return mean / (sd_safe / np.sqrt(n))


def components_from_edge_mask(supra: np.ndarray, ij: np.ndarray, n_rois: int):
    """Return (component_node_label_per_edge, edge_counts_per_component_value)."""
    if not supra.any():
        return np.full(supra.size, -1, dtype=np.int32), np.zeros(0, dtype=np.int32)
    sel = np.where(supra)[0]
    rows = ij[0, sel]
    cols = ij[1, sel]
    data = np.ones(sel.size, dtype=np.int8)
    adj = csr_matrix(
        (np.concatenate([data, data]),
         (np.concatenate([rows, cols]), np.concatenate([cols, rows]))),
        shape=(n_rois, n_rois),
    )
    _, node_labels = connected_components(adj, directed=False)
    edge_label = np.full(supra.size, -1, dtype=np.int32)
    edge_label[sel] = node_labels[rows]
    label_vals, counts = np.unique(node_labels[rows], return_counts=True)
    return edge_label, dict(zip(label_vals.tolist(), counts.tolist()))


def summarise_components(edge_label: np.ndarray, ij: np.ndarray) -> list[dict]:
    """Build [{label, edges, nodes}] sorted by edge count descending."""
    sel = np.where(edge_label >= 0)[0]
    if sel.size == 0:
        return []
    label_vals, counts = np.unique(edge_label[sel], return_counts=True)
    out = []
    for lab, cnt in zip(label_vals.tolist(), counts.tolist()):
        edges_in_c = sel[edge_label[sel] == lab]
        nodes_in_c = np.unique(np.concatenate([ij[0, edges_in_c], ij[1, edges_in_c]]))
        out.append({'label': int(lab), 'edges': int(cnt), 'nodes': int(nodes_in_c.size)})
    return sorted(out, key=lambda d: -d['edges'])


def run_paired_nbs(diff_edges: np.ndarray, ij: np.ndarray, n_rois: int,
                   thresh: float, k: int, rng: np.random.Generator,
                   tag: str) -> dict:
    """One-tailed paired NBS; tests positive paired-t (mean(diff) > 0)."""
    t0 = time.time()
    n_subj = diff_edges.shape[0]

    obs_t = paired_t(diff_edges)
    obs_supra = obs_t > thresh
    obs_label, _ = components_from_edge_mask(obs_supra, ij, n_rois)
    obs_components = summarise_components(obs_label, ij)
    obs_max = obs_components[0]['edges'] if obs_components else 0
    print(f"  [{tag}] observed: suprathreshold edges={int(obs_supra.sum())} "
          f"components={len(obs_components)} max-component={obs_max} edges",
          flush=True)

    null_max = np.zeros(k, dtype=np.int32)
    flips = rng.choice([-1.0, 1.0], size=(k, n_subj))
    chunk = 100
    for start in range(0, k, chunk):
        end = min(start + chunk, k)
        for p in range(start, end):
            t_perm = paired_t(diff_edges * flips[p][:, None])
            supra = t_perm > thresh
            _, counts = components_from_edge_mask(supra, ij, n_rois)
            null_max[p] = max(counts.values()) if counts else 0
        elapsed = time.time() - t0
        rate = end / elapsed if elapsed else 0
        eta = (k - end) / rate if rate else 0
        print(f"  [{tag}] perm {end}/{k}  null_max so far={null_max[:end].max() if end else 0}  "
              f"rate={rate:.1f}/s eta={eta:.0f}s", flush=True)

    comp_results = []
    for c in obs_components:
        p_corr = float((null_max >= c['edges']).sum() + 1) / (k + 1)
        comp_results.append({
            'edges': c['edges'], 'nodes': c['nodes'], 'p_corrected': p_corr,
        })
    return {'tag': tag, 'n_components': len(obs_components),
            'components': comp_results, 'null_max': null_max}


def run_conjunction_nbs(diff_a: np.ndarray, diff_b: np.ndarray,
                        ij: np.ndarray, n_rois: int,
                        thresh: float, k: int, rng: np.random.Generator,
                        tag: str) -> dict:
    """Conjunction NBS: edges with t>thresh in both contrasts.
    diff_a, diff_b: (n_subj, n_edges); same subject ordering across both.
    Permutation null: independent sign-flips on diff_a and diff_b per subject.
    """
    t0 = time.time()
    n_subj = diff_a.shape[0]

    t_a = paired_t(diff_a)
    t_b = paired_t(diff_b)
    obs_supra = (t_a > thresh) & (t_b > thresh)
    obs_label, _ = components_from_edge_mask(obs_supra, ij, n_rois)
    obs_components = summarise_components(obs_label, ij)
    obs_max = obs_components[0]['edges'] if obs_components else 0
    print(f"  [{tag}] observed: edges_a={int((t_a>thresh).sum())} "
          f"edges_b={int((t_b>thresh).sum())} conjunction={int(obs_supra.sum())} "
          f"max-component={obs_max} edges", flush=True)

    null_max = np.zeros(k, dtype=np.int32)
    flips_a = rng.choice([-1.0, 1.0], size=(k, n_subj))
    flips_b = rng.choice([-1.0, 1.0], size=(k, n_subj))
    chunk = 100
    for start in range(0, k, chunk):
        end = min(start + chunk, k)
        for p in range(start, end):
            ta = paired_t(diff_a * flips_a[p][:, None])
            tb = paired_t(diff_b * flips_b[p][:, None])
            supra = (ta > thresh) & (tb > thresh)
            _, counts = components_from_edge_mask(supra, ij, n_rois)
            null_max[p] = max(counts.values()) if counts else 0
        elapsed = time.time() - t0
        rate = end / elapsed if elapsed else 0
        eta = (k - end) / rate if rate else 0
        print(f"  [{tag}] perm {end}/{k}  null_max so far={null_max[:end].max() if end else 0}  "
              f"rate={rate:.1f}/s eta={eta:.0f}s", flush=True)

    comp_results = []
    for c in obs_components:
        p_corr = float((null_max >= c['edges']).sum() + 1) / (k + 1)
        comp_results.append({
            'edges': c['edges'], 'nodes': c['nodes'], 'p_corrected': p_corr,
        })
    return {'tag': tag, 'n_components': len(obs_components),
            'components': comp_results, 'null_max': null_max,
            'edges_a_supra': int((t_a > thresh).sum()),
            'edges_b_supra': int((t_b > thresh).sum()),
            'conjunction_edges': int(obs_supra.sum())}


def load_paired(awake_subj_to_idx: dict, mats_a: np.ndarray, subjects_a: list,
                mats_b: np.ndarray, subjects_b: list,
                iu: tuple[np.ndarray, np.ndarray]) -> tuple:
    paired = sorted(set(subjects_a) & set(subjects_b))
    a_idx = [subjects_a.index(s) for s in paired]
    b_idx = [subjects_b.index(s) for s in paired]
    a_edges = edge_stack(mats_a[a_idx], iu)
    b_edges = edge_stack(mats_b[b_idx], iu)
    return paired, a_edges, b_edges


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--measure', default='phi_r',
                        choices=['phi_r', 'phi', 'synergy', 'redundancy'])
    parser.add_argument('--dataset', default=DATASET,
                        help='dataset id (default: %(default)s)')
    args = parser.parse_args()
    measure = args.measure
    dataset = args.dataset

    project_dir = Path(__file__).resolve().parent.parent.parent
    out_dir = project_dir / 'results' / 'group_analysis' / dataset / 'nbs' / measure
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading per-scan {measure} matrices for {dataset}", flush=True)
    awake_mats, awake_subj = load_group_matrices(dataset, 'restawake', measure)
    deep_mats, deep_subj = load_group_matrices(dataset, 'restdeep', measure)
    rec_mats, rec_subj = load_group_matrices(dataset, 'restrecovery', measure)
    print(f"  restawake N={len(awake_subj)}", flush=True)
    print(f"  restdeep  N={len(deep_subj)}", flush=True)
    print(f"  restrecovery N={len(rec_subj)}", flush=True)

    n_rois = awake_mats.shape[-1]
    iu = np.triu_indices(n_rois, k=1)
    ij = np.array(iu)
    print(f"  matrices: {n_rois}×{n_rois}, edges in upper triangle: {ij.shape[1]}",
          flush=True)

    rng = np.random.default_rng(SEED)
    all_results: list[dict] = []

    def save_partial():
        for r in all_results:
            np.save(out_dir / f'nbs_{r["tag"]}_null_max.npy',
                    np.asarray(r['null_max'], dtype=np.int32))
        meta_partial = {
            'dataset': dataset, 'measure': measure,
            't_thresh': T_THRESH, 'k_permutations': N_PERMUTATIONS, 'seed': SEED,
            'tails': [{k: v for k, v in r.items() if k != 'null_max'}
                      for r in all_results],
        }
        with open(out_dir / f'nbs_{measure}_propofol_partial.yaml', 'w') as fh:
            yaml.safe_dump(meta_partial, fh, sort_keys=False)

    # ---- Figure 4B: deep vs awake (paired) ----
    paired_ad, awake_e_ad, deep_e_ad = load_paired({}, awake_mats, awake_subj,
                                                    deep_mats, deep_subj, iu)
    print(f"\n{'='*70}\nFig 4B: NBS deep vs awake (paired N={len(paired_ad)})\n{'='*70}",
          flush=True)
    diff_dec = awake_e_ad - deep_e_ad
    diff_inc = deep_e_ad - awake_e_ad
    res_4b_dec = run_paired_nbs(diff_dec, ij, n_rois, T_THRESH, N_PERMUTATIONS,
                                rng, 'fig4B_decreases')
    all_results.append(res_4b_dec)
    save_partial()
    res_4b_inc = run_paired_nbs(diff_inc, ij, n_rois, T_THRESH, N_PERMUTATIONS,
                                rng, 'fig4B_increases')
    all_results.append(res_4b_inc)
    save_partial()

    # ---- Figure 4C: deep vs recovery (paired) ----
    paired_dr, deep_e_dr, rec_e_dr = load_paired({}, deep_mats, deep_subj,
                                                  rec_mats, rec_subj, iu)
    print(f"\n{'='*70}\nFig 4C: NBS deep vs recovery (paired N={len(paired_dr)})\n{'='*70}",
          flush=True)
    diff_dec_dr = rec_e_dr - deep_e_dr
    diff_inc_dr = deep_e_dr - rec_e_dr
    res_4c_dec = run_paired_nbs(diff_dec_dr, ij, n_rois, T_THRESH, N_PERMUTATIONS,
                                rng, 'fig4C_decreases')
    all_results.append(res_4c_dec)
    save_partial()
    res_4c_inc = run_paired_nbs(diff_inc_dr, ij, n_rois, T_THRESH, N_PERMUTATIONS,
                                rng, 'fig4C_increases')
    all_results.append(res_4c_inc)
    save_partial()

    # ---- Figure 4D: conjunction (subjects with awake AND deep AND recovery) ----
    paired_3 = sorted(set(awake_subj) & set(deep_subj) & set(rec_subj))
    print(f"\n{'='*70}\nFig 4D: conjunction NBS (awake-deep) INTERSECT (rec-deep), "
          f"paired N={len(paired_3)}\n{'='*70}", flush=True)
    a_idx = [awake_subj.index(s) for s in paired_3]
    d_idx = [deep_subj.index(s) for s in paired_3]
    r_idx = [rec_subj.index(s) for s in paired_3]
    a_e = edge_stack(awake_mats[a_idx], iu)
    d_e = edge_stack(deep_mats[d_idx], iu)
    r_e = edge_stack(rec_mats[r_idx], iu)
    diff_ad = a_e - d_e   # awake > deep
    diff_rd = r_e - d_e   # recovery > deep
    res_4d = run_conjunction_nbs(diff_ad, diff_rd, ij, n_rois, T_THRESH,
                                  N_PERMUTATIONS, rng, 'fig4D_decreases_conj')
    all_results.append(res_4d)
    save_partial()

    # ---- Report ----
    print(f"\n{'='*70}\nRESULTS\n{'='*70}", flush=True)
    for r in all_results:
        head = f"\n[{r['tag']}] components={r['n_components']}"
        if 'conjunction_edges' in r:
            head += (f"  edges_a_supra={r['edges_a_supra']} "
                     f"edges_b_supra={r['edges_b_supra']} "
                     f"conjunction={r['conjunction_edges']}")
        print(head, flush=True)
        for c in r['components'][:5]:
            star = '*' if c['p_corrected'] < 0.05 else ' '
            print(f"   [{star}] edges={c['edges']:5d} nodes={c['nodes']:3d} "
                  f"p_corr={c['p_corrected']:.4f}", flush=True)

    meta = {
        'dataset': dataset, 'measure': measure,
        't_thresh': T_THRESH, 'k_permutations': N_PERMUTATIONS, 'seed': SEED,
        'method': 'NBS (Zalesky 2010); paired t; sign-flip null; conjunction = '
                  'edges suprathreshold in BOTH contrasts independently. '
                  't>3 ~= F=9 (Luppi 2024).',
        'paired_subjects': {
            'fig4B_awake_deep': paired_ad,
            'fig4C_deep_recovery': paired_dr,
            'fig4D_conjunction': paired_3,
        },
        'tails': [{k: v for k, v in r.items() if k != 'null_max'}
                  for r in all_results],
    }
    yaml_path = out_dir / f'nbs_{measure}_propofol_full.yaml'
    with open(yaml_path, 'w') as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    print(f"\nSaved metadata: {yaml_path}", flush=True)

    for r in all_results:
        np.save(out_dir / f'nbs_{r["tag"]}_null_max.npy',
                np.asarray(r['null_max'], dtype=np.int32))


if __name__ == '__main__':
    main()
