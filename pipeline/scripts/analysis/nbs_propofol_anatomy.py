"""
Anatomical interpretation of the conjunction-NBS Phi_R components for
Cambridge propofol (Luppi 2024 Figure 4D replication).

For each significant conjunction component, list its ROIs with:
  - Yeo network membership
  - Structure (Schaefer cortical label or Tian subcortical)
  - HCP synergistic workspace membership

Then run an enrichment test: vs a uniform random subset of the same size
drawn from all 454 ROIs, are the conjunction ROIs over-represented in
- the HCP synergistic workspace
- the Default Mode Network
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from analysis_loader import load_group_matrices
from load_roi_mapping import load_roi_network_mapping

DATASET = 'cambridge_propofol'
T_THRESH = 3.0
N_BOOT = 10_000
SEED = 42


def symmetrise(mat: np.ndarray) -> np.ndarray:
    out = np.where(np.isnan(mat), 0.0, mat)
    out = 0.5 * (out + out.T)
    np.fill_diagonal(out, 0.0)
    return out


def edge_stack(matrices, iu):
    s = np.stack([symmetrise(m) for m in matrices], axis=0)
    return s[:, iu[0], iu[1]]


def paired_t(diffs):
    n = diffs.shape[0]
    mean = diffs.mean(axis=0)
    sd = diffs.std(axis=0, ddof=1)
    sd_safe = np.where(sd > 0, sd, np.inf)
    return mean / (sd_safe / np.sqrt(n))


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

    rsn_labels, rsn_names, roi_df = load_roi_network_mapping(project_dir)
    roi_df = roi_df.reset_index(drop=True)

    workspace_mask_path = (project_dir / 'results' / 'group_analysis'
                           / 'hcp' / 'rest' / 'hcp_rest_workspace_mask.npy')
    workspace = np.load(workspace_mask_path).astype(bool)
    print(f"HCP workspace ROIs: {int(workspace.sum())}/{workspace.size}", flush=True)

    print(f"Loading {measure} matrices for {dataset}", flush=True)
    awake_mats, awake_subj = load_group_matrices(dataset, 'restawake', measure)
    deep_mats, deep_subj = load_group_matrices(dataset, 'restdeep', measure)
    rec_mats, rec_subj = load_group_matrices(dataset, 'restrecovery', measure)

    paired = sorted(set(awake_subj) & set(deep_subj) & set(rec_subj))
    a_idx = [awake_subj.index(s) for s in paired]
    d_idx = [deep_subj.index(s) for s in paired]
    r_idx = [rec_subj.index(s) for s in paired]
    n_rois = awake_mats.shape[-1]
    iu = np.triu_indices(n_rois, k=1)
    ij = np.array(iu)

    a_e = edge_stack(awake_mats[a_idx], iu)
    d_e = edge_stack(deep_mats[d_idx], iu)
    r_e = edge_stack(rec_mats[r_idx], iu)

    t_ad = paired_t(a_e - d_e)
    t_rd = paired_t(r_e - d_e)
    conj = (t_ad > T_THRESH) & (t_rd > T_THRESH)
    n_conj = int(conj.sum())
    print(f"Paired N={len(paired)} | conjunction edges={n_conj}", flush=True)

    sel = np.where(conj)[0]
    rows = ij[0, sel]
    cols = ij[1, sel]
    adj = csr_matrix(
        (np.ones(2 * sel.size, dtype=np.int8),
         (np.concatenate([rows, cols]), np.concatenate([cols, rows]))),
        shape=(n_rois, n_rois),
    )
    n_comp, node_labels = connected_components(adj, directed=False)
    edge_label = node_labels[rows]
    label_vals, edge_counts = np.unique(edge_label, return_counts=True)
    order = np.argsort(-edge_counts)
    label_vals = label_vals[order]
    edge_counts = edge_counts[order]

    components: list[dict] = []
    sig_threshold = 5
    for rank, (lab, n_e) in enumerate(zip(label_vals, edge_counts), start=1):
        edges_in_c = sel[edge_label == lab]
        nodes = np.unique(np.concatenate([ij[0, edges_in_c], ij[1, edges_in_c]]))
        components.append({
            'rank': rank, 'label': int(lab), 'n_edges': int(n_e),
            'n_nodes': int(nodes.size), 'node_indices': nodes.tolist(),
            'edge_indices': edges_in_c.tolist(),
        })

    print(f"\n{'='*70}\nConjunction components (top {sig_threshold})\n{'='*70}", flush=True)
    for c in components[:sig_threshold]:
        nodes = np.array(c['node_indices'])
        sub = roi_df.iloc[nodes].copy()
        sub['ws'] = workspace[nodes]
        n_in_ws = int(workspace[nodes].sum())
        net_counts = Counter(sub['Network_Name'].tolist())
        print(f"\nComponent #{c['rank']}: {c['n_edges']} edges, {c['n_nodes']} nodes "
              f"(workspace: {n_in_ws}/{c['n_nodes']})", flush=True)
        for net, cnt in sorted(net_counts.items(), key=lambda x: -x[1]):
            print(f"   {net:<14s} {cnt}", flush=True)
        for _, row in sub.iterrows():
            ws_tag = 'WS' if row['ws'] else '  '
            print(f"     ROI{int(row['ROI']):3d} [{ws_tag}] {row['Network_Name']:<14s} "
                  f"{row['Hemisphere']} {row['Structure_Name']}", flush=True)

    sig_nodes = sorted({n for c in components[:sig_threshold] for n in c['node_indices']})
    print(f"\n{'='*70}\nUnion of top-5 component nodes: {len(sig_nodes)} ROIs\n{'='*70}", flush=True)
    sig_arr = np.array(sig_nodes)
    sig_in_ws = int(workspace[sig_arr].sum())
    sig_in_dmn = int((rsn_labels[sig_arr] == 6).sum())
    print(f"  in HCP workspace : {sig_in_ws}/{len(sig_nodes)} "
          f"({100*sig_in_ws/len(sig_nodes):.1f}%)", flush=True)
    print(f"  in Default       : {sig_in_dmn}/{len(sig_nodes)} "
          f"({100*sig_in_dmn/len(sig_nodes):.1f}%)", flush=True)

    print(f"  Network breakdown:", flush=True)
    for net_id, net_name in enumerate(rsn_names):
        cnt = int((rsn_labels[sig_arr] == net_id).sum())
        total = int((rsn_labels == net_id).sum())
        print(f"    {net_name:<14s} {cnt:3d}/{total:3d} ({100*cnt/max(total,1):.1f}% of network's ROIs)",
              flush=True)

    print(f"\n{'='*70}\nEnrichment test (K={N_BOOT})\n{'='*70}", flush=True)
    rng = np.random.default_rng(SEED)
    boot_ws = np.zeros(N_BOOT, dtype=np.int32)
    boot_dmn = np.zeros(N_BOOT, dtype=np.int32)
    for b in range(N_BOOT):
        idx = rng.choice(n_rois, size=len(sig_nodes), replace=False)
        boot_ws[b] = int(workspace[idx].sum())
        boot_dmn[b] = int((rsn_labels[idx] == 6).sum())
    p_ws = float((boot_ws >= sig_in_ws).sum() + 1) / (N_BOOT + 1)
    p_dmn = float((boot_dmn >= sig_in_dmn).sum() + 1) / (N_BOOT + 1)
    print(f"  workspace overlap: observed={sig_in_ws}, "
          f"null mean={boot_ws.mean():.1f}, p_one_tail={p_ws:.4f}", flush=True)
    print(f"  DMN overlap      : observed={sig_in_dmn}, "
          f"null mean={boot_dmn.mean():.1f}, p_one_tail={p_dmn:.4f}", flush=True)

    rows_csv = []
    for c in components[:sig_threshold]:
        for n in c['node_indices']:
            rec = roi_df.iloc[n]
            rows_csv.append({
                'component': c['rank'],
                'roi_index_0based': int(n),
                'roi_index_1based': int(rec['ROI']),
                'network_id': int(rec['Network_ID']),
                'network_name': rec['Network_Name'],
                'hemisphere': rec['Hemisphere'],
                'structure': rec['Structure_Name'],
                'in_hcp_workspace': bool(workspace[n]),
            })
    out_csv = out_dir / f'nbs_{measure}_conjunction_top5_anatomy.csv'
    pd.DataFrame(rows_csv).to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}", flush=True)


if __name__ == '__main__':
    main()
