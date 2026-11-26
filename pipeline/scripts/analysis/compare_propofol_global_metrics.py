"""
Subject-level test: do PhiID metrics (synergy, redundancy, phi, phi_r)
decrease under propofol sedation in the Cambridge dataset?

Loads per-subject PhiID matrices (already on disk) via analysis_loader, which
applies the ART exclusion gate from cambridge_propofol.yaml. Computes whole-
brain mean over the upper triangle for each scan, pairs subjects across
conditions, and reports paired t-tests against restawake.

Luppi et al. (2024) prediction (eLife 12:e88173):
  synergy   : decrease under deep propofol
  phi       : decrease under deep propofol
  phi_r     : decrease under deep propofol
  redundancy: roughly preserved or slight increase
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel, ttest_ind, wilcoxon

sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from analysis_loader import load_group_all_measures

DATASET = 'cambridge_propofol'
CONTROL = 'restawake'
CONDITIONS = ['restlight', 'restdeep', 'restrecovery']
MEASURES = ['synergy', 'redundancy', 'phi', 'phi_r']

WORKSPACE_DATASET = 'hcp'
WORKSPACE_TASK = 'rest'


def upper_tri_mean(mat: np.ndarray, roi_mask: np.ndarray | None = None) -> float:
    n = mat.shape[0]
    if roi_mask is None:
        iu = np.triu_indices(n, k=1)
    else:
        ws_idx = np.where(roi_mask)[0]
        ii, jj = np.meshgrid(ws_idx, ws_idx, indexing='ij')
        keep = ii < jj
        iu = (ii[keep], jj[keep])
    vals = mat[iu]
    return float(np.nanmean(vals))


def per_scan_means(matrices: np.ndarray, roi_mask: np.ndarray | None = None) -> np.ndarray:
    return np.array([upper_tri_mean(m, roi_mask) for m in matrices])


def load_hcp_workspace_mask(project_dir: Path) -> np.ndarray:
    path = (project_dir / 'results' / 'group_analysis' / WORKSPACE_DATASET
            / WORKSPACE_TASK / f'{WORKSPACE_DATASET}_{WORKSPACE_TASK}_workspace_mask.npy')
    mask = np.load(path).astype(bool)
    print(f"Loaded HCP workspace mask: {path}")
    print(f"  Workspace ROIs: {int(mask.sum())} / {mask.size}")
    return mask


def build_rows(by_task: dict, scope_label: str, roi_mask: np.ndarray | None) -> list[dict]:
    rows: list[dict] = []
    for task in CONDITIONS:
        ctrl_subs = by_task[CONTROL]['subjects']
        cond_subs = by_task[task]['subjects']
        paired_subs = sorted(set(ctrl_subs) & set(cond_subs))
        ctrl_idx = [ctrl_subs.index(s) for s in paired_subs]
        cond_idx = [cond_subs.index(s) for s in paired_subs]

        for measure in MEASURES:
            ctrl_vals = per_scan_means(
                by_task[CONTROL]['matrices'][measure][ctrl_idx], roi_mask
            )
            cond_vals = per_scan_means(
                by_task[task]['matrices'][measure][cond_idx], roi_mask
            )
            delta = cond_vals - ctrl_vals
            t, p_t = ttest_rel(cond_vals, ctrl_vals)
            try:
                w_stat, p_w = wilcoxon(cond_vals, ctrl_vals)
            except ValueError:
                w_stat, p_w = (np.nan, np.nan)
            cohen_dz = (float(np.mean(delta) / np.std(delta, ddof=1))
                        if np.std(delta, ddof=1) > 0 else np.nan)

            rows.append({
                'scope': scope_label,
                'comparison': f'{task} vs {CONTROL}',
                'measure': measure,
                'n_paired': len(paired_subs),
                'mean_ctrl': float(np.mean(ctrl_vals)),
                'mean_cond': float(np.mean(cond_vals)),
                'mean_delta': float(np.mean(delta)),
                'pct_change': float(100 * np.mean(delta) / np.mean(ctrl_vals)),
                't_stat': float(t),
                'p_paired_t': float(p_t),
                'wilcoxon_W': float(w_stat) if w_stat is not None else np.nan,
                'p_wilcoxon': float(p_w) if p_w is not None else np.nan,
                'cohen_dz': cohen_dz,
                'direction': 'decrease' if np.mean(delta) < 0 else 'increase',
            })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', default=DATASET,
                        help='dataset id (default: %(default)s)')
    args = parser.parse_args()
    dataset = args.dataset

    project_dir = Path(__file__).resolve().parent.parent.parent
    workspace_mask = load_hcp_workspace_mask(project_dir)

    print(f"\n{'='*70}\nLoading per-scan PhiID matrices for {dataset}\n{'='*70}")
    by_task: dict[str, dict] = {}
    for task in [CONTROL] + CONDITIONS:
        print(f"\n-- {task} --")
        mats, subjects = load_group_all_measures(dataset, task)
        by_task[task] = {'subjects': subjects, 'matrices': mats, 'n': len(subjects)}
        print(f"   N = {len(subjects)} scans")

    rows = []
    rows += build_rows(by_task, scope_label='whole_brain', roi_mask=None)
    rows += build_rows(by_task, scope_label='hcp_workspace', roi_mask=workspace_mask)

    df = pd.DataFrame(rows)
    out_dir = Path('results/group_analysis') / dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / 'global_metric_comparison.csv'
    df.to_csv(csv_path, index=False, float_format='%.6f')

    pd.set_option('display.width', 200)
    pd.set_option('display.max_columns', None)

    for scope in ('whole_brain', 'hcp_workspace'):
        print(f"\n{'='*80}\nScope: {scope}\n{'='*80}")
        sub = df[df['scope'] == scope].drop(columns='scope')
        print(sub.to_string(index=False, float_format=lambda x: f'{x:.4f}'))

    print(f"\nSaved: {csv_path}")

    print(f"\n{'='*80}\nLuppi 2024 direction check (restdeep) — whole-brain vs HCP-workspace\n{'='*80}")
    expected = {'synergy': 'decrease', 'phi': 'decrease', 'phi_r': 'decrease',
                'redundancy': 'preserved/increase'}
    for measure, predicted in expected.items():
        line = f"  {measure:<10s}: predicted={predicted:<20s}"
        for scope in ('whole_brain', 'hcp_workspace'):
            row = df[(df['scope'] == scope) & (df['measure'] == measure)
                     & (df['comparison'] == 'restdeep vs restawake')].iloc[0]
            tag = ('OK' if (predicted == 'decrease' and row['direction'] == 'decrease')
                   or (predicted == 'preserved/increase')
                   else 'MISS')
            line += (f"  [{scope:<13s}] {row['direction']:<8s} "
                     f"({row['pct_change']:+.2f}%, p={row['p_paired_t']:.3g}) [{tag}]")
        print(line)


if __name__ == '__main__':
    main()
