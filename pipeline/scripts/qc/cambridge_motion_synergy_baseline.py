#!/usr/bin/env python3
"""
Cambridge propofol — per-scan motion ↔ synergy/redundancy/phi_r baseline.

Snapshots the relationship between motion and PhiID measures at the current
state of the pipeline so we can compare pre-ART vs post-ART side-by-side.

Run BEFORE deleting old cambridge results:
    python scripts/qc/cambridge_motion_synergy_baseline.py --label pre_art

Run AGAIN after the ART rerun completes:
    python scripts/qc/cambridge_motion_synergy_baseline.py --label post_art

Outputs (under results/qc/cambridge_propofol/motion_synergy_<label>/):
    per_scan.csv          — one row per scan with motion + synergy summary stats
    correlations.csv      — per-task Pearson r between every (motion, synergy) pair
                            with paired-cohort sample size and p-value
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

PROJECT_DIR = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_DIR / "scripts" / "utils"))

from dataset_handler import load_dataset_handler  # noqa: E402

DATASET = "cambridge_propofol"
N_ROIS = 454
PHIID_MEASURES = ("synergy", "redundancy", "phi", "phi_r")
WORKSPACE_MASK_PATH = PROJECT_DIR / "results" / "group_analysis" / "hcp" / "rest" / "hcp_rest_workspace_mask.npy"

# Power 2012 framewise displacement: |Δtrans_x| + |Δtrans_y| + |Δtrans_z|
# + 50 * (|Δrot_x| + |Δrot_y| + |Δrot_z|), rotations converted from rad to mm
# on a 50 mm sphere.
def compute_fd_power(confounds: pd.DataFrame) -> np.ndarray:
    cols = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]
    miss = [c for c in cols if c not in confounds.columns]
    if miss:
        raise KeyError(f"missing motion columns: {miss}")
    motion = confounds[cols].to_numpy(dtype=float)
    motion[:, 3:6] *= 50.0
    delta = np.diff(motion, axis=0, prepend=motion[:1])
    return np.abs(delta).sum(axis=1)


def scan_summary(scan: dict, workspace: np.ndarray, dataset: str) -> dict | None:
    subject = scan["subject"]
    task = scan["task"]
    run = scan.get("run", 1) or 1
    ses = scan.get("session")
    ses_str = f"_ses-{ses}" if ses else ""
    scan_id = f"sub-{subject}{ses_str}_task-{task}_run-{run}"

    confounds_path = Path(scan["confounds_file"]) if scan.get("confounds_file") else None
    if confounds_path is None or not confounds_path.exists():
        return None
    confounds = pd.read_csv(confounds_path, sep="\t")
    fd = compute_fd_power(confounds)
    fd_mean = float(np.nanmean(fd))
    fd_max = float(np.nanmax(fd))
    fd_frac_05 = float(np.mean(fd > 0.5))
    fd_frac_02 = float(np.mean(fd > 0.2))

    trans = confounds[["trans_x", "trans_y", "trans_z"]].to_numpy(dtype=float)
    rot = confounds[["rot_x", "rot_y", "rot_z"]].to_numpy(dtype=float)
    trans_max = float(np.max(np.abs(trans)))
    rot_max = float(np.max(np.abs(rot)) * 180.0 / np.pi)

    out = dict(
        dataset=dataset,
        scan_id=scan_id,
        subject=subject,
        task=task,
        run=run,
        n_volumes=int(len(confounds)),
        fd_mean=fd_mean,
        fd_max=fd_max,
        fd_frac_gt_0p2=fd_frac_02,
        fd_frac_gt_0p5=fd_frac_05,
        trans_max_mm=trans_max,
        rot_max_deg=rot_max,
    )

    phiid_dir = PROJECT_DIR / "results" / "phiid_matrices"
    for measure in PHIID_MEASURES:
        f = phiid_dir / f"dataset-{dataset}_{scan_id}_{measure}.npy"
        if not f.exists():
            out[f"{measure}_brainwide_mean"] = np.nan
            out[f"{measure}_workspace_mean"] = np.nan
            out[f"{measure}_workspace_node_strength_mean"] = np.nan
            continue
        m = np.load(f)
        if m.shape != (N_ROIS, N_ROIS):
            out[f"{measure}_brainwide_mean"] = np.nan
            out[f"{measure}_workspace_mean"] = np.nan
            out[f"{measure}_workspace_node_strength_mean"] = np.nan
            continue
        triu = np.triu_indices(N_ROIS, k=1)
        out[f"{measure}_brainwide_mean"] = float(np.nanmean(m[triu]))

        ws_idx = np.where(workspace)[0]
        ws_block = m[np.ix_(ws_idx, ws_idx)]
        ws_triu = np.triu_indices(len(ws_idx), k=1)
        out[f"{measure}_workspace_mean"] = float(np.nanmean(ws_block[ws_triu]))

        node_strength = np.nansum(m, axis=1)
        out[f"{measure}_workspace_node_strength_mean"] = float(np.nanmean(node_strength[ws_idx]))

    return out


def correlations_per_task(per_scan: pd.DataFrame) -> pd.DataFrame:
    motion_cols = ["fd_mean", "fd_max", "fd_frac_gt_0p2", "fd_frac_gt_0p5", "trans_max_mm", "rot_max_deg"]
    synergy_cols = [c for c in per_scan.columns if any(c.startswith(m + "_") for m in PHIID_MEASURES)]

    rows = []
    for task, group in per_scan.groupby("task", sort=False):
        for mcol in motion_cols:
            for scol in synergy_cols:
                pair = group[[mcol, scol]].dropna()
                n = len(pair)
                if n < 3:
                    rows.append(dict(task=task, motion_var=mcol, synergy_var=scol,
                                     n=n, pearson_r=np.nan, p_value=np.nan))
                    continue
                if pair[mcol].std(ddof=0) == 0 or pair[scol].std(ddof=0) == 0:
                    rows.append(dict(task=task, motion_var=mcol, synergy_var=scol,
                                     n=n, pearson_r=np.nan, p_value=np.nan))
                    continue
                r, p = pearsonr(pair[mcol], pair[scol])
                rows.append(dict(task=task, motion_var=mcol, synergy_var=scol,
                                 n=n, pearson_r=float(r), p_value=float(p)))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True,
                        help="suffix for output dir, e.g. 'pre_art' or 'post_art'")
    parser.add_argument('--dataset', default=DATASET,
                        help='dataset id (default: %(default)s)')
    args = parser.parse_args()
    dataset = args.dataset

    if not WORKSPACE_MASK_PATH.exists():
        sys.exit(f"workspace mask missing: {WORKSPACE_MASK_PATH}")
    workspace = np.load(WORKSPACE_MASK_PATH).astype(bool)
    if workspace.shape != (N_ROIS,):
        sys.exit(f"workspace mask wrong shape: {workspace.shape}")
    print(f"workspace size: {int(workspace.sum())} ROIs")

    handler = load_dataset_handler(dataset)
    scans = handler.discover_all_scans()
    print(f"discovered {len(scans)} scans for {dataset}")

    rows = []
    for s in scans:
        try:
            row = scan_summary(s, workspace, dataset)
        except Exception as e:
            print(f"  skip sub-{s['subject']} task-{s['task']}: {e}")
            continue
        if row is None:
            continue
        rows.append(row)

    if not rows:
        sys.exit("no scans summarised — nothing to write")

    per_scan = pd.DataFrame(rows).sort_values(["task", "subject"]).reset_index(drop=True)
    correlations = correlations_per_task(per_scan)

    out_dir = PROJECT_DIR / "results" / "qc" / dataset / f"motion_synergy_{args.label}"
    out_dir.mkdir(parents=True, exist_ok=True)
    per_scan_path = out_dir / "per_scan.csv"
    corr_path = out_dir / "correlations.csv"
    per_scan.to_csv(per_scan_path, index=False)
    correlations.to_csv(corr_path, index=False)

    print(f"\nwrote {len(per_scan)} per-scan rows -> {per_scan_path}")
    print(f"wrote {len(correlations)} correlation rows -> {corr_path}")

    print("\nstrongest |r| per task (top 5):")
    top = (correlations.dropna(subset=["pearson_r"])
           .assign(abs_r=lambda d: d["pearson_r"].abs())
           .sort_values(["task", "abs_r"], ascending=[True, False])
           .groupby("task")
           .head(5))
    print(top.to_string(index=False))


if __name__ == "__main__":
    main()
