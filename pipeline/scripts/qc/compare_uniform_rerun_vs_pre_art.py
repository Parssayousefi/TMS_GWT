"""Gate 6 — fidelity check: does the 2026-08 uniform rerun reproduce the ORIGINAL
uniform (no-ART) Cambridge run that the ART variant overwrote on 2026-04-28?

The only surviving numeric trace of the original run is the per-scan summary written
on 2026-04-27 to results/qc/cambridge_propofol/motion_synergy_pre_art/per_scan.csv
(see docs/PLAN_cambridge_uniform_rerun.md, "fidelity oracle"). This script compares
that snapshot against the equivalent summary computed from the rerun.

Inputs
  --pre-art   per_scan.csv from the April pre-ART snapshot (the oracle)
  --rerun     per_scan.csv produced by
              scripts/qc/cambridge_motion_synergy_baseline.py
                  --dataset cambridge_propofol_uniform --label uniform_rerun

Motion columns (fd_*, trans_max_mm, rot_max_deg) are denoising-independent and act as
a control: they must match essentially exactly. PhiID columns are the actual test —
CONN/SPM segmentation is not bit-deterministic, so we judge on correlation and
relative magnitude rather than equality.

Exit code 0 = PASS, 1 = FAIL.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MOTION_COLS = [
    "fd_mean", "fd_max", "fd_frac_gt_0p2", "fd_frac_gt_0p5",
    "trans_max_mm", "rot_max_deg",
]
PHIID_COLS = [
    "synergy_brainwide_mean", "synergy_workspace_mean", "synergy_workspace_node_strength_mean",
    "redundancy_brainwide_mean", "redundancy_workspace_mean", "redundancy_workspace_node_strength_mean",
    "phi_brainwide_mean", "phi_workspace_mean", "phi_workspace_node_strength_mean",
    "phi_r_brainwide_mean", "phi_r_workspace_mean", "phi_r_workspace_node_strength_mean",
]

# Gate thresholds (plan Gate 6)
MIN_R_PHIID = 0.99
MAX_MEDIAN_RELDIFF = 0.01      # 1 % typical deviation
MAX_MOTION_ABSDIFF = 1e-6      # motion must be reproduced exactly


def compare_column(a: pd.Series, b: pd.Series) -> dict:
    mask = np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 3:
        return {"n": n, "r": np.nan, "median_reldiff": np.nan, "max_absdiff": np.nan}
    x, y = a[mask].to_numpy(float), b[mask].to_numpy(float)
    r = float(np.corrcoef(x, y)[0, 1]) if np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan
    denom = np.where(np.abs(x) > 0, np.abs(x), np.nan)
    reldiff = np.abs(y - x) / denom
    return {
        "n": n,
        "r": r,
        "median_reldiff": float(np.nanmedian(reldiff)),
        "max_absdiff": float(np.max(np.abs(y - x))),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pre-art", default="results/qc/cambridge_propofol/motion_synergy_pre_art/per_scan.csv")
    ap.add_argument("--rerun", default="results/qc/cambridge_propofol_uniform/motion_synergy_uniform_rerun/per_scan.csv")
    ap.add_argument("--out-dir", default="results/qc/cambridge_propofol_uniform/fidelity_gate6")
    args = ap.parse_args()

    pre = pd.read_csv(args.pre_art)
    new = pd.read_csv(args.rerun)
    print(f"pre-ART snapshot : {args.pre_art}  ({len(pre)} rows)")
    print(f"uniform rerun    : {args.rerun}  ({len(new)} rows)")

    merged = pre.merge(new, on="scan_id", suffixes=("_pre", "_new"), how="inner")
    print(f"matched scans    : {len(merged)}")

    rows = []
    for col in MOTION_COLS + PHIID_COLS:
        cp, cn = f"{col}_pre", f"{col}_new"
        if cp not in merged or cn not in merged:
            print(f"  WARNING: column missing, skipped: {col}")
            continue
        stats = compare_column(merged[cp], merged[cn])
        stats["column"] = col
        stats["kind"] = "motion" if col in MOTION_COLS else "phiid"
        rows.append(stats)

    res = pd.DataFrame(rows)[["kind", "column", "n", "r", "median_reldiff", "max_absdiff"]]
    pd.set_option("display.width", 160)
    print("\n" + res.to_string(index=False, float_format=lambda v: f"{v:.6g}"))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_dir / "gate6_fidelity_summary.csv", index=False)
    merged.to_csv(out_dir / "gate6_merged_per_scan.csv", index=False)

    motion = res[res.kind == "motion"]
    phiid = res[res.kind == "phiid"]

    motion_ok = bool((motion["max_absdiff"] <= MAX_MOTION_ABSDIFF).all())
    r_ok = bool((phiid["r"] >= MIN_R_PHIID).all())
    mag_ok = bool((phiid["median_reldiff"] <= MAX_MEDIAN_RELDIFF).all())

    print(f"\nmotion reproduced exactly (<= {MAX_MOTION_ABSDIFF:g}): {motion_ok}")
    print(f"all PhiID r >= {MIN_R_PHIID}:                      {r_ok}"
          f"   (min r = {phiid['r'].min():.6f})")
    print(f"all PhiID median reldiff <= {MAX_MEDIAN_RELDIFF:.1%}:        {mag_ok}"
          f"   (max median reldiff = {phiid['median_reldiff'].max():.4%})")

    verdict = motion_ok and r_ok and mag_ok
    print(f"\nGATE 6: {'PASS' if verdict else 'FAIL'}")
    print(f"written: {out_dir / 'gate6_fidelity_summary.csv'}")
    if not verdict:
        print("\nDo NOT proceed to the publication analyses (plan Phase 7) until this is "
              "understood: the rerun does not reproduce the original uniform run.")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
