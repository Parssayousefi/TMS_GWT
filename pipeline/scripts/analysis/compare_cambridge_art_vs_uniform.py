"""Phase 8 — ART vs uniform Cambridge preprocessing: how much did ART scrubbing change?

Compares the two coexisting Cambridge variants scan-by-scan and group-by-group:
  cambridge_propofol          ART denoise (functional_art FD>0.5mm / GS-Z>3 + scrubbing
                              regressors), 5 scans further dropped at >=20% outliers
  cambridge_propofol_uniform  standard denoise, identical to the other 5 datasets

Outputs (results/group_analysis/cambridge_comparison/):
  scan_level_similarity.csv    per scan x measure: upper-triangle Pearson r, mean shift,
                               joined to that scan's ART outlier fraction
  group_level_similarity.csv   per task x measure: r of group-average matrices, node-strength
                               rank correlation, workspace/gateway/broadcaster Jaccard
  summary.md                   human-readable verdict

The scan-level view answers the mechanistic question (does divergence track how many
volumes ART flagged?); the group-level view answers the publication question (do the
workspace and its gateway/broadcaster split survive the preprocessing change?).

Usage:
  python scripts/analysis/compare_cambridge_art_vs_uniform.py
  python scripts/analysis/compare_cambridge_art_vs_uniform.py --restrict-to-art-included
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

# Stored matrices carry a NaN diagonal, so nanmean over the scan stack sees all-NaN
# slices there. Expected, not a data problem.
warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_DIR / "scripts" / "utils"))

ART_DS = "cambridge_propofol"
UNI_DS = "cambridge_propofol_uniform"
TASKS = ["restawake", "restlight", "restdeep", "restrecovery"]
MEASURES = ["synergy", "redundancy", "phi", "phi_r"]
N_ROIS = 454

PHIID_DIR = PROJECT_DIR / "results" / "phiid_matrices"
ART_SUMMARY = PROJECT_DIR / "results" / "timeseries_intermediate" / ART_DS / "art_outlier_summary.csv"
ART_EXCL = PROJECT_DIR / "results" / "qc" / ART_DS / "exclusions_art.csv"
WORKSPACE_MASK = (PROJECT_DIR / "results" / "group_analysis" / "hcp" / "rest"
                  / "hcp_rest_workspace_mask.npy")


def matrix_path(dataset: str, scan_id: str, measure: str) -> Path:
    return PHIID_DIR / f"dataset-{dataset}_{scan_id}_{measure}.npy"


def discover_scan_ids() -> list[str]:
    """Scan ids present for BOTH variants (phi_r used as the sentinel measure)."""
    art = {p.name.split(f"dataset-{ART_DS}_")[1].rsplit("_phi_r.npy", 1)[0]
           for p in PHIID_DIR.glob(f"dataset-{ART_DS}_sub-*_phi_r.npy")}
    uni = {p.name.split(f"dataset-{UNI_DS}_")[1].rsplit("_phi_r.npy", 1)[0]
           for p in PHIID_DIR.glob(f"dataset-{UNI_DS}_sub-*_phi_r.npy")}
    return sorted(art & uni)


def upper_tri(m: np.ndarray) -> np.ndarray:
    iu = np.triu_indices_from(m, k=1)
    return m[iu]


def node_strength(m: np.ndarray) -> np.ndarray:
    m = m.copy()
    np.fill_diagonal(m, np.nan)
    with np.errstate(invalid="ignore"):
        s = np.nansum(m, axis=1)
    s[np.all(~np.isfinite(m), axis=1)] = np.nan
    return s


def workspace_from(syn: np.ndarray, red: np.ndarray) -> np.ndarray:
    """Reproduce 05_group_average_workspace.identify_workspace: rank(strength) sign test."""
    sr = rankdata(node_strength(syn), nan_policy="omit")
    rr = rankdata(node_strength(red), nan_policy="omit")
    return (sr - rr) > 0


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    inter = np.logical_and(a, b).sum()
    union = np.logical_or(a, b).sum()
    return float(inter / union) if union else np.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--restrict-to-art-included", action="store_true",
                    help="drop the 5 scans ART excluded (>=20%% outliers) so the comparison "
                         "isolates the denoising change from the sample change")
    ap.add_argument("--out-dir", default=str(PROJECT_DIR / "results" / "group_analysis"
                                             / "cambridge_comparison"))
    args = ap.parse_args()

    scan_ids = discover_scan_ids()
    print(f"scans present in both variants: {len(scan_ids)}")

    art_excluded: set[str] = set()
    if ART_EXCL.exists():
        # Same read as qc_loader._load_file: '#' provenance header, subject kept as str.
        # Schema is subject/task, not scan_id, so rebuild the scan id (always run-1 here).
        ex = pd.read_csv(ART_EXCL, dtype={"subject": str}, comment="#")
        ex = ex[ex["status"].astype(str).str.upper().eq("EXCLUDE")]
        art_excluded = {f"sub-{r.subject}_task-{r.task}_run-1" for r in ex.itertuples()}
        print(f"ART-excluded scans: {len(art_excluded)} -> {sorted(art_excluded)}")
    if args.restrict_to_art_included:
        scan_ids = [s for s in scan_ids if s not in art_excluded]
        print(f"restricted to ART-included: {len(scan_ids)} scans")

    outlier_frac = {}
    if ART_SUMMARY.exists():
        s = pd.read_csv(ART_SUMMARY)
        outlier_frac = dict(zip(s["scan_id"], s["fraction_outliers"]))

    # ---------------- scan level ----------------
    rows = []
    for scan_id in scan_ids:
        task = scan_id.split("task-")[1].split("_")[0]
        for measure in MEASURES:
            pa, pu = matrix_path(ART_DS, scan_id, measure), matrix_path(UNI_DS, scan_id, measure)
            if not (pa.exists() and pu.exists()):
                continue
            a, u = upper_tri(np.load(pa)), upper_tri(np.load(pu))
            m = np.isfinite(a) & np.isfinite(u)
            r = float(np.corrcoef(a[m], u[m])[0, 1]) if m.sum() > 2 else np.nan
            rows.append({
                "scan_id": scan_id, "task": task, "measure": measure,
                "r_upper_tri": r,
                "mean_art": float(np.nanmean(a)), "mean_uniform": float(np.nanmean(u)),
                "mean_shift_pct": float(100 * (np.nanmean(u) - np.nanmean(a)) / np.nanmean(a))
                if np.nanmean(a) else np.nan,
                "art_fraction_outliers": outlier_frac.get(scan_id, np.nan),
                "art_excluded": scan_id in art_excluded,
            })
    scan_df = pd.DataFrame(rows)

    # ---------------- group level ----------------
    workspace_hcp = np.load(WORKSPACE_MASK).astype(bool) if WORKSPACE_MASK.exists() else None
    grows = []
    for task in TASKS:
        ids = [s for s in scan_ids if f"task-{task}_" in s]
        if not ids:
            continue
        avg = {}
        for ds in (ART_DS, UNI_DS):
            for measure in MEASURES:
                stack = [np.load(matrix_path(ds, s, measure)) for s in ids
                         if matrix_path(ds, s, measure).exists()]
                avg[(ds, measure)] = np.nanmean(np.stack(stack), axis=0) if stack else None
        for measure in MEASURES:
            A, U = avg[(ART_DS, measure)], avg[(UNI_DS, measure)]
            if A is None or U is None:
                continue
            a, u = upper_tri(A), upper_tri(U)
            m = np.isfinite(a) & np.isfinite(u)
            sa, su = node_strength(A), node_strength(U)
            ms = np.isfinite(sa) & np.isfinite(su)
            grows.append({
                "task": task, "measure": measure, "n_scans": len(ids),
                "r_group_matrix": float(np.corrcoef(a[m], u[m])[0, 1]),
                "spearman_node_strength": float(spearmanr(sa[ms], su[ms]).statistic),
            })
        # workspace / gateway-broadcaster overlap needs synergy+redundancy together
        if avg[(ART_DS, "synergy")] is not None and avg[(UNI_DS, "synergy")] is not None:
            wa = workspace_from(avg[(ART_DS, "synergy")], avg[(ART_DS, "redundancy")])
            wu = workspace_from(avg[(UNI_DS, "synergy")], avg[(UNI_DS, "redundancy")])
            grows.append({
                "task": task, "measure": "workspace_mask", "n_scans": len(ids),
                "jaccard_dataset_defined_workspace": jaccard(wa, wu),
                "n_workspace_art": int(wa.sum()), "n_workspace_uniform": int(wu.sum()),
            })
            if workspace_hcp is not None:
                grows[-1]["jaccard_vs_hcp_art"] = jaccard(wa, workspace_hcp)
                grows[-1]["jaccard_vs_hcp_uniform"] = jaccard(wu, workspace_hcp)
    group_df = pd.DataFrame(grows)

    out_dir = Path(args.out_dir)
    if args.restrict_to_art_included:
        out_dir = out_dir / "art_included_only"
    out_dir.mkdir(parents=True, exist_ok=True)
    scan_df.to_csv(out_dir / "scan_level_similarity.csv", index=False)
    group_df.to_csv(out_dir / "group_level_similarity.csv", index=False)

    # ---------------- report ----------------
    pd.set_option("display.width", 170)
    print("\n=== scan-level r (upper triangle), by measure ===")
    print(scan_df.groupby("measure")["r_upper_tri"]
          .agg(["count", "min", "median", "max"]).to_string(float_format=lambda v: f"{v:.4f}"))

    print("\n=== does divergence track ART outlier fraction? ===")
    lines = []
    for measure in MEASURES:
        d = scan_df[(scan_df.measure == measure)].dropna(subset=["r_upper_tri", "art_fraction_outliers"])
        if len(d) > 3:
            rho = spearmanr(d.art_fraction_outliers, d.r_upper_tri).statistic
            line = (f"{measure:12s} spearman(outlier_frac, r) = {rho:+.3f}  "
                    f"(n={len(d)}; negative = more scrubbing -> more divergence)")
            print(line)
            lines.append(line)

    print("\n=== group-level ===")
    print(group_df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    with open(out_dir / "summary.md", "w", encoding="utf-8") as f:
        f.write("# Cambridge ART vs uniform preprocessing comparison\n\n")
        f.write(f"- scans compared: {len(scan_ids)}"
                f"{' (ART-included only)' if args.restrict_to_art_included else ''}\n")
        f.write(f"- ART-excluded scans present on disk: {sorted(art_excluded)}\n\n")
        f.write("## Scan-level similarity (upper-triangle Pearson r)\n\n```\n")
        f.write(scan_df.groupby("measure")["r_upper_tri"]
                .agg(["count", "min", "median", "max"]).to_string())
        f.write("\n```\n\n## Divergence vs ART outlier fraction\n\n```\n")
        f.write("\n".join(lines))
        f.write("\n```\n\n## Group level\n\n```\n")
        f.write(group_df.to_string(index=False))
        f.write("\n```\n")
    print(f"\nwritten: {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
