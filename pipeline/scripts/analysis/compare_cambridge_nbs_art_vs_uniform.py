"""Phase 8 (NBS part) — do the network-level propofol findings survive the switch
from ART-scrubbed to uniform preprocessing?

Reads the nbs_{measure}_propofol_full.yaml written by nbs_propofol_phir.py for both
Cambridge variants and reports, per measure and contrast: how many components survive
FWE correction, the best corrected p, the size of the largest surviving component, and
the paired N. A contrast is called SAME / ONLY-ART / ONLY-UNIFORM on the p<.05 verdict,
but the printed p-values matter more than the label — a .0398 vs .0504 pair is the same
evidence either side of an arbitrary line, not a changed finding.

Usage:
  python scripts/analysis/compare_cambridge_nbs_art_vs_uniform.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
BASE = PROJECT_DIR / "results" / "group_analysis"
VARIANTS = {"ART": "cambridge_propofol", "UNIFORM": "cambridge_propofol_uniform"}
MEASURES = ["synergy", "redundancy", "phi_r", "phi"]
ALPHA = 0.05


def load_meta(dataset: str, measure: str) -> dict | None:
    p = BASE / dataset / "nbs" / measure / f"nbs_{measure}_propofol_full.yaml"
    if not p.exists():
        return None
    with open(p) as f:
        return yaml.safe_load(f)


def paired_n(meta: dict) -> dict:
    return {k: len(v) for k, v in (meta.get("paired_subjects") or {}).items()}


def rows_for(meta: dict, variant: str, measure: str) -> list[dict]:
    npair = paired_n(meta)
    # map a tail tag to whichever paired-subject list it came from
    def n_for(tag: str):
        for key, n in npair.items():
            stem = key.split("_")[0]          # fig4B / fig4C / fig4D
            if tag.startswith(stem):
                return n
        return None

    out = []
    for tail in meta.get("tails", []):
        comps = tail.get("components") or []
        sig = [c for c in comps if c.get("p_corrected", 1.0) < ALPHA]
        best = min(comps, key=lambda c: c.get("p_corrected", 1.0)) if comps else {}
        largest_sig = max(sig, key=lambda c: c.get("edges", 0)) if sig else {}
        out.append({
            "measure": measure,
            "contrast": tail.get("tag", "?"),
            "variant": variant,
            "n_paired": n_for(tail.get("tag", "")),
            "n_components": tail.get("n_components", len(comps)),
            "n_significant": len(sig),
            "best_p_corrected": best.get("p_corrected"),
            "largest_sig_edges": largest_sig.get("edges", 0),
            "largest_sig_nodes": largest_sig.get("nodes", 0),
        })
    return out


def main() -> int:
    records = []
    for measure in MEASURES:
        for variant, ds in VARIANTS.items():
            meta = load_meta(ds, measure)
            if meta is None:
                continue
            records.append(pd.DataFrame(rows_for(meta, variant, measure)))
    if not records:
        print("no NBS metadata found")
        return 1
    df = pd.concat(records, ignore_index=True)

    pd.set_option("display.width", 200)
    for measure in MEASURES:
        sub = df[df.measure == measure]
        if sub.empty:
            continue
        print(f"\n{'=' * 96}\n{measure.upper()}\n{'=' * 96}")
        wide = sub.pivot(index="contrast", columns="variant",
                         values=["n_paired", "n_significant", "best_p_corrected",
                                 "largest_sig_edges"])
        print(wide.to_string(float_format=lambda v: f"{v:.4g}"))
        if {"ART", "UNIFORM"}.issubset(set(sub.variant)):
            print("\n  verdicts:")
            for contrast in sorted(sub.contrast.unique()):
                a = sub[(sub.contrast == contrast) & (sub.variant == "ART")]
                u = sub[(sub.contrast == contrast) & (sub.variant == "UNIFORM")]
                if a.empty or u.empty:
                    continue
                sa = int(a.n_significant.iloc[0]) > 0
                su = int(u.n_significant.iloc[0]) > 0
                verdict = "SAME" if sa == su else ("ONLY-UNIFORM" if su else "ONLY-ART")
                pa, pu = a.best_p_corrected.iloc[0], u.best_p_corrected.iloc[0]
                note = ""
                if verdict != "SAME" and pa is not None and pu is not None:
                    if min(pa, pu) > 0.01 and max(pa, pu) < 0.15:
                        note = "  <- borderline: same evidence, opposite sides of .05"
                print(f"    {contrast:24s} ART p={pa:.4f} / UNIFORM p={pu:.4f}"
                      f"   -> {verdict}{note}")
        else:
            missing = {"ART", "UNIFORM"} - set(sub.variant)
            print(f"\n  (no run for: {', '.join(sorted(missing))} — comparison incomplete)")

    out_dir = BASE / "cambridge_comparison"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "nbs_art_vs_uniform.csv", index=False)
    print(f"\nwritten: {out_dir / 'nbs_art_vs_uniform.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
