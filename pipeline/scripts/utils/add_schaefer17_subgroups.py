#!/usr/bin/env python3
"""
Augment roi_to_network_mapping.csv with paper-style label columns.
=================================================================

Adds three columns to atlases/roi_to_network_mapping.csv WITHOUT changing
the existing Network_ID / Network_Name / Structure_Name columns:

    Network17_ID         Schaefer 17-net parcel ID (1..400); NA for subcortical
    Network17_Subgroup   17-net subgroup token (VisCent, ContA, DefaultB,
                         TempPar, ...). NA for subcortical
    Substructure         Fine subcortical label (HIP, AMY, NAc, PUT, CAU,
                         GP, vTHA, dTHA). NA for cortical

The 7-net analysis pipeline (PhiID, workspace, gateways/broadcasters) is
unaffected: the new columns are for visualization only. Do NOT use
Network17_* for participation-coefficient module partition — Luppi locks
PC computation to the 8-module partition (Network_ID 0..7).

Idempotent. Rerunning a second time is a no-op.

Method (cortical):
    1. Fetch Schaefer 400/7-net and 400/17-net NIfTIs from nilearn (same
       call create_roi_network_mapping.py uses for 7-net, so parcel IDs in
       the existing ROI column match atlas7 by construction).
    2. For each 7-net parcel pid7 in 1..400, find the MODAL 17-net parcel
       ID at the voxels where data7 == pid7. Modal handles the few
       boundary voxels where 7-net and 17-net parcellations disagree.
    3. Parse the subgroup token from atlas17.labels[pid17 - 1] as
       label.split('_')[2] (robust to LimbicA_TempPole vs SomMotA forms).

Method (subcortical):
    Pattern-match on ROI_Label (Tian S4) to derive Substructure.
"""

from __future__ import annotations

import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
CSV_PATH = BASE_DIR / "atlases" / "roi_to_network_mapping.csv"

N_CORTICAL = 400
N_SUBCORTICAL = 54

EXPECTED_SUBSTRUCTURE_COUNTS = {
    "HIP": 10, "AMY": 4, "NAc": 4, "PUT": 8,
    "CAU": 8, "GP": 4, "vTHA": 10, "dTHA": 6,
}


def parse_substructure(roi_label: str) -> str | None:
    """Map a Tian S4 label like 'PUT-VA-rh' to a coarse paper-style name."""
    if not isinstance(roi_label, str):
        return None
    # Strip hemisphere suffix (-rh / -lh) if present.
    stem = roi_label[:-3] if roi_label.endswith(("-rh", "-lh")) else roi_label
    if stem.startswith("HIP"):
        return "HIP"
    if stem in ("lAMY", "mAMY"):
        return "AMY"
    if stem.startswith("NAc"):
        return "NAc"
    if stem.startswith("PUT"):
        return "PUT"
    if stem.startswith("CAU"):
        return "CAU"
    if stem in ("pGP", "aGP"):
        return "GP"
    if stem.startswith("THA-V"):
        return "vTHA"
    if stem.startswith("THA-D"):
        return "dTHA"
    return None


def build_7to17_mapping() -> tuple[dict[int, int], list[str]]:
    """Return {pid7: pid17_modal} and the 17-net label list (index pid17-1)."""
    try:
        from nilearn import datasets, image
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "nilearn is required; install with `pip install nilearn`"
        ) from exc

    print("Fetching Schaefer 400 / 7-net via nilearn ...")
    atlas7 = datasets.fetch_atlas_schaefer_2018(
        n_rois=N_CORTICAL, yeo_networks=7, resolution_mm=2, verbose=0
    )
    print("Fetching Schaefer 400 / 17-net via nilearn ...")
    atlas17 = datasets.fetch_atlas_schaefer_2018(
        n_rois=N_CORTICAL, yeo_networks=17, resolution_mm=2, verbose=0
    )

    img7 = image.load_img(atlas7.maps)
    img17 = image.load_img(atlas17.maps)
    if img7.shape != img17.shape or not np.allclose(img7.affine, img17.affine):
        print("  affine/shape mismatch — resampling 17-net onto 7-net grid")
        img17 = image.resample_to_img(img17, img7, interpolation="nearest")

    data7 = np.asarray(img7.dataobj).astype(np.int32)
    data17 = np.asarray(img17.dataobj).astype(np.int32)

    print(f"  data7  unique IDs (excl 0): {len(np.unique(data7)) - 1}")
    print(f"  data17 unique IDs (excl 0): {len(np.unique(data17)) - 1}")

    pid7_to_pid17: dict[int, int] = {}
    for pid7 in range(1, N_CORTICAL + 1):
        mask = data7 == pid7
        if not mask.any():
            raise RuntimeError(
                f"7-net parcel {pid7} has no voxels in the atlas NIfTI"
            )
        ids17 = data17[mask]
        ids17 = ids17[ids17 > 0]
        if ids17.size == 0:
            raise RuntimeError(
                f"7-net parcel {pid7} has zero 17-net coverage (all background)"
            )
        counts = Counter(ids17.tolist())
        pid7_to_pid17[pid7] = counts.most_common(1)[0][0]

    labels17 = [
        lbl.decode() if isinstance(lbl, (bytes, bytearray)) else lbl
        for lbl in atlas17.labels
    ]
    if labels17 and labels17[0] == "Background":
        labels17 = labels17[1:]
    if len(labels17) != N_CORTICAL:
        raise RuntimeError(
            f"Expected {N_CORTICAL} 17-net labels, got {len(labels17)}"
        )

    return pid7_to_pid17, labels17


def is_already_done(df: pd.DataFrame) -> bool:
    """Idempotency guard: True if all 3 columns exist with expected coverage."""
    required = {"Network17_ID", "Network17_Subgroup", "Substructure"}
    if not required.issubset(df.columns):
        return False
    cortical = df[(df["ROI"] >= 1) & (df["ROI"] <= N_CORTICAL)]
    subcort = df[(df["ROI"] > N_CORTICAL) & (df["ROI"] <= N_CORTICAL + N_SUBCORTICAL)]
    return (
        cortical["Network17_ID"].notna().sum() == N_CORTICAL
        and cortical["Network17_Subgroup"].notna().sum() == N_CORTICAL
        and subcort["Substructure"].notna().sum() == N_SUBCORTICAL
    )


def main() -> int:
    if not CSV_PATH.exists():
        print(f"ERROR: CSV not found at {CSV_PATH}", file=sys.stderr)
        print(
            "Run scripts/utils/create_roi_network_mapping.py first.",
            file=sys.stderr,
        )
        return 1

    df = pd.read_csv(CSV_PATH)
    print(f"Loaded {CSV_PATH}  ({len(df)} rows, {len(df.columns)} cols)")

    if is_already_done(df):
        print("Network17_ID / Network17_Subgroup / Substructure already populated.")
        print("Already up to date, no changes written.")
        return 0

    pid7_to_pid17, labels17 = build_7to17_mapping()
    print(f"Built 7->17 parcel mapping for {len(pid7_to_pid17)} / {N_CORTICAL} cortical parcels")

    df["Network17_ID"] = pd.array([pd.NA] * len(df), dtype="Int64")
    df["Network17_Subgroup"] = pd.array([pd.NA] * len(df), dtype="string")
    df["Substructure"] = pd.array([pd.NA] * len(df), dtype="string")

    cortical_mask = df["ROI"].between(1, N_CORTICAL)
    subcort_mask = df["ROI"].between(N_CORTICAL + 1, N_CORTICAL + N_SUBCORTICAL)

    for idx, row in df[cortical_mask].iterrows():
        pid7 = int(row["ROI"])
        pid17 = pid7_to_pid17[pid7]
        label17 = labels17[pid17 - 1]
        tokens = label17.split("_")
        if len(tokens) < 3 or tokens[0] != "17Networks":
            raise RuntimeError(f"Unexpected 17-net label format: {label17!r}")
        subgroup = tokens[2]
        df.at[idx, "Network17_ID"] = pid17
        df.at[idx, "Network17_Subgroup"] = subgroup

    n_parsed = 0
    for idx, row in df[subcort_mask].iterrows():
        sub = parse_substructure(row["ROI_Label"])
        if sub is None:
            print(
                f"WARNING: unparsed subcortical label at ROI {row['ROI']}: "
                f"{row['ROI_Label']!r}",
                file=sys.stderr,
            )
            continue
        df.at[idx, "Substructure"] = sub
        n_parsed += 1
    print(f"Parsed Substructure for {n_parsed} / {N_SUBCORTICAL} subcortical ROIs")

    print("\nCortical 17-net subgroup distribution:")
    cort_counts = (
        df.loc[cortical_mask, "Network17_Subgroup"].value_counts().sort_index()
    )
    for name, n in cort_counts.items():
        print(f"  {name:<15s} {n:>4d}")
    print(f"  {'TOTAL':<15s} {cort_counts.sum():>4d}")

    print("\nSubcortical substructure distribution:")
    sub_counts = (
        df.loc[subcort_mask, "Substructure"].value_counts().sort_index()
    )
    for name, n in sub_counts.items():
        exp = EXPECTED_SUBSTRUCTURE_COUNTS.get(name)
        marker = "" if exp is None or n == exp else f"  (expected {exp})"
        print(f"  {name:<6s} {n:>3d}{marker}")
    print(f"  {'TOTAL':<6s} {sub_counts.sum():>3d}")

    tmp_path = CSV_PATH.with_suffix(".csv.tmp")
    df.to_csv(tmp_path, index=False)
    try:
        os.replace(tmp_path, CSV_PATH)
    except PermissionError:
        side_path = CSV_PATH.with_name(CSV_PATH.stem + ".augmented.csv")
        os.replace(tmp_path, side_path)
        print(
            f"\nCSV is locked (likely open in Excel). Wrote augmented copy to:"
            f"\n  {side_path}"
            f"\nClose the locker and run:"
            f"\n  move /Y {side_path.name} {CSV_PATH.name}",
            file=sys.stderr,
        )
        return 2
    print(f"\nWrote {CSV_PATH}  ({len(df)} rows x {len(df.columns)} columns)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
