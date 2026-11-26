"""
Create a master subjective effects CSV combining 5D-ASC questionnaire data
across 4 psychedelic neuroimaging studies.

Datasets:
- London LSD: already has drug-placebo difference scores (0-1 scale)
- Berlin Ketamine: already has drug-placebo difference scores (0-1 scale)
- Basel LPM: raw per-condition scores (0-100 mean VAS), needs differencing
- Basel LAM: raw per-condition scores (sum VAS), needs differencing + normalization

Outputs:
- data/master_subjective_effects.csv (all datasets, all drugs, with demographics)
- data/basel_lpm/basel_lpm_subjective.csv (LSD-only, backward-compatible format)
- data/basel_lam/basel_lam_subjective.csv (LSD-only, backward-compatible format)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = PROJECT_ROOT / "data"

# Canonical 11 5D-ASC subscale names
SUBSCALES = [
    "unity", "spiritual", "bliss", "insight", "disembodiment",
    "impaired_control", "anxiety", "complex_imagery", "elementary_imagery",
    "audio_visual_syn", "changed_meaning",
]

# Short names used in London/Berlin CSVs
SHORT_TO_CANONICAL = {
    "UNITY": "unity",
    "SPIRIT": "spiritual",
    "BLISS": "bliss",
    "INSIGHT": "insight",
    "DISEM": "disembodiment",
    "CONT&COG": "impaired_control",
    "ANXIETY": "anxiety",
    "COMPLEX": "complex_imagery",
    "ELEMENT": "elementary_imagery",
    "AUD/VIS": "audio_visual_syn",
    "MEANING_PERC": "changed_meaning",
}

# LPM Excel column names
LPM_TO_CANONICAL = {
    "Experience_of_Unity": "unity",
    "Spiritual_Experience": "spiritual",
    "Blissful_State": "bliss",
    "Insightfulness": "insight",
    "Disembodiment": "disembodiment",
    "Impaired_Control_and_Cognition": "impaired_control",
    "Anxiety": "anxiety",
    "Complex_Imagery": "complex_imagery",
    "Elementary_Imagery": "elementary_imagery",
    "Audio_Visual_Synesthesiae": "audio_visual_syn",
    "Changed_Meaning_of_Percepts": "changed_meaning",
}

# LAM Excel column names (without condition prefix)
LAM_TO_CANONICAL = {
    "Experiance of Unity": "unity",
    "Spiritual Experience": "spiritual",
    "Blissful State": "bliss",
    "Insightfulness": "insight",
    "Disembodiment": "disembodiment",
    "Impaired Control and Coginition": "impaired_control",
    "Anxiety": "anxiety",
    "Complex Imagery": "complex_imagery",
    "Elementary Imagery": "elementary_imagery",
    "Audio- Visual Synesthesiae": "audio_visual_syn",
    "Changed Meaning of Percepts": "changed_meaning",
}

# 5D-ASC item counts per subscale (Studerus et al., 2010)
# Used to normalize LAM raw sum scores
ITEM_COUNTS = {
    "unity": 5,
    "spiritual": 3,
    "bliss": 3,
    "insight": 3,
    "disembodiment": 3,
    "impaired_control": 7,
    "anxiety": 6,
    "complex_imagery": 3,
    "elementary_imagery": 3,
    "audio_visual_syn": 3,
    "changed_meaning": 3,
}

# 3D-ASCr higher-order dimensions (Stocker et al., 2025, J Psychopharmacol).
# Factor = unweighted mean of its constituent 11-ASC subscale means.
# DisE uses the strict CFA single-indicator construct (impaired_control only):
# the paper drops Anxiety from DisE because of floor effects (42% zeros in their
# healthy-volunteer raw scores), and our pooled drug-placebo diffs are even more
# floored (50% <= 0, median 0.001). Anxiety is retained as the standalone
# `anxiety` subscale per the paper's recommendation.
ASCR_DIMENSIONS = {
    "pose": ["unity", "spiritual", "bliss", "insight", "disembodiment", "changed_meaning"],
    "dise": ["impaired_control"],
    "pere": ["complex_imagery", "elementary_imagery", "audio_visual_syn"],
}

CANONICAL_TO_SHORT = {v: k for k, v in SHORT_TO_CANONICAL.items()}

MASTER_COLUMNS = [
    "dataset", "subject_id", "subject_num", "drug",
    "sex", "age", "psychedelic_experience",
] + SUBSCALES

# Final column order written to the master CSV (3D-ASCr factors appended after subscales).
OUTPUT_COLUMNS = MASTER_COLUMNS + list(ASCR_DIMENSIONS.keys())


def add_3d_ascr_factors(df: pd.DataFrame) -> pd.DataFrame:
    """Append PosE / DisE / PerE columns (Stocker et al., 2025) to a master frame.

    Each higher-order factor is the unweighted nanmean of its constituent 11-ASC
    subscale means (which are already per-subject means on the 0-1 scale in our
    master file). Because our subscales are drug-placebo difference scores,
    mean-of-differences = difference-of-means, so equal-weight aggregation is
    valid and yields drug-placebo differences in the 3D-ASCr factor space.
    """
    for factor, subscales in ASCR_DIMENSIONS.items():
        df[factor] = df[subscales].mean(axis=1, skipna=True)
    return df


def load_london_lsd() -> pd.DataFrame:
    """Load London LSD dataset (already difference scores on 0-1 scale)."""
    path = DATA_DIR / "london_lsd" / "london_lsd_subejctive.csv"
    df = pd.read_csv(path)

    df = df.rename(columns=SHORT_TO_CANONICAL)
    df["dataset"] = "london_lsd"
    df["subject_num"] = df["Subject"].astype(int)
    df["subject_id"] = df["subject_num"].apply(lambda x: f"london_lsd_{x:02d}")
    df["drug"] = "LSD"
    df["sex"] = np.nan
    df["age"] = np.nan
    df["psychedelic_experience"] = np.nan

    return df[MASTER_COLUMNS]


def load_berlin_ketamine() -> pd.DataFrame:
    """Load Berlin Ketamine dataset (already difference scores on 0-1 scale).
    Merge demographics from participants.tsv."""
    path = DATA_DIR / "berlin_ketamine" / "berlin_ketamine_subjective.csv"
    df = pd.read_csv(path)

    df = df.rename(columns=SHORT_TO_CANONICAL)
    df["dataset"] = "berlin_ketamine"
    df["subject_num"] = df["Subject"].astype(int)
    df["subject_id"] = df["subject_num"].apply(lambda x: f"berlin_ketamine_{x:02d}")
    df["drug"] = "ketamine"

    # Merge demographics
    participants_path = DATA_DIR / "berlin_ketamine" / "participants.tsv"
    participants = pd.read_csv(participants_path, sep="\t")
    participants["subject_num"] = participants["participant_id"].str.extract(r"(\d+)").astype(int)

    df = df.merge(participants[["subject_num", "age", "sex"]], on="subject_num", how="left")
    df["sex"] = df["sex"].str.upper()
    df["psychedelic_experience"] = np.nan

    return df[MASTER_COLUMNS]


def load_basel_lpm(drugs=("LSD", "Psilocybin", "Meskalin")) -> pd.DataFrame:
    """Load Basel LPM dataset. Compute drug-placebo difference scores.

    LPM data is in long format (one row per subject-condition) with scores
    on a 0-100 mean VAS scale. Normalize to 0-1 by dividing by 100.
    """
    path = DATA_DIR / "basel_lpm" / "Demographics and 5D-ASC after 24 h LPM study.xlsx"
    df = pd.read_excel(path, sheet_name="R")

    lpm_subscale_cols = list(LPM_TO_CANONICAL.keys())

    # Get placebo scores per subject
    placebo = df[df["Condition"] == "Placebo"].set_index("Study_ID")[lpm_subscale_cols]

    rows = []
    for drug in drugs:
        drug_df = df[df["Condition"] == drug].copy()
        for _, row in drug_df.iterrows():
            sid = row["Study_ID"]
            subject_num = int(sid.split("-")[1])

            # Compute difference (drug - placebo) and normalize to 0-1
            pla_scores = placebo.loc[sid]
            diff = {}
            for col, canonical in LPM_TO_CANONICAL.items():
                diff[canonical] = (row[col] - pla_scores[col]) / 100.0

            # Drug name normalization
            drug_name_map = {"LSD": "LSD", "Psilocybin": "psilocybin", "Meskalin": "mescaline"}
            drug_name = drug_name_map[drug]

            rows.append({
                "dataset": "basel_lpm",
                "subject_id": f"basel_lpm_{subject_num:03d}",
                "subject_num": subject_num,
                "drug": drug_name,
                "sex": row["Sex"].upper() if isinstance(row["Sex"], str) else np.nan,
                "age": row["Age"],
                "psychedelic_experience": row["Lifetime_psychedelic_Experience"]
                    if pd.notna(row["Lifetime_psychedelic_Experience"]) else np.nan,
                **diff,
            })

    return pd.DataFrame(rows)[MASTER_COLUMNS]


def load_basel_lam(drugs=("LSD", "Amph", "MDMA")) -> pd.DataFrame:
    """Load Basel LAM dataset. Compute drug-placebo difference scores.

    LAM data is in wide format (condition-prefixed columns) with raw sum
    scores. Normalize by dividing by (n_items * 100) per subscale.
    """
    path = DATA_DIR / "basel_lam" / "Demographics and 5D-ASC after 11 h LAM study..xlsx"
    df = pd.read_excel(path, sheet_name="analysis")

    # Drug name mapping for output
    drug_names = {"LSD": "LSD", "Amph": "amphetamine", "MDMA": "MDMA"}

    rows = []
    for drug_prefix in drugs:
        for _, row in df.iterrows():
            pid = row["PatientID"]
            subject_num = int(pid.split("-")[1])

            # Build column names for this drug and placebo
            diff = {}
            any_nan = False
            for lam_col, canonical in LAM_TO_CANONICAL.items():
                drug_col = f"{drug_prefix}_{lam_col}"
                pla_col = f"Pla_{lam_col}"
                drug_val = row[drug_col]
                pla_val = row[pla_col]

                if pd.isna(drug_val):
                    any_nan = True
                    break

                max_score = ITEM_COUNTS[canonical] * 100
                diff[canonical] = (drug_val - pla_val) / max_score

            if any_nan:
                continue

            sex = "F" if row["Gender"] == 1 else "M"

            rows.append({
                "dataset": "basel_lam",
                "subject_id": f"basel_lam_{subject_num:03d}",
                "subject_num": subject_num,
                "drug": drug_names[drug_prefix],
                "sex": sex,
                "age": row["Age"],
                "psychedelic_experience": np.nan,
                **diff,
            })

    return pd.DataFrame(rows)[MASTER_COLUMNS]


# Monash PsiConnect ASC11 column names -> canonical subscales.
# Verified against the ASC11.json Descriptions (ds006110 phenotype):
#   PERCEPTS = "Changing Meaning of Percepts", COGNITION = "Impaired Cognition
#   and Control". Item counts match Studerus et al. 2010 (ITEM_COUNTS above).
MONASH_TO_CANONICAL = {
    "ASC11_UNITY": "unity",
    "ASC11_SPIRITUAL": "spiritual",
    "ASC11_BLISSFUL": "bliss",
    "ASC11_INSIGHTFUL": "insight",
    "ASC11_DISEMBODY": "disembodiment",
    "ASC11_COGNITION": "impaired_control",
    "ASC11_ANXIETY": "anxiety",
    "ASC11_COMPLEX": "complex_imagery",
    "ASC11_ELEMENTARY": "elementary_imagery",
    "ASC11_AUDIOVISUAL": "audio_visual_syn",
    "ASC11_PERCEPTS": "changed_meaning",
}


def load_monash_psilo() -> pd.DataFrame:
    """Load Monash PsiConnect (ds006110). Post-psilocybin ASC scores only.

    Scored ASC11 subscales (phenotype/scored/ses-02/ASC11.tsv) are
    mean-of-items on a 0-100 VAS scale (sum / n_items, per the JSON sidecar
    Operations), so normalisation is /100 — same scale convention as Basel LPM.

    DECISION (user, 2026-08-25): use the raw post-drug (ses-02) session
    scores with NO baseline subtraction. Context: the study is open-label
    pre/post with no placebo arm; the ses-01 (no-drug) ASC is missing
    entirely for 16 subjects and near floor for the rest (medians 0.3-4 on
    the 0-100 scale), so subtracting it would cost ~1/3 of the cohort for a
    <=0.05 shift in the 0-1 scores. Note the Basel loaders in this file do
    subtract the post-placebo ASC; their placebo ASC is likewise near floor.

    Demographics: the released demographics.tsv is privacy-redacted ("nd"),
    so sex is NaN; age comes from participants.tsv. The cohort is
    operationally psychedelic-naive (Novelli et al. 2026, Participants), so
    psychedelic_experience = 0.
    """
    pheno = DATA_DIR / "monash_psilo" / "derivatives" / "phenotype" / "scored"
    ses02 = pd.read_csv(pheno / "ses-02" / "ASC11.tsv", sep="\t")

    participants = pd.read_csv(DATA_DIR / "monash_psilo" / "participants.tsv", sep="\t")
    age_map = participants.set_index("participant_id")["age"]

    rows = []
    for _, row in ses02.iterrows():
        pid = row["participant_id"]           # e.g. "sub-PC001"
        diff = {}
        any_nan = False
        for col, canonical in MONASH_TO_CANONICAL.items():
            if pd.isna(row[col]):
                any_nan = True
                break
            diff[canonical] = row[col] / 100.0
        if any_nan:
            continue

        subject_num = int(pid.replace("sub-PC", ""))
        rows.append({
            "dataset": "monash_psilo",
            "subject_id": f"monash_psilo_PC{subject_num:03d}",
            "subject_num": subject_num,
            "drug": "psilocybin",
            "sex": np.nan,
            "age": age_map.get(pid, np.nan),
            "psychedelic_experience": 0,
            **diff,
        })

    return pd.DataFrame(rows)[MASTER_COLUMNS]


def generate_compat_csv(master_df: pd.DataFrame, dataset: str, drug_filter: str):
    """Generate backward-compatible per-dataset CSV (LSD-only, old column names)."""
    subset = master_df[
        (master_df["dataset"] == dataset) & (master_df["drug"] == drug_filter)
    ].copy()

    # Rename to short column names
    rename = {v: k for k, v in SHORT_TO_CANONICAL.items()}
    subset = subset.rename(columns=rename)

    # Select only Subject + 11 subscales (old format)
    short_cols = ["Subject"] + list(SHORT_TO_CANONICAL.keys())
    subset["Subject"] = subset["subject_num"]
    subset = subset[short_cols].sort_values("Subject").reset_index(drop=True)

    out_path = DATA_DIR / dataset / f"{dataset}_subjective.csv"
    subset.to_csv(out_path, index=False)
    print(f"  Wrote {out_path} ({len(subset)} rows)")


def validate(df: pd.DataFrame):
    """Print validation summary."""
    print("\n=== Validation ===")

    # Row counts per dataset/drug
    print("\nRow counts:")
    print(df.groupby(["dataset", "drug"]).size().to_string())

    # Value ranges
    print("\nSubscale value ranges:")
    for s in SUBSCALES:
        vals = df[s].dropna()
        print(f"  {s}: [{vals.min():.3f}, {vals.max():.3f}]  mean={vals.mean():.3f}")

    # 3D-ASCr higher-order factor ranges
    print("\n3D-ASCr factor value ranges:")
    for factor in ASCR_DIMENSIONS:
        if factor in df.columns:
            vals = df[factor].dropna()
            print(f"  {factor}: [{vals.min():.3f}, {vals.max():.3f}]  mean={vals.mean():.3f}")

    # Check for out-of-range values
    for s in SUBSCALES:
        outliers = df[(df[s] < -0.15) | (df[s] > 1.05)]
        if len(outliers) > 0:
            print(f"  WARNING: {len(outliers)} out-of-range values in {s}")

    # NaN audit
    nan_counts = df[SUBSCALES].isna().sum()
    if nan_counts.any():
        print(f"\nNaN counts: {nan_counts[nan_counts > 0].to_dict()}")
    else:
        print("\nNo NaN values in subscales.")

    # Demographics coverage
    print("\nDemographics coverage:")
    for col in ["sex", "age", "psychedelic_experience"]:
        coverage = df.groupby("dataset")[col].apply(lambda x: x.notna().sum())
        print(f"  {col}: {coverage.to_dict()}")

    # Duplicate check
    dupes = df.duplicated(subset=["subject_id", "drug"], keep=False)
    if dupes.any():
        print(f"\nWARNING: {dupes.sum()} duplicate subject_id+drug rows!")
    else:
        print("\nNo duplicate subject_id+drug combinations.")

    print(f"\nTotal rows: {len(df)}")


def main():
    parser = argparse.ArgumentParser(description="Create master subjective effects CSV")
    parser.add_argument("--no-compat", action="store_true",
                        help="Skip generating backward-compatible per-dataset CSVs")
    args = parser.parse_args()

    print("Loading datasets...")

    london = load_london_lsd()
    print(f"  London LSD: {len(london)} rows")

    berlin = load_berlin_ketamine()
    print(f"  Berlin Ketamine: {len(berlin)} rows")

    lpm = load_basel_lpm(drugs=("LSD", "Psilocybin", "Meskalin"))
    print(f"  Basel LPM: {len(lpm)} rows")

    lam = load_basel_lam(drugs=("LSD", "Amph", "MDMA"))
    print(f"  Basel LAM: {len(lam)} rows")

    monash = load_monash_psilo()
    print(f"  Monash PsiConnect: {len(monash)} rows")

    # Combine
    master = pd.concat([london, berlin, lpm, lam, monash], ignore_index=True)
    master = master.sort_values(["dataset", "drug", "subject_num"]).reset_index(drop=True)

    # Append 3D-ASCr higher-order factors (PosE, DisE, PerE).
    master = add_3d_ascr_factors(master)
    master = master[OUTPUT_COLUMNS]

    # Save master (write to a tmp file first so we don't clobber if it's locked)
    out_path = DATA_DIR / "master_subjective_effects.csv"
    tmp_path = out_path.with_suffix(".csv.tmp")
    master.to_csv(tmp_path, index=False)
    try:
        tmp_path.replace(out_path)
    except PermissionError:
        print(f"\nWARNING: {out_path} is locked (probably open in another app).")
        print(f"         Wrote to {tmp_path} instead — close the file and rename manually.")
    else:
        print(f"\nWrote master CSV: {out_path}")

    # Backward-compatible CSVs
    if not args.no_compat:
        print("\nGenerating backward-compatible CSVs...")
        generate_compat_csv(master, "basel_lpm", "LSD")
        generate_compat_csv(master, "basel_lam", "LSD")

    validate(master)


if __name__ == "__main__":
    main()
