import pandas as pd
import numpy as np
import warnings
import logging
import os
import sys
import yaml
from datetime import datetime
warnings.filterwarnings("ignore")

# ── Load config ───────────────────────────────────────────────────────────────
with open("config.yaml", "r") as f:
    cfg = yaml.safe_load(f)

# ── Setup logging ─────────────────────────────────────────────────────────────
os.makedirs('logs', exist_ok=True)
 
log_file = cfg["wrangling_log"]
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[
        logging.FileHandler(log_file, mode='w'),
        logging.StreamHandler(sys.stdout)
    ]
)
 
logging.info(f"Dataset Preparation Log - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
logging.info("=" * 80)

# ============================================================
# 0. LOAD DATA
# ============================================================
logging.info("\n[0] LOAD DATA")
logging.info("-" * 40)

df = pd.read_csv(cfg["updrs_path"])


df = df.rename(columns={
    "PATNO": "id",
    "EVENT_ID": "visit_name",
    "NHY": "hy",
    "PDSTATE": "state",
    'PDMEDYN': 'med_for_pd',
    'HRPOSTMED': 'post_hr',
    'NP3TOT': 'updrs3',
})

df["state"] = df["state"].str.upper().str.strip().fillna("Unknown")
df["hy"]    = pd.to_numeric(df["hy"], errors="coerce")

logging.info(f"Raw rows loaded: {len(df)}")
logging.info(f"State value counts (raw):\n{df['state'].value_counts(dropna=False).to_string()}")

# Keep only valid visits
df = df[df["visit_name"].str.startswith(("V", "BL"), na=False)].copy()
logging.info(f"Rows after visit filter (V*/BL only): {len(df)}")

# Drop rows with missing or out-of-range HY stages
n_before = len(df)
df = df[df['hy'].between(1, 5)].copy()
logging.info(f"Rows dropped for missing/out-of-range HY: {n_before - len(df)}")
logging.info(f"Rows remaining: {len(df)}")

# Merge metadata (from master key)
ms = pd.read_csv(cfg["metadata_path"], low_memory=False).drop_duplicates("GP2ID")
ms = ms[ms["study"].str.contains("PPMI")].rename(columns={"clinical_id": "id"})
ms["id"] = ms["id"].astype(int)

n_before = len(df)
df = (
    df.merge(
        ms[["id", "GP2_phenotype_for_qc", "age", "age_of_onset",
            "age_at_diagnosis", "biological_sex_for_qc"]],
        on="id", how="inner"
    )
    .rename(columns={"biological_sex_for_qc": "sex"})
)
logging.info(f"Rows after metadata merge: {len(df)} (dropped {n_before - len(df)})")

n_before = len(df)
df = df[df["GP2_phenotype_for_qc"] == "PD"].copy()
logging.info(f"Rows after PD filter: {len(df)} (dropped {n_before - len(df)})")
logging.info(f"Unique participants: {df['id'].nunique()}")
logging.info(f"Visits present: {sorted(df['visit_name'].unique())}")

df["state_upper"] = df["state"].str.upper()

# ============================================================
# 1. COMPUTE VISIT_MONTH FROM VISIT SCHEDULE
# ============================================================
logging.info("\n[1] COMPUTE VISIT_MONTH")
logging.info("-" * 40)

visit_schedule = {
    'BL': 0,   'V01': 3,   'V02': 6,   'V03': 9,   'V04': 12,
    'V05': 18,  'V06': 24,  'V07': 30,  'V08': 36,  'V09': 42,
    'V10': 48,  'V11': 54,  'V12': 60,  'V13': 72,  'V14': 84,
    'V15': 96,  'V16': 108, 'V17': 120, 'V18': 132, 'V19': 144,
    'V20': 156, 'V21': 168, 'V22': 180
}
df['visit_month'] = df['visit_name'].map(visit_schedule)
n_before = len(df)
df = df[df['visit_month'].notna()].copy()
logging.info(f"Rows dropped for unmapped visit names: {n_before - len(df)}")
logging.info(f"Visit month distribution:\n{df['visit_month'].value_counts().sort_index().to_string()}")

# ============================================================
# 2. SELECT HY PER VISIT BY CLINIC DEFINITION
# ============================================================
logging.info("\n[2] CLINIC DEFINITIONS (ON/OFF)")
logging.info("-" * 40)

n_unknown        = int((df["state_upper"] == "UNKNOWN").sum())
n_unknown_on_med = int(((df["state_upper"] == "UNKNOWN") & (df["med_for_pd"] == 1)).sum())
n_unknown_no_med_status = int(((df["state_upper"] == "UNKNOWN") & (df["med_for_pd"].isna())).sum())
logging.info(f"Unknown-state rows: {n_unknown} total -- "
             f"{n_unknown_on_med} on PD meds (excluded from fallback), "
             f"{n_unknown_no_med_status} with missing med_for_pd (excluded from fallback)")

def select_hy(group, prefer):
    preferred = group[group["state_upper"] == prefer]
    if len(preferred) > 0:
        return preferred.iloc[0]["hy"]
    # Unknown state is only usable as a stand-in for a visit before the
    # participant started PD medication -- med_for_pd == 0. A visit where
    # meds had already started, or where med_for_pd itself is missing, is
    # excluded rather than guessed at (med_for_pd NaN == 0 is False, so
    # missing values are dropped here too without a separate check).
    unknown = group[(group["state_upper"] == "UNKNOWN") & (group["med_for_pd"] == 0)]
    if len(unknown) > 0:
        return unknown.iloc[0]["hy"]
    return np.nan

on_hy = (
    df.groupby(["id", "visit_name"])
    .apply(lambda g: select_hy(g, "ON"))
    .reset_index().rename(columns={0: "hy"})
)
off_hy = (
    df.groupby(["id", "visit_name"])
    .apply(lambda g: select_hy(g, "OFF"))
    .reset_index().rename(columns={0: "hy"})
)

n_on_before, n_off_before = len(on_hy), len(off_hy)
# select_hy() returns NaN for visits with no usable HY under this clinic
# state (e.g. med_for_pd==1 with state=="Unknown", or a visit with only the
# other state's assessment recorded) -- these must be dropped here, not left
# in, since (hy >= 3).astype(int) below would silently coerce NaN to 0 and
# count an unusable visit as a confirmed HY<3 negative.
on_hy  = on_hy.dropna(subset=["hy"])
off_hy = off_hy.dropna(subset=["hy"])
logging.info(f"On Clinic  — visits: {len(on_hy)} (dropped {n_on_before - len(on_hy)} with no usable ON HY)")
logging.info(f"Off Clinic — visits: {len(off_hy)} (dropped {n_off_before - len(off_hy)} with no usable OFF HY)")

covars = (
    df.sort_values("state_upper")
    .groupby(["id", "visit_name"]).first().reset_index()
    [["id", "visit_name", "visit_month", "age", "age_of_onset",
      "age_at_diagnosis", "INFODT"]]
)
sex_col = df.groupby("id")["sex"].first().reset_index()

on_df  = on_hy.merge(covars, on=["id", "visit_name"]).merge(sex_col, on="id")
off_df = off_hy.merge(covars, on=["id", "visit_name"]).merge(sex_col, on="id")

# ============================================================
# 3. DISEASE DURATION
# ============================================================
logging.info("\n[3] DISEASE DURATION")
logging.info("-" * 40)

def add_disease_duration(df, label):
    onset_n = df["age_of_onset"].notna().sum()
    diag_n  = df["age_at_diagnosis"].notna().sum()
    ref_col = "age_of_onset" if onset_n >= diag_n else "age_at_diagnosis"
    logging.info(f"  [{label}] Using '{ref_col}' "
                 f"(age_of_onset n={onset_n}, age_at_diagnosis n={diag_n})")
    df = df.copy()
    df["disease_duration"] = df["age"] - df[ref_col]
    neg = (df["disease_duration"] < 0).sum()
    if neg > 0:
        logging.info(f"  [{label}] Dropping {neg} rows with negative disease duration")
    df = df[df["disease_duration"] >= 0]
    logging.info(f"  [{label}] Rows remaining: {len(df)}, IDs: {df['id'].nunique()}")
    return df

on_df  = add_disease_duration(on_df,  "On Clinic")
off_df = add_disease_duration(off_df, "Off Clinic")

# ============================================================
# 4. DERIVED COLUMNS
# ============================================================
logging.info("\n[4] DERIVED COLUMNS")
logging.info("-" * 40)

on_df["hy3"]       = (on_df["hy"]  >= 3).astype(int)
off_df["hy3"]      = (off_df["hy"] >= 3).astype(int)
on_df["sex_male"]  = (on_df["sex"] == "Male").astype(int)
off_df["sex_male"] = (off_df["sex"] == "Male").astype(int)
on_df["clinic"]    = "On Clinic"
off_df["clinic"]   = "Off Clinic"

for label, d in [("On Clinic", on_df), ("Off Clinic", off_df)]:
    bl = d[d["visit_month"] == 0]
    logging.info(f"  [{label}] HY3 prevalence at baseline: "
                 f"{bl['hy3'].mean():.1%} ({bl['hy3'].sum()}/{len(bl)})")
    logging.info(f"  [{label}] Sex distribution:\n"
                 f"{d.drop_duplicates('id')['sex'].value_counts().to_string()}")

# ============================================================
# 5. SAVE OUTPUTS
# ============================================================
logging.info("\n[5] SAVE OUTPUTS")
logging.info("-" * 40)

on_df.to_csv(cfg["on_clinic_path"], index=False)
off_df.to_csv(cfg["off_clinic_path"], index=False)

raw_hy = (
    df.groupby(["id", "visit_name", "visit_month", "state_upper"])[["hy", "med_for_pd", "updrs3"]]
    .first().reset_index()
    .rename(columns={"state_upper": "state"})
)
raw_hy.to_csv(cfg["raw_hy_path"], index=False)

logging.info(f"Saved on_clinic.csv  — {len(on_df)} rows, {on_df['id'].nunique()} IDs")
logging.info(f"Saved off_clinic.csv — {len(off_df)} rows, {off_df['id'].nunique()} IDs")
logging.info(f"Saved raw_hy.csv     — {len(raw_hy)} rows")
logging.info("\n" + "=" * 80)
logging.info(f"Done - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")