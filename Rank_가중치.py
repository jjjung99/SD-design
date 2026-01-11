import os
import numpy as np
import pandas as pd

# ============================================================
# INPUT / OUTPUT
# ============================================================
INPUT_CSV = r"C:\Users\rhtmd\orthogonal_sd_results_ALL_only\final_strict.csv"

OUT_DIR = os.path.dirname(INPUT_CSV)
RANKED_ALL_CSV = os.path.join(OUT_DIR, "ranked_all_weighted_B.csv")
TOP50_CSV = os.path.join(OUT_DIR, "top50_weighted_B.csv")

# ============================================================
# Weights (사용자 확정: 균형형 B)  -- 합=1.00
# ============================================================
W_UP15   = 0.35  # dg_WT_ASD12__UP15_O_SDSp27 (bigger better)
W_CORE1  = 0.12  # dg_WT_ASDcore__O_SD (bigger better)
W_CORE2  = 0.12  # dg_WT_SD__O_ASDcore (bigger better)
W_WT12   = 0.08  # dg_WT_ASD12__O_SDSp12 (bigger better)
W_FUNC   = 0.20  # dg_O_SD__O_ASD_core6 (closer to -9.25 is better)
W_STRUCT = 0.13  # hairpin_MFE_context (closer to 0 is better)

TARGET_FUNC = -9.25  # midpoint of [-10, -8.5]

# ============================================================
# Helpers
# ============================================================
def minmax(arr: np.ndarray) -> np.ndarray:
    """Min-max normalize to [0,1]. If constant -> zeros."""
    mn = np.nanmin(arr)
    mx = np.nanmax(arr)
    if mx - mn == 0:
        return np.zeros_like(arr)
    return (arr - mn) / (mx - mn)

def higher_is_better(series: pd.Series) -> np.ndarray:
    """Normalize so higher is better -> [0,1]."""
    return minmax(series.to_numpy(dtype=float))

def closer_is_better(series: pd.Series, target: float) -> np.ndarray:
    """
    Normalize so closer to target is better -> [0,1].
    score = 1 - minmax(|x-target|)
    """
    dist = (series.to_numpy(dtype=float) - target)
    dist = np.abs(dist)
    return 1.0 - minmax(dist)

# ============================================================
# Load
# ============================================================
df = pd.read_csv(INPUT_CSV).copy()

required_cols = [
    "SD_DNA", "Spacer_DNA",
    "dg_WT_ASD12__UP15_O_SDSp27",
    "dg_WT_ASDcore__O_SD",
    "dg_WT_SD__O_ASDcore",
    "dg_WT_ASD12__O_SDSp12",
    "dg_O_SD__O_ASD_core6",
    "hairpin_MFE_context",
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns in CSV: {missing}")

# Deduplicate by (SD, Spacer)
df["key"] = df["SD_DNA"].astype(str) + "_" + df["Spacer_DNA"].astype(str)
before = len(df)
df = df.drop_duplicates(subset=["key"]).reset_index(drop=True)
after = len(df)
print(f"[INFO] Loaded: {before} rows -> Unique (SD,Spacer): {after} rows")

# ============================================================
# Component scores (0~1, higher is better)
# ============================================================
df["score_up15"]  = higher_is_better(df["dg_WT_ASD12__UP15_O_SDSp27"])
df["score_core1"] = higher_is_better(df["dg_WT_ASDcore__O_SD"])
df["score_core2"] = higher_is_better(df["dg_WT_SD__O_ASDcore"])
df["score_wt12"]  = higher_is_better(df["dg_WT_ASD12__O_SDSp12"])

# Functionality stability: closer to -9.25 is better
df["score_func_-9.25"] = closer_is_better(df["dg_O_SD__O_ASD_core6"], target=TARGET_FUNC)

# Structure: closer to 0 is better (abs MFE)
df["abs_hairpin_MFE"] = df["hairpin_MFE_context"].astype(float).abs()
df["score_struct_0"] = closer_is_better(df["abs_hairpin_MFE"], target=0.0)

# ============================================================
# Final Score (weighted sum)
# ============================================================
df["Score"] = (
    W_UP15   * df["score_up15"] +
    W_CORE1  * df["score_core1"] +
    W_CORE2  * df["score_core2"] +
    W_WT12   * df["score_wt12"] +
    W_FUNC   * df["score_func_-9.25"] +
    W_STRUCT * df["score_struct_0"]
)

# ============================================================
# Rank + Save
# ============================================================
df_ranked = df.sort_values(by="Score", ascending=False).reset_index(drop=True)
df_ranked.insert(0, "Rank", df_ranked.index + 1)

# 보기 좋은 컬럼 순서
front_cols = [
    "Rank", "SD_DNA", "Spacer_DNA", "Score",
    "dg_WT_ASD12__UP15_O_SDSp27",
    "dg_WT_ASDcore__O_SD",
    "dg_WT_SD__O_ASDcore",
    "dg_WT_ASD12__O_SDSp12",
    "dg_O_SD__O_ASD_core6",
    "hairpin_MFE_context",
    "abs_hairpin_MFE",
    "score_up15", "score_core1", "score_core2", "score_wt12",
    "score_func_-9.25", "score_struct_0",
]
final_cols = front_cols + [c for c in df_ranked.columns if c not in front_cols]
df_ranked = df_ranked[final_cols]

df_ranked.to_csv(RANKED_ALL_CSV, index=False, encoding="utf-8-sig")
df_ranked.head(50).to_csv(TOP50_CSV, index=False, encoding="utf-8-sig")

print("[SAVED]")
print(" -", RANKED_ALL_CSV)
print(" -", TOP50_CSV)

print("\n[TOP5 Preview]")
print(df_ranked[["Rank","SD_DNA","Spacer_DNA","Score"]].head(5))
