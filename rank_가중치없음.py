import os
import numpy as np
import pandas as pd

# ============================================================
# INPUT / OUTPUT
# ============================================================
INPUT_CSV = r"C:\Users\rhtmd\orthogonal_sd_results_ALL_only\final_strict.csv"

OUT_DIR = os.path.dirname(INPUT_CSV)
RANKED_ALL_CSV = os.path.join(OUT_DIR, "ranked_all_lexicographic.csv")
TOP50_CSV = os.path.join(OUT_DIR, "top50_lexicographic.csv")

TARGET_STRONG = -9.25  # midpoint of [-10, -8.5]

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
# Derived tie-break columns (NO weights)
# ============================================================
# core margin: 둘 중 더 약한 쪽(min)을 기준으로 "최악의 경우에도 직교성 확보"를 우선
df["core_margin_min"] = np.minimum(
    df["dg_WT_ASDcore__O_SD"].astype(float),
    df["dg_WT_SD__O_ASDcore"].astype(float)
)

# strong stability: -9.25에 가까울수록 좋음 (거리 작을수록)
df["strong_dist_to_-9.25"] = (df["dg_O_SD__O_ASD_core6"].astype(float) - TARGET_STRONG).abs()

# structure: 0에 가까울수록 좋음 (절댓값 작을수록)
df["abs_hairpin_MFE"] = df["hairpin_MFE_context"].astype(float).abs()

# ============================================================
# Lexicographic ranking (priority order)
# ============================================================
# 1) UP15 cross-binding 억제: bigger better
# 2) core 직교성: bigger better (min margin)
# 3) WT-ASD12 vs SDSp12: bigger better
# 4) strong 안정성: smaller distance better
# 5) 구조: smaller abs better
df_ranked = df.sort_values(
    by=[
        "dg_WT_ASD12__UP15_O_SDSp27",
        "core_margin_min",
        "dg_WT_ASD12__O_SDSp12",
        "strong_dist_to_-9.25",
        "abs_hairpin_MFE",
        # 마지막 안전장치(완전 동률 방지): key
        "key",
    ],
    ascending=[
        False,  # 1
        False,  # 2
        False,  # 3
        True,   # 4
        True,   # 5
        True,   # key
    ]
).reset_index(drop=True)

df_ranked.insert(0, "Rank", df_ranked.index + 1)

# 보기 좋은 컬럼 순서
front_cols = [
    "Rank", "SD_DNA", "Spacer_DNA",
    "dg_WT_ASD12__UP15_O_SDSp27",
    "dg_WT_ASDcore__O_SD",
    "dg_WT_SD__O_ASDcore",
    "core_margin_min",
    "dg_WT_ASD12__O_SDSp12",
    "dg_O_SD__O_ASD_core6",
    "strong_dist_to_-9.25",
    "hairpin_MFE_context",
    "abs_hairpin_MFE",
]
final_cols = front_cols + [c for c in df_ranked.columns if c not in front_cols]

df_ranked = df_ranked[final_cols]

# Save
df_ranked.to_csv(RANKED_ALL_CSV, index=False, encoding="utf-8-sig")
df_ranked.head(50).to_csv(TOP50_CSV, index=False, encoding="utf-8-sig")

print("[SAVED]")
print(" -", RANKED_ALL_CSV)
print(" -", TOP50_CSV)

print("\n[TOP5 Preview]")
print(df_ranked[["Rank", "SD_DNA", "Spacer_DNA",
                 "dg_WT_ASD12__UP15_O_SDSp27", "core_margin_min",
                 "dg_WT_ASD12__O_SDSp12", "strong_dist_to_-9.25",
                 "abs_hairpin_MFE"]].head(5))
