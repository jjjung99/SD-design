"""
Top 50 SD (6bp) selection from final_strict.csv

What this script does
- Reads final_strict.csv (already passed: strong + ortho, violations=0)
- Builds a ranking score per (SD, spacer) row
- Aggregates to SD-level by taking the BEST spacer per SD (max score)
- Outputs:
  1) top50_SD_best_spacer.csv  : Top 50 SD, with the best spacer and metrics
  2) top50_rows.csv            : Top 50 (SD, spacer) rows by score (sanity)
  3) sd_summary_all.csv        : SD-level summary table (optional)

Ranking (default, tunable)
- Prefer larger cross-binding energies (more positive = less WT cross-binding):
    cross_sum = dg_WT_ASDcore__O_SD + dg_WT_SD__O_ASDcore
- Prefer weaker secondary structure near SD region:
    mfe_penalty = abs(hairpin_MFE_context)  (closer to 0 is better)
    extra penalty if MFE is too stable (e.g., < -5)
- Prefer strong ΔG centered (optional):
    dg_center_bonus = -abs(dg_O_SD__O_ASD_full10 - target)
"""

import os
import pandas as pd

# =========================
# Input / Output
# =========================
IN_CSV = r"C:\Users\rhtmd\orthogonal_sd_results_v2\final_strict.csv"
   # <- change if needed (e.g., r"orthogonal_sd_results_v2\final_strict.csv")
OUT_DIR = r"top50_sd_outputs"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_TOP50_SD = os.path.join(OUT_DIR, "top50_SD_best_spacer.csv")
OUT_TOP50_ROWS = os.path.join(OUT_DIR, "top50_rows.csv")
OUT_SD_SUMMARY = os.path.join(OUT_DIR, "sd_summary_all.csv")

# =========================
# Column names (must match your CSV)
# =========================
COL_SD      = "SD_DNA"
COL_SPACER  = "Spacer_DNA"

COL_STRONG  = "dg_O_SD__O_ASD_full10"     # strong criterion already satisfied in final_strict
COL_CROSS1  = "dg_WT_ASDcore__O_SD"       # WT ASD core (CTCCTT) vs O-SD
COL_CROSS2  = "dg_WT_SD__O_ASDcore"       # WT SD (AGGAGG) vs O-ASD core
COL_MFE     = "hairpin_MFE_context"

# =========================
# Tunable ranking weights
# =========================
W_CROSS = 1.0        # weight for cross_sum
W_MFE   = 1.0        # weight for abs(MFE) penalty
W_DG    = 0.4        # weight for DG closeness-to-target bonus (optional)

DG_TARGET = -9.25    # center of [-10, -8.5]
MFE_STABLE_THRESHOLD = -5.0  # below this, structure might be too stable -> extra penalty
MFE_EXTRA_PENALTY = 2.0      # added penalty when MFE < threshold

TOP_N_SD = 50
TOP_N_ROWS = 50

# =========================
# Load
# =========================
df = pd.read_csv(IN_CSV)

# Basic sanity checks
required_cols = [COL_SD, COL_SPACER, COL_STRONG, COL_CROSS1, COL_CROSS2, COL_MFE]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in {IN_CSV}: {missing}")

# =========================
# Build ranking score per row
# =========================
df = df.copy()

# Cross-binding should be large positive (less chance of WT cross-binding)
df["cross_sum"] = df[COL_CROSS1] + df[COL_CROSS2]

# Hairpin: closer to 0 is better (avoid strong structure)
df["mfe_abs"] = df[COL_MFE].abs()

# Extra penalty if MFE too stable (very negative)
df["mfe_extra"] = (df[COL_MFE] < MFE_STABLE_THRESHOLD).astype(int) * MFE_EXTRA_PENALTY

# Strong ΔG preference: closer to target is better (optional)
df["dg_dev"] = (df[COL_STRONG] - DG_TARGET).abs()

# Final score (higher is better)
# score = +cross_sum  - abs(MFE) - extraStablePenalty - DG_deviation*weight
df["score"] = (
    W_CROSS * df["cross_sum"]
    - W_MFE * df["mfe_abs"]
    - df["mfe_extra"]
    - W_DG * df["dg_dev"]
)

# =========================
# 1) Top rows (SD+spacer) for sanity
# =========================
top_rows = df.sort_values("score", ascending=False).head(TOP_N_ROWS)
top_rows.to_csv(OUT_TOP50_ROWS, index=False)

# =========================
# 2) SD-level aggregation: pick BEST spacer per SD
# =========================
# For each SD, keep the row with maximum score
idx_best = df.groupby(COL_SD)["score"].idxmax()
best_per_sd = df.loc[idx_best].copy()

# Optional: SD-level summary metrics (how robust across spacers?)
sd_summary = (
    df.groupby(COL_SD)
      .agg(
          n_spacers=(COL_SPACER, "nunique"),
          score_max=("score", "max"),
          score_mean=("score", "mean"),
          score_std=("score", "std"),
          cross_sum_max=("cross_sum", "max"),
          cross_sum_mean=("cross_sum", "mean"),
          mfe_abs_min=("mfe_abs", "min"),
          dg_full10_mean=(COL_STRONG, "mean"),
      )
      .reset_index()
      .sort_values("score_max", ascending=False)
)
sd_summary.to_csv(OUT_SD_SUMMARY, index=False)

# Top 50 SD by BEST spacer score
top50_sd = best_per_sd.sort_values("score", ascending=False).head(TOP_N_SD)

# Keep a clean set of columns for the report
report_cols = [
    COL_SD, COL_SPACER,
    COL_STRONG, COL_CROSS1, COL_CROSS2, "cross_sum",
    COL_MFE, "mfe_abs", "mfe_extra",
    "dg_dev",
    "score"
]
top50_sd[report_cols].to_csv(OUT_TOP50_SD, index=False)

print(f"[DONE] Loaded rows: {len(df)}")
print(f"[DONE] Unique SDs: {df[COL_SD].nunique()}")
print(f"[SAVED] {OUT_TOP50_SD}")
print(f"[SAVED] {OUT_TOP50_ROWS}")
print(f"[SAVED] {OUT_SD_SUMMARY}")

# Print quick preview
print("\nTop 10 SD (with best spacer):")
print(top50_sd[report_cols].head(10).to_string(index=False))
