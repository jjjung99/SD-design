import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count
import RNA

# ==================================================
# Fixed sequences (불변)
# ==================================================

UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"
AU_RICH_DNA = "TTAATTAA"

WT_SD_DNA = "AGGAGG"
WT_ASD_CORE_DNA = "CTCCTT"
WT_ASD_12_RNA = "AUCACCUCCUUA"

UPSTREAM_LEN = 15

# ==================================================
# Design space (설계 대상)
# ==================================================

SD_LEN = 6
SPACER_LEN = 6

SPACER_BASES = ["A", "T"]  # spacer는 AT-only

# ==================================================
# Energy thresholds (교수님 조건 그대로)
# ==================================================

DG_STRONG_MIN = -10.0
DG_STRONG_MAX = -8.5
ORTHO_CUTOFF = 0.0

# ==================================================
# Output
# ==================================================

OUT_DIR = "final_SD_screening_with_None"
os.makedirs(OUT_DIR, exist_ok=True)

ALL_CSV = os.path.join(OUT_DIR, "all_results.csv")
FINAL_CSV = os.path.join(OUT_DIR, "final_candidates.csv")
RANKED_CSV = os.path.join(OUT_DIR, "ranked_candidates.csv")

# ==================================================
# Utilities
# ==================================================

def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")

def revcomp_dna(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

def duplex_dG_preserve(rna_a: str, rna_b: str) -> float:
    """
    ViennaRNA duplex ΔG
    결합이 형성되지 않으면 큰 양수값(1000.0)으로 보존
    """
    dG = float(RNA.duplexfold(rna_a, rna_b).energy)
    if dG > 1000:
        return 1000.0
    return dG

def upstream_15_dna():
    return (UTR_DNA + AU_RICH_DNA)[-UPSTREAM_LEN:]

# ==================================================
# Candidate generation
# ==================================================

def generate_sd_list_random():
    bases = ["A", "T", "C", "G"]
    for tup in product(bases, repeat=SD_LEN):
        yield "".join(tup)

def generate_spacer_list():
    for tup in product(SPACER_BASES, repeat=SPACER_LEN):
        yield "".join(tup)

# ==================================================
# Core evaluation
# ==================================================

def score_one(args):
    sd_dna, spacer_dna = args

    sd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)
    sdsp12_rna = sd_rna + spacer_rna

    oasd_core_rna = dna_to_rna(revcomp_dna(sd_dna))
    up15_rna = dna_to_rna(upstream_15_dna())
    context_rna = up15_rna + sdsp12_rna

    # --------------------------------------------------
    # A. Translation-competent (기능성, 수치 필수)
    # --------------------------------------------------
    dg_A = duplex_dG_preserve(sd_rna, oasd_core_rna)
    A_ok = (DG_STRONG_MIN <= dg_A <= DG_STRONG_MAX)

    # --------------------------------------------------
    # B. Core orthogonality (결합 여부만)
    # --------------------------------------------------
    dg_B1 = duplex_dG_preserve(dna_to_rna(WT_ASD_CORE_DNA), sd_rna)
    B1_ok = (dg_B1 > ORTHO_CUTOFF)

    dg_B2 = duplex_dG_preserve(dna_to_rna(WT_SD_DNA), oasd_core_rna)
    B2_ok = (dg_B2 > ORTHO_CUTOFF)

    # --------------------------------------------------
    # C. Extended orthogonality
    # --------------------------------------------------
    dg_C3 = duplex_dG_preserve(WT_ASD_12_RNA, sdsp12_rna)
    C3_ok = (dg_C3 > ORTHO_CUTOFF)

    dg_C4 = duplex_dG_preserve(WT_ASD_12_RNA, context_rna)
    C4_ok = (dg_C4 > ORTHO_CUTOFF)

    final_ok = A_ok and B1_ok and B2_ok and C3_ok and C4_ok

    return {
        "SD_DNA": sd_dna,
        "Spacer_DNA": spacer_dna,

        "dg_O_SD__O_ASDcore": dg_A,
        "dg_WT_ASDcore__O_SD": dg_B1,
        "dg_WT_SD__O_ASDcore": dg_B2,
        "dg_WT_ASD12__SDSp12": dg_C3,
        "dg_WT_ASD12__UP15_SDSp12": dg_C4,

        "A_ok": int(A_ok),
        "B1_ok": int(B1_ok),
        "B2_ok": int(B2_ok),
        "C3_ok": int(C3_ok),
        "C4_ok": int(C4_ok),

        "final_ok": int(final_ok),
    }

# ==================================================
# Ranking (Lexicographic)
# ==================================================

def rank_key(row):
    return (
        abs(row["dg_O_SD__O_ASDcore"] + 9.0),   # 기능성 최우선
        row["dg_WT_ASD12__UP15_SDSp12"],        # context 직교성
        row["dg_WT_ASD12__SDSp12"],             # intrinsic 직교성
    )

# ==================================================
# Main
# ==================================================

def main():
    sds = list(generate_sd_list_random())
    spacers = list(generate_spacer_list())
    combos = [(sd, sp) for sd in sds for sp in spacers]

    print(f"[INFO] SD candidates: {len(sds)}")
    print(f"[INFO] Spacer candidates: {len(spacers)}")
    print(f"[INFO] Total combinations: {len(combos)}")

    nproc = max(1, cpu_count() - 1)

    all_rows = []
    final_rows = []

    with Pool(processes=nproc) as pool:
        for row in pool.imap_unordered(score_one, combos, chunksize=800):
            all_rows.append(row)
            if row["final_ok"] == 1:
                final_rows.append(row)

    # Save ALL
    with open(ALL_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_rows[0].keys())
        w.writeheader()
        w.writerows(all_rows)

    # Save FINAL
    with open(FINAL_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=final_rows[0].keys())
        w.writeheader()
        w.writerows(final_rows)

    # Ranking
    ranked = sorted(final_rows, key=rank_key)

    with open(RANKED_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ranked[0].keys())
        w.writeheader()
        w.writerows(ranked)

    print(f"[RESULT] ALL rows: {len(all_rows)}")
    print(f"[RESULT] FINAL candidates: {len(final_rows)}")
    print(f"[SAVED] {ALL_CSV}")
    print(f"[SAVED] {FINAL_CSV}")
    print(f"[SAVED] {RANKED_CSV}")

if __name__ == "__main__":
    main()






