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

GC_MIN = 2
GC_MAX = 4
SPACER_BASES = ["A", "T"]

# ==================================================
# Energy thresholds (필터 기준)
# ==================================================

# A. Translation-competent
DG_STRONG_MIN = -10.0
DG_STRONG_MAX = -8.5

# Orthogonality cutoff
ORTHO_CUTOFF = 0.0


# ==================================================
# Output
# ==================================================

OUT_DIR = "all_filters_applied"
os.makedirs(OUT_DIR, exist_ok=True)

ALL_CSV = os.path.join(OUT_DIR, "all_results.csv")
FINAL_CSV = os.path.join(OUT_DIR, "final_candidates.csv")


# ==================================================
# Utilities
# ==================================================

def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")

def revcomp_dna(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

def duplex_dG_safe(rna_a: str, rna_b: str):
    """
    ViennaRNA duplex ΔG
    - 결합이 형성되지 않으면 None 반환
    """
    dG = float(RNA.duplexfold(rna_a, rna_b).energy)
    if dG > 1000:   # ViennaRNA 관행적 처리
        return None
    return dG

def upstream_15_dna():
    return (UTR_DNA + AU_RICH_DNA)[-UPSTREAM_LEN:]


# ==================================================
# Candidate generation
# ==================================================

def generate_sd_list():
    bases = ["A", "T", "C", "G"]
    for tup in product(bases, repeat=SD_LEN):
        sd = "".join(tup)
        gc = sum(b in "GC" for b in sd)
        if GC_MIN <= gc <= GC_MAX:
            yield sd

def generate_spacer_list():
    for tup in product(SPACER_BASES, repeat=SPACER_LEN):
        yield "".join(tup)


# ==================================================
# Core evaluation (모든 필터 반영)
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
    # A. Translation-competent (필수)
    # --------------------------------------------------
    dg_A = duplex_dG_safe(sd_rna, oasd_core_rna)
    A_ok = (
        dg_A is not None and
        DG_STRONG_MIN <= dg_A <= DG_STRONG_MAX
    )

    # --------------------------------------------------
    # B1. WT-ASD(core) : O-SD 
    # --------------------------------------------------
    dg_B1 = duplex_dG_safe(dna_to_rna(WT_ASD_CORE_DNA), sd_rna)
    B1_ok = (dg_B1 is not None and dg_B1 > ORTHO_CUTOFF)

    # --------------------------------------------------
    # B2. WT-SD : O-ASD(core) 
    # --------------------------------------------------
    dg_B2 = duplex_dG_safe(dna_to_rna(WT_SD_DNA), oasd_core_rna)
    B2_ok = (dg_B2 is not None and dg_B2 > ORTHO_CUTOFF)

    # --------------------------------------------------
    # C3. WT-ASD12 : SD+Spacer 
    # --------------------------------------------------
    dg_C3 = duplex_dG_safe(WT_ASD_12_RNA, sdsp12_rna)
    C3_ok = (dg_C3 is not None and dg_C3 > ORTHO_CUTOFF)

    # --------------------------------------------------
    # C4. WT-ASD12 : UP15 + SD+Spacer 
    # --------------------------------------------------
    dg_C4 = duplex_dG_safe(WT_ASD_12_RNA, context_rna)
    C4_ok = (dg_C4 is not None and dg_C4 > ORTHO_CUTOFF)

    # --------------------------------------------------
    # FINAL FILTER (모든 필수 조건)
    # --------------------------------------------------
    final_ok = A_ok  and B1_ok and B2_ok and C3_ok and C4_ok


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
# Main
# ==================================================

def main():
    sds = list(generate_sd_list())
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

    print(f"[RESULT] ALL rows: {len(all_rows)}")
    print(f"[RESULT] FINAL candidates: {len(final_rows)}")
    print(f"[SAVED] {ALL_CSV}")
    print(f"[SAVED] {FINAL_CSV}")


if __name__ == "__main__":
    main()
