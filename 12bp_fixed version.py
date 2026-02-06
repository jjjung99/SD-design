import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count
import RNA

# ==================================================
# Fixed sequences
# ==================================================

UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"
AU_RICH_DNA = "TTAATTAA"

WT_SD_DNA = "AGGAGG"
WT_ASD_CORE_DNA = "CTCCTT"
WT_ASD_12_RNA = "AUCACCUCCUUA"

UPSTREAM_LEN = 15

# ==================================================
# Design space
# ==================================================

SD_LEN = 6
SPACER_LEN = 6

GC_MIN = 2
GC_MAX = 4
SPACER_BASES = ["A", "T"]

# ==================================================
# Energy thresholds
# ==================================================

DG_STRONG_MIN = -10.0
DG_STRONG_MAX = -8.5
ORTHO_CUTOFF = 0.0

# ==================================================
# Output
# ==================================================

OUT_DIR = "csv_B_extended_softcore"
os.makedirs(OUT_DIR, exist_ok=True)

FINAL_CSV = os.path.join(OUT_DIR, "csv_B_candidates.csv")

# ==================================================
# Utilities
# ==================================================

def dna_to_rna(seq):
    return seq.replace("T", "U")

def revcomp_dna(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

def duplex_dG(rna_a, rna_b):
    return float(RNA.duplexfold(rna_a, rna_b).energy)

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
# Core evaluation (CSV_B logic)
# ==================================================

def score_one(args):
    sd_dna, spacer_dna = args

    # RNA forms
    sd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)
    sdsp12_rna = sd_rna + spacer_rna

    # O-ASD definitions
    oasd12_dna = revcomp_dna(sd_dna + spacer_dna)[:12]
    oasd12_rna = dna_to_rna(oasd12_dna)

    oasd_core_dna = revcomp_dna(sd_dna)
    oasd_core_rna = dna_to_rna(oasd_core_dna)

    # Context
    up15_rna = dna_to_rna(upstream_15_dna())
    context_rna = up15_rna + sdsp12_rna

    # =========================
    # A. Functional (HARD)
    # =========================
    dg_A = duplex_dG(oasd12_rna, sdsp12_rna)
    A_ok = (DG_STRONG_MIN <= dg_A <= DG_STRONG_MAX)

    # =========================
    # C3. WT-ASD12 vs SD+Spacer (HARD)
    # =========================
    dg_C3 = duplex_dG(WT_ASD_12_RNA, sdsp12_rna)
    C3_ok = (dg_C3 > ORTHO_CUTOFF)

    # =========================
    # C4. WT-ASD12 vs UP15+SD+Spacer (HARD)
    # =========================
    dg_C4 = duplex_dG(WT_ASD_12_RNA, context_rna)
    C4_ok = (dg_C4 > ORTHO_CUTOFF)

    # =========================
    # Soft evaluation only
    # =========================

    # B0. WT-ASD12 vs O-SD6
    dg_B0 = duplex_dG(WT_ASD_12_RNA, sd_rna)

    # B2. WT-SD vs O-ASD(core)
    dg_B2 = duplex_dG(dna_to_rna(WT_SD_DNA), oasd_core_rna)

    final_ok = A_ok and C3_ok and C4_ok

    return {
        "SD_DNA": sd_dna,
        "Spacer_DNA": spacer_dna,

        # Hard filter energies
        "dg_OASD12__SDSp12": dg_A,
        "dg_WT_ASD12__SDSp12": dg_C3,
        "dg_WT_ASD12__UP15_SDSp12": dg_C4,

        # Soft evaluation
        "dg_WT_ASD12__O_SD6": dg_B0,
        "dg_WT_SD__O_ASDcore": dg_B2,

        # Flags
        "A_ok": int(A_ok),
        "C3_ok": int(C3_ok),
        "C4_ok": int(C4_ok),
        "final_ok": int(final_ok),
    }

# ==================================================
# Main
# ==================================================

def main():
    combos = [
        (sd, sp)
        for sd in generate_sd_list()
        for sp in generate_spacer_list()
    ]

    print(f"[INFO] Total combinations: {len(combos)}")

    nproc = max(1, cpu_count() - 1)
    final_rows = []

    with Pool(processes=nproc) as pool:
        for row in pool.imap_unordered(score_one, combos, chunksize=800):
            if row["final_ok"] == 1:
                final_rows.append(row)

    with open(FINAL_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=final_rows[0].keys())
        w.writeheader()
        w.writerows(final_rows)

    print(f"[RESULT] 12bp_final candidates: {len(final_rows)}")
    print(f"[SAVED] {FINAL_CSV}")

if __name__ == "__main__":
    main()
