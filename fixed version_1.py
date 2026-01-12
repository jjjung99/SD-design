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
AUG_DNA = "ATG"

EMPTY_ORF_PREFIX = (
    "CTGCTGGGTGAGCTTTCTCCGTAAACTTAAAGGAAAAGATTCCGTTGAAAGATT"
    "CAAAGCTATCGTTCAGCGTATACAAGAGACTTCCTCCTGAGACTCGTGTTCCC"
    "GTACCGAACTCT"
)

WT_SD_DNA = "AGGAGG"
WT_ASD_CORE_DNA = "CTCCTT"
WT_ASD_12_RNA = "AUCACCUCCUUA"


# ==================================================
# Design parameters
# ==================================================

SD_LEN = 6
SPACER_LEN = 6

GC_MIN = 2
GC_MAX = 4
SPACER_BASES = ["A", "T"]

DG_STRONG_MIN = -10.0
DG_STRONG_MAX = -8.5

ORTHO_POSITIVE_CUTOFF = 0.0
UPSTREAM_LEN = 15


# ==================================================
# Output
# ==================================================

OUT_DIR = "final_sd_screening"
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

def duplex_dG_safe(rna_a: str, rna_b: str, fail_cutoff=1000):
    """
    Returns:
        float ΔG if duplex formation is valid
        None if duplex formation failed (ViennaRNA returns huge value)
    """
    dG = float(RNA.duplexfold(rna_a, rna_b).energy)
    if dG > fail_cutoff:
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
# Core evaluation
# ==================================================

def score_one(args):
    sd_dna, spacer_dna = args

    sd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)
    sdsp12_rna = sd_rna + spacer_rna

    oasd_core_rna = dna_to_rna(revcomp_dna(sd_dna))
    context_rna = dna_to_rna(upstream_15_dna()) + sdsp12_rna

    # --- A. translation functionality ---
    dg_strong = duplex_dG_safe(sd_rna, oasd_core_rna)
    strong_ok = (
        dg_strong is not None
        and DG_STRONG_MIN <= dg_strong <= DG_STRONG_MAX
    )

    # --- B. core-level cross binding (record only) ---
    dg_WT_ASDcore__O_SD = duplex_dG_safe(
        dna_to_rna(WT_ASD_CORE_DNA), sd_rna
    )
    dg_WT_SD__O_ASDcore = duplex_dG_safe(
        dna_to_rna(WT_SD_DNA), oasd_core_rna
    )

    # --- C1. SD+spacer vs WT-ASD12 ---
    dg_WT_ASD12__SDSp12 = duplex_dG_safe(
        WT_ASD_12_RNA, sdsp12_rna
    )
    sdsp12_ok = (
        dg_WT_ASD12__SDSp12 is not None
        and dg_WT_ASD12__SDSp12 > ORTHO_POSITIVE_CUTOFF
    )

    # --- C2. context safety ---
    dg_WT_ASD12__context = duplex_dG_safe(
        WT_ASD_12_RNA, context_rna
    )
    context_ok = (
        dg_WT_ASD12__context is not None
        and dg_WT_ASD12__context > ORTHO_POSITIVE_CUTOFF
    )

    # --- FINAL FILTER ---
    final_ok = strong_ok and sdsp12_ok and context_ok

    return {
        "SD_DNA": sd_dna,
        "Spacer_DNA": spacer_dna,

        "dg_O_SD__O_ASDcore": dg_strong,
        "dg_WT_ASDcore__O_SD": dg_WT_ASDcore__O_SD,
        "dg_WT_SD__O_ASDcore": dg_WT_SD__O_ASDcore,

        "dg_WT_ASD12__SDSp12": dg_WT_ASD12__SDSp12,
        "dg_WT_ASD12__UP15_SDSp12": dg_WT_ASD12__context,

        "strong_ok": int(strong_ok),
        "sdsp12_ok": int(sdsp12_ok),
        "context_ok": int(context_ok),
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