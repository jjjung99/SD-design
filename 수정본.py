# ============================================================
# Orthogonal SD design (SD 6bp is the ONLY design target)
# - Fixed: UTR(50), AU-rich(8 fixed), AUG, empty ORF prefix, coding DNA
# - Variable: SD(6, GC count 2~4), spacer(6, AT-only)
# - Key idea:
#   * WT-ASD core is CTCCTT (core motif) -> 중심으로 cross-binding(orthogonality) 평가
#   * Do NOT optimize "SD 6bp alone energy" as a standalone objective.
#   * Instead, evaluate:
#       (1) translation-competent binding: ΔG(O-SD : O-ASD) in [-10, -8.5]
#       (2) orthogonality vs WT core: ΔG(WT_ASD_core : O-SD) > 0
#       (3) orthogonality reverse: ΔG(WT_SD : O_ASD_core) > 0
#       (4) hairpin MFE (context folding) for structure avoidance
# - Compute ALL combos first (no filtering during generation)
# - Save CSV + plots
# ============================================================

import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count

import RNA
import matplotlib.pyplot as plt


# ===============================
# 1) Fixed sequences (DNA; unchanged)
# ===============================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"  # 50bp fixed
AU_RICH_DNA = "TTAATTAA"  # 8bp fixed
AUG_DNA = "ATG"

EMPTY_ORF_PREFIX = (
    "CTGCTGGGTGAGCTTTCTCCGTAAACTTAAAGGAAAAGATTCCGTTGAAAGATT"
    "CAAAGCTATCGTTCAGCGTATACAAGAGACTTCCTCCTGAGACTCGTGTTCCC"
    "GTACCGAACTCT"
)

CODING_DNA = (
    "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAA"
    "TGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAA"
    "CCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCC"
    "GCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGG"
    "CATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTA"
    "TGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTAC"
    "CTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGAT"
    "GGCAGGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGC"
    "TGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTG"
    "CCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCGCTAAAGTTCTTCACCC"
    "CCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATACCGGAAATCCTCAAGCACCAGGTACGC"
    "TCATTGGTGCCAGCCGTGATGAAGACGAATTACCGGTCAAGGGCATTTCCAATCTGAATAACATGGCAATGTTCAGCGTT"
    "TCTGGTCCGGGGATGAAAGGGATGGTCGGCATGGCGGCGCGCGTCTTTGCAGCGATGTCACGCGCCCGTATTTCCGTGGT"
    "GCTGATTACGCAATCATCTTCCGAATACAGCATCAGTTTCTGCGTTCCACAAAGCGACTGTGTGCGAGCTGAACGGGCAA"
    "TGCAGGAAGAGTTCTACCTGGAACTGAAAGAAGGCTTACTGGAGCCGCTGGCAGTGACGGAACGGCTGGCCATTATCTCG"
    "GTGGTAGGTGATGGTATGCGCACCTTGCGTGGGATCTCGGCGAAATTCTTTGCCGCACTGGCCCGCGCCAATATCAACAT"
    "TGTCGCCATTGCTCAGGGATCTTCTGAACGCTCAATCTCTGTCGTGGTAAATAACGATGATGCGACCACTGGCGTGCGCG"
    "TTACTCATCAGATGCTGTTCAATACCGATCAGGTTATCGAAGTGTTTGTGATTGGCGTCGGTGGCGTTGGCGGTGCGCTG"
    "CTGGAGCAACTGAAGCGTCAGCAAAGCTGGCTGAAGAATAAACATATCGACTTACGTGTCTGCGGTGTTGCCAACTCGAA"
    "GGCTCTGCTCACCAATGTACATGGCCTTAATCTGGAAAACTGGCAGGAAGAACTGGCGCAAGCCAAAGAGCCGTTTAATC"
    "TCGGGCGCTTAATTCGCCTCGTGAAAGAATATCATCTGCTGAACCCGGTCATTGTTGACTGCACTTCCAGCCAGGCAGTG"
    "GCGGATCAATATGCCGACTTCCTGCGCGAAGGTTTCCACGTTGTCACGCCGAACAAAAAGGCCAACACCTCGTCGATGGA"
    "TTACTACCATCAGTTGCGTTATGCGGCGGAAAAATCGCGGCGTAAATTCCTCTATGACACCAACGTTGGGGCTGGATTAC"
    "CGGTTATTGAGAACCTGCAAAATCTGCTCAATGCAGGTGATGAATTGATGAAGTTCTCCGGCATTCTTTCTGGTTCGCTT"
    "TCTTATATCTTCGGCAAGTTAGACGAAGGCATGAGTTTCTCCGAGGCGACCACGCTGGCGCGGGAAATGGGTTATACCGA"
    "ACCGGACCCGCGAGATGATCTTTCTGGTATGGATGTGGCGCGTAAACTATTGATTCTCGCTCGTGAAACGGGACGTGAAC"
    "TGGAGCTGGCGGATATTGAAATTGAACCTGTGCTGCCCGCAGAGTTTAACGCCGAGGGTGATGTTGCCGCTTTTATGGCG"
    "AATCTGTCACAACTCGACGATCTCTTTGCCGCGCGCGTGGCGAAGGCCCGTGATGAAGGAAAAGTTTTGCGCTATGTTGG"
    "CAATATTGATGAAGATGGCGTCTGCCGCGTGAAGATTGCCGAAGTGGATGGTAATGATCCGCTGTTCAAAGTGAAAAATG"
    "GCGAAAACGCCCTGGCCTTCTATAGCCACTATTATCAGCCGCTGCCGTTGGTACTGCGCGGATATGGTGCGGGCAATGAC"
    "GTTACAGCTGCCGGTGTCTTTGCTGATCTGCTACGTACCCTCTCATGGAAGTTAGGAGTCTGA"
)

# WT motifs (DNA)
WT_SD_DNA = "AGGAGG"
WT_ASD_CORE_DNA = "CTCCTT"
WT_ASD_FULL_DNA = "CACCTCCTTA"


# ===============================
# 2) Variable design rules
# ===============================
SD_LEN = 6
SPACER_LEN = 6


# ===============================
# 3) Criteria
# ===============================
DG_STRONG_MIN = -10.0
DG_STRONG_MAX = -8.5
ORTHO_POSITIVE_CUTOFF = 0.0

# context folding length (avoid too heavy compute)
FOLD_ORF_PREFIX_N = 60


# ===============================
# Output (Windows-safe)
# ===============================
OUT_DIR = "orthogonal_sd_results_v2"
os.makedirs(OUT_DIR, exist_ok=True)

ALL_CSV = os.path.join(OUT_DIR, "all_combinations.csv")
FINAL_STRICT_CSV = os.path.join(OUT_DIR, "final_strict.csv")          # strong + ortho all pass
FINAL_RELAXED_CSV = os.path.join(OUT_DIR, "final_relaxed_viol_le_1.csv")  # strong + viol<=1

PLOT_STRONG = os.path.join(OUT_DIR, "hist_strong_dG.png")
PLOT_CROSS1 = os.path.join(OUT_DIR, "hist_cross_WT_ASDcore_vs_O_SD.png")
PLOT_CROSS2 = os.path.join(OUT_DIR, "hist_cross_WT_SD_vs_O_ASDcore.png")
PLOT_MFE = os.path.join(OUT_DIR, "hist_hairpin_MFE.png")
PLOT_VIOL = os.path.join(OUT_DIR, "bar_orthogonality_violations.png")


# ===============================
# Utilities
# ===============================
def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")

def revcomp_dna(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

def duplex_dG(rna_a: str, rna_b: str) -> float:
    d = RNA.duplexfold(rna_a, rna_b)
    return float(d.energy)

def fold_mfe(rna_seq: str) -> float:
    _, mfe = RNA.fold(rna_seq)
    return float(mfe)

def build_context_rna(sd_dna: str, spacer_dna: str) -> str:
    ctx_dna = (
        UTR_DNA
        + AU_RICH_DNA
        + sd_dna
        + spacer_dna
        + AUG_DNA
        + EMPTY_ORF_PREFIX[:FOLD_ORF_PREFIX_N]
    )
    return dna_to_rna(ctx_dna)


# ===============================
# Candidate generation
# ===============================
def generate_sd_list():
    # SD: 6bp, GC count 2~4 allowed
    bases = ["A", "T", "C", "G"]
    sds = []
    for tup in product(bases, repeat=SD_LEN):
        sd = "".join(tup)
        gc = sum(1 for b in sd if b in ("G", "C"))
        if 2 <= gc <= 4:
            sds.append(sd)
    return sds

def generate_spacer_list():
    # Spacer: 6bp, AT-only
    bases = ["A", "T"]
    return ["".join(tup) for tup in product(bases, repeat=SPACER_LEN)]


# ===============================
# Precompute WT RNAs
# ===============================
WT_SD_RNA = dna_to_rna(WT_SD_DNA)
WT_ASD_CORE_RNA = dna_to_rna(WT_ASD_CORE_DNA)
WT_ASD_FULL_RNA = dna_to_rna(WT_ASD_FULL_DNA)


# ===============================
# Scoring one combo
# ===============================
def score_one(args):
    sd_dna, spacer_dna = args

    # IMPORTANT: Only SD is "designed/mutated"
    # O-ASD is derived from SD (complementary pairing partner)
    o_asd_core_dna = revcomp_dna(sd_dna)              # 6bp (core)
    o_asd_full_dna = "AA" + o_asd_core_dna + "AA"     # 10bp (simple context)

    osd_rna = dna_to_rna(sd_dna)
    oasd_core_rna = dna_to_rna(o_asd_core_dna)
    oasd_full_rna = dna_to_rna(o_asd_full_dna)

    # (1) Translation-competent binding (absolute range)
    # Use full10 for the absolute criterion; keep core6 as reference only
    dg_osd_oasd_full10 = duplex_dG(osd_rna, oasd_full_rna)
    dg_osd_oasd_core6  = duplex_dG(osd_rna, oasd_core_rna)

    # (2) Orthogonality focus: WT-ASD core (CTCCTT) vs designed O-SD
    dg_wt_asdcore_osd = duplex_dG(osd_rna, WT_ASD_CORE_RNA)

    # (3) Reverse cross-binding: WT-SD vs O-ASD core
    dg_wt_sd_oasdcore = duplex_dG(WT_SD_RNA, oasd_core_rna)

    # (4) Hairpin / structure in context
    mfe_ctx = fold_mfe(build_context_rna(sd_dna, spacer_dna))

    # Conditions
    strong_ok = (DG_STRONG_MIN <= dg_osd_oasd_full10 <= DG_STRONG_MAX)

    v1 = (dg_wt_asdcore_osd <= ORTHO_POSITIVE_CUTOFF)
    v2 = (dg_wt_sd_oasdcore <= ORTHO_POSITIVE_CUTOFF)
    viol_count = int(v1) + int(v2)
    ortho_ok = (viol_count == 0)

    return {
        "SD_DNA": sd_dna,
        "Spacer_DNA": spacer_dna,
        "O_ASD_core_DNA": o_asd_core_dna,
        "O_ASD_full10_DNA": o_asd_full_dna,

        # Binding energies
        "dg_O_SD__O_ASD_full10": dg_osd_oasd_full10,   # used for strong criterion
        "dg_O_SD__O_ASD_core6": dg_osd_oasd_core6,     # reference only (not standalone objective)

        # Cross-binding (core-focused)
        "dg_WT_ASDcore__O_SD": dg_wt_asdcore_osd,
        "dg_WT_SD__O_ASDcore": dg_wt_sd_oasdcore,

        # Structure
        "hairpin_MFE_context": mfe_ctx,

        # Flags
        "strong_ok_full10": int(strong_ok),
        "ortho_ok_core": int(ortho_ok),
        "violations_count": viol_count,
    }


# ===============================
# Plot helper
# ===============================
def save_hist(vals, path, title, xlabel):
    plt.figure()
    plt.hist(vals, bins=60)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


# ===============================
# Main
# ===============================
def main():
    sds = generate_sd_list()
    spacers = generate_spacer_list()
    total = len(sds) * len(spacers)

    print(f"SD candidates: {len(sds)}")
    print(f"Spacer candidates: {len(spacers)}")
    print(f"Total combos: {total}")

    combos = [(sd, sp) for sd in sds for sp in spacers]

    nproc = max(1, cpu_count() - 1)
    print(f"Using processes: {nproc}")

    # columns
    fieldnames = list(score_one((sds[0], spacers[0])).keys())

    strong_vals = []
    cross1_vals = []
    cross2_vals = []
    mfe_vals = []
    viol_hist = [0, 0, 0]

    strict_rows = []
    relaxed_rows = []

    with open(ALL_CSV, "w", newline="") as f_all:
        w_all = csv.DictWriter(f_all, fieldnames=fieldnames)
        w_all.writeheader()

        with Pool(processes=nproc) as pool:
            for row in pool.imap_unordered(score_one, combos, chunksize=800):
                w_all.writerow(row)

                strong_vals.append(row["dg_O_SD__O_ASD_full10"])
                cross1_vals.append(row["dg_WT_ASDcore__O_SD"])
                cross2_vals.append(row["dg_WT_SD__O_ASDcore"])
                mfe_vals.append(row["hairpin_MFE_context"])

                vc = row["violations_count"]
                if vc < 0:
                    vc = 0
                if vc > 2:
                    vc = 2
                viol_hist[vc] += 1

                # strict final: strong + ortho(pass all)
                if row["strong_ok_full10"] == 1 and row["ortho_ok_core"] == 1:
                    strict_rows.append(row)

                # relaxed: strong + violations <= 1
                if row["strong_ok_full10"] == 1 and row["violations_count"] <= 1:
                    relaxed_rows.append(row)

    # save finals
    with open(FINAL_STRICT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(strict_rows)

    with open(FINAL_RELAXED_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(relaxed_rows)

    print("Saved:")
    print(f"- {ALL_CSV}")
    print(f"- {FINAL_STRICT_CSV} (n={len(strict_rows)})")
    print(f"- {FINAL_RELAXED_CSV} (n={len(relaxed_rows)})")

    # plots
    save_hist(strong_vals, PLOT_STRONG, "Strong binding ΔG (O-SD vs O-ASD_full10)", "ΔG (kcal/mol)")
    save_hist(cross1_vals, PLOT_CROSS1, "Cross ΔG (WT_ASD_core=CTCCTT vs Designed O-SD)", "ΔG (kcal/mol)")
    save_hist(cross2_vals, PLOT_CROSS2, "Cross ΔG (WT_SD=AGGAGG vs O-ASD_core)", "ΔG (kcal/mol)")
    save_hist(mfe_vals, PLOT_MFE, "Hairpin MFE (context fold)", "MFE (kcal/mol)")

    plt.figure()
    plt.bar([0, 1, 2], viol_hist)
    plt.title("Orthogonality violations count (0~2)")
    plt.xlabel("violations_count")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(PLOT_VIOL, dpi=200)
    plt.close()

    print("Saved plots:")
    print(f"- {PLOT_STRONG}")
    print(f"- {PLOT_CROSS1}")
    print(f"- {PLOT_CROSS2}")
    print(f"- {PLOT_MFE}")
    print(f"- {PLOT_VIOL}")


if __name__ == "__main__":
    main()
