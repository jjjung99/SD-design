import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count
import RNA


# Fixed sequences

UPSTREAM = "TTAATTAA"   # fixed upstream 8bp
AUG = "ATG"

WT_SD = "AGGAGG"
WT_ASD_CORE = "CTCCTT"
WT_ASD_FULL = "CACCTCCTTA"


# CDS window generator

def load_cds_windows(cds: str, k: int = 6):
    windows = []
    for i in range(len(cds) - k):
        windows.append(cds[i:i+k])
    return windows


# Utility functions

def dna_to_rna(seq):
    return seq.replace("T", "U")

def revcomp(seq):
    table = {"A":"T", "T":"A", "C":"G", "G":"C"}
    return "".join(table[b] for b in reversed(seq))

def duplex_dG(a, b):
    return float(RNA.duplexfold(a, b).energy)

def make_O_ASD(sd_dna):
    rc = revcomp(sd_dna)
    return "AA" + rc + "AA"


# SD & Spacer generation

def generate_SD_list():
    bases = ["A","T","C","G"]
    out = []
    for tup in product(bases, repeat=6):
        sd = "".join(tup)
        gc = sum(1 for c in sd if c in ("G","C"))
        if 2 <= gc <= 4:
            out.append(sd)
    return out

def generate_spacer_list():
    return ["".join(p) for p in product(["A","T"], repeat=6)]


# Step2: SD-level scoring

def evaluate_UTR20(args):
    sd, sp = args

    o_asd = make_O_ASD(sd)
    o_asd_core = o_asd[2:8]

    osd_rna = dna_to_rna(sd)
    oasd_full_rna = dna_to_rna(o_asd)
    oasd_core_rna = dna_to_rna(o_asd_core)

    wt_asd_core_rna = dna_to_rna(WT_ASD_CORE)
    wt_sd_rna = dna_to_rna(WT_SD)

    # ΔG calculations
    dg1 = duplex_dG(osd_rna, oasd_full_rna)         # strong binding
    dg2 = duplex_dG(osd_rna, wt_asd_core_rna)       # WT-ASD core vs SD
    dg3 = duplex_dG(wt_sd_rna, oasd_core_rna)       # WT-SD vs O-ASD core


    return {
        "SD": sd,
        "Spacer": sp,
        "UTR20": UPSTREAM + sd + sp,
        "O_ASD": o_asd,
        "dg_O_SD__O_ASD": dg1,
        "dg_WT_ASDcore__O_SD": dg2,
        "dg_WT_SD__O_ASDcore": dg3
    }


# Step3: CDS cross-binding

def count_harmful_windows(sd_dna, cds_windows):
    sd_rna = dna_to_rna(sd_dna)
    count = 0
    for w in cds_windows:
        rc = revcomp(w)
        dg = duplex_dG(sd_rna, dna_to_rna(rc))
        if dg <= -4.0:
            count += 1
    return count


# Main pipeline

def run_pipeline(CDS):

    #Generate all UTR20
    SDs = generate_SD_list()
    SPs = generate_spacer_list()
    combos = [(sd, sp) for sd in SDs for sp in SPs]

    print(f"Total UTR20 = {len(combos)}")

    #SD-level orthogonality
    nproc = max(1, cpu_count() - 1)
    print(f"Using {nproc} CPUs")

    with Pool(nproc) as pool:
        sd_results = list(pool.imap_unordered(evaluate_UTR20, combos, chunksize=500))

    sd_candidates = [r for r in sd_results if r["strong_ok"] and r["ortho_ok"]]
    print(f"Orthogonal SD count = {len(sd_candidates)}")

    #CDS scanning
    cds_windows = load_cds_windows(CDS)
    print(f"CDS windows = {len(cds_windows)}")

    final_rows = []
    for row in sd_candidates:
        hc = count_harmful_windows(row["SD"], cds_windows)
        row["harmful_windows"] = hc
        final_rows.append(row)

    # Save results
    os.makedirs("results_professor_version", exist_ok=True)
    out_csv = "results_professor_version/final_candidates.csv"

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=final_rows[0].keys())
        w.writeheader()
        w.writerows(final_rows)

    print(f"Saved final candidates → {out_csv}")
    return final_rows


# Entry point

if __name__ == "__main__":

    CDS = (
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

    CDS = CDS.replace("\n","").replace(" ","")

    run_pipeline(CDS)
