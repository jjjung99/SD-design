# orthogonal_sd_fullscan_vscode.py
# VSCode + Python 3.10 + ViennaRNA(RNA) + multiprocessing
# Full combinatorial scan (no pre-filter during generation): generate all combos then compute dG.
# Writes hits.csv + stats.csv + hist arrays, and progress logs.

import os
import csv
import time
import math
import argparse
from itertools import product
from multiprocessing import Pool, cpu_count

# ----------------------------
# ViennaRNA import (must exist)
# ----------------------------
try:
    import RNA
except Exception as e:
    raise SystemExit(
        "[FATAL] import RNA failed. ViennaRNA Python bindings are not installed in this Python.\n"
        "In VSCode terminal, run:\n"
        "  python -m pip install viennarna\n"
        "Then test:\n"
        "  python -c \"import RNA; print(RNA.__version__)\"\n\n"
        f"Original error: {e}"
    )

# ----------------------------
# Fixed sequences (DNA)
# ----------------------------
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

# 16S ASD reference (DNA; core is CTCCTT)
ASD10_DNA = "CACCTCCTTA"

# ----------------------------
# utils
# ----------------------------
def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")

def gc_percent(seq: str) -> float:
    gc = sum(1 for b in seq if b in ("G", "C"))
    return 100.0 * gc / len(seq)

def all_sd_6_with_gc_40_60():
    # 6-mer over A/C/G/T with GC% in [40, 60] => GC count must be 3 (50%) only
    # because 6 * 0.4 = 2.4 and 6 * 0.6 = 3.6, integer GC count => 3
    bases = "ACGT"
    out = []
    for tup in product(bases, repeat=6):
        s = "".join(tup)
        gc = sum(1 for c in s if c in "GC")
        if gc == 3:
            out.append(s)
    return out  # size = C(6,3)*2^3*2^3 = 20*8*8 = 1280

def all_spacer_6_AT_only():
    out = []
    for tup in product("AT", repeat=6):
        out.append("".join(tup))
    return out  # 64

def all_oasd_core_6():
    bases = "ACGT"
    return ["".join(t) for t in product(bases, repeat=6)]  # 4096

def duplex_dG_kcal(rna1: str, rna2: str) -> float:
    # ViennaRNA duplexfold returns energy in kcal/mol
    # use RNA.duplexfold for intermolecular binding
    d = RNA.duplexfold(rna1, rna2)
    return float(d.energy)

def hairpin_dG_kcal(rna: str) -> float:
    # RNA.fold returns (structure, mfe)
    _, mfe = RNA.fold(rna)
    return float(mfe)

# ----------------------------
# worker
# ----------------------------
def eval_combo(args):
    sd_dna, spacer_dna, oasd_core_dna, compute_hairpin, hairpin_window = args

    # Build O-mRNA segment around SD
    # Full construct (DNA): UTR + AU + SD + spacer + AUG + empty + coding
    # Convert to RNA for folding
    mrna_dna = (
        UTR_DNA + AU_RICH_DNA + sd_dna + spacer_dna + AUG_DNA + EMPTY_ORF_PREFIX + CODING_DNA
    )
    mrna_rna = dna_to_rna(mrna_dna)

    # Define O-ASD (use core 6 only; if you later extend to 10, do it here)
    oasd_rna = dna_to_rna(oasd_core_dna)

    # Energies:
    # 1) O-SD : O-ASD (functional strong binding target [-10,-8.5])
    dg_o = duplex_dG_kcal(dna_to_rna(sd_dna), oasd_rna)

    # 2) WT-ASD : O-SD should be > 0
    dg_wtasd_osd = duplex_dG_kcal(dna_to_rna(sd_dna), dna_to_rna(WT_ASD_CORE_DNA))

    # 3) O-ASD : WT-SD should be > 0
    dg_oasd_wtsd = duplex_dG_kcal(dna_to_rna(WT_SD_DNA), oasd_rna)

    # Hairpin/secondary structure:
    if compute_hairpin:
        if hairpin_window > 0:
            # window around SD start: take region UTR+AU+SD+spacer+AUG (+ extra)
            # Find SD start index in mrna_dna: len(UTR)+len(AU)
            sd_start = len(UTR_DNA) + len(AU_RICH_DNA)
            left = max(0, sd_start - hairpin_window)
            right = min(len(mrna_rna), sd_start + 6 + hairpin_window)
            window_rna = mrna_rna[left:right]
            hp = hairpin_dG_kcal(window_rna)
        else:
            hp = hairpin_dG_kcal(mrna_rna)
    else:
        hp = float("nan")

    return (sd_dna, spacer_dna, oasd_core_dna, dg_o, dg_wtasd_osd, dg_oasd_wtsd, hp)

# ----------------------------
# main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results", help="output directory")
    ap.add_argument("--nproc", type=int, default=max(1, cpu_count()-1), help="processes")
    ap.add_argument("--hits_only", action="store_true", help="write only passing hits")
    ap.add_argument("--dg_min", type=float, default=-10.0)
    ap.add_argument("--dg_max", type=float, default=-8.5)
    ap.add_argument("--orth_threshold", type=float, default=0.0, help="WT cross-binding must be > this")
    ap.add_argument("--compute_hairpin", action="store_true", help="also compute secondary structure mfe")
    ap.add_argument("--hairpin_window", type=int, default=60, help="window radius around SD; 0=full RNA")
    ap.add_argument("--progress_every", type=int, default=5, help="print progress every N OASD cores")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    hits_path = os.path.join(args.outdir, "hits.csv")
    stats_path = os.path.join(args.outdir, "stats.csv")

    sds = all_sd_6_with_gc_40_60()
    spacers = all_spacer_6_AT_only()
    oasds = all_oasd_core_6()

    total = len(sds) * len(spacers) * len(oasds)

    print("[INFO] SD count:", len(sds))
    print("[INFO] Spacer count:", len(spacers))
    print("[INFO] O-ASD core count:", len(oasds))
    print("[INFO] Total combos:", total)

    # streaming write
    with open(hits_path, "w", newline="", encoding="utf-8") as f_hits:
        w = csv.writer(f_hits)
        w.writerow([
            "SD_DNA", "Spacer_DNA", "O_ASD_CORE_DNA",
            "dG_O_SD_O_ASD",
            "dG_WT_ASD_O_SD",
            "dG_O_ASD_WT_SD",
            "hairpin_mfe"
        ])

        start = time.time()
        N = 0
        hits = 0
        neg_count_0 = 0  # count cases where either cross binding <=0
        neg_count_1 = 0  # count where exactly one cross binding <=0
        neg_count_2 = 0  # count where both cross bindings <=0

        # For average stats
        sum_dg_o = 0.0
        sum_hp = 0.0
        hp_N = 0

        # Iterate per OASD to show progress
        # We will use multiprocessing over SD*spacer for each OASD to limit memory
        for i, oasd in enumerate(oasds):
            tasks = (
                (sd, sp, oasd, args.compute_hairpin, args.hairpin_window)
                for sd in sds
                for sp in spacers
            )

            # chunked pool map
            with Pool(processes=args.nproc) as pool:
                for res in pool.imap_unordered(eval_combo, tasks, chunksize=500):
                    sd_dna, spacer_dna, oasd_core_dna, dg_o, dg_wtasd_osd, dg_oasd_wtsd, hp = res

                    N += 1
                    sum_dg_o += dg_o
                    if not math.isnan(hp):
                        sum_hp += hp
                        hp_N += 1

                    cross_flags = [
                        dg_wtasd_osd <= args.orth_threshold,
                        dg_oasd_wtsd <= args.orth_threshold
                    ]
                    cross_bad = sum(1 for x in cross_flags if x)
                    if cross_bad == 0:
                        pass
                    elif cross_bad == 1:
                        neg_count_1 += 1
                    else:
                        neg_count_2 += 1
                    if cross_bad > 0:
                        neg_count_0 += 1

                    # hit definition: strong O binding + both cross > threshold
                    is_hit = (args.dg_min <= dg_o <= args.dg_max) and \
                             (dg_wtasd_osd > args.orth_threshold) and \
                             (dg_oasd_wtsd > args.orth_threshold)

                    if (not args.hits_only) or is_hit:
                        w.writerow([sd_dna, spacer_dna, oasd_core_dna, dg_o, dg_wtasd_osd, dg_oasd_wtsd, hp])

                    if is_hit:
                        hits += 1

            if (i + 1) % args.progress_every == 0:
                elapsed = (time.time() - start) / 60.0
                print(f"[PROGRESS] OASD {i+1}/{len(oasds)}  evals={N:,}/{total:,}  hits={hits:,}  elapsed={elapsed:.1f} min")

        # write stats
        avg_dg = sum_dg_o / max(1, N)
        avg_hp = (sum_hp / hp_N) if hp_N > 0 else float("nan")

        with open(stats_path, "w", newline="", encoding="utf-8") as f_stat:
            w2 = csv.writer(f_stat)
            w2.writerow(["total_expected_evals", total])
            w2.writerow(["total_computed_evals", N])
            w2.writerow(["hits_written", hits])
            w2.writerow(["avg_dG_O_SD_O_ASD", avg_dg])
            w2.writerow(["avg_hairpin_mfe", avg_hp])
            w2.writerow(["crossbinding_le_threshold_any", neg_count_0])
            w2.writerow(["crossbinding_le_threshold_exactly1", neg_count_1])
            w2.writerow(["crossbinding_le_threshold_both", neg_count_2])

    print("[DONE] saved:", hits_path)
    print("[DONE] saved:", stats_path)

if __name__ == "__main__":
    main()
