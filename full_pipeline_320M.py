import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count
import RNA


# ------------------------------
# FIXED SEQUENCES
# ------------------------------
UPSTREAM = "TTAATTAA"   # fixed 8bp AU-rich
WT_ASD_CORE = "CTCCTT"
WT_SD = "AGGAGG"

# DNA→RNA
def dna_to_rna(x): return x.replace("T","U")
WT_ASD_CORE_RNA = dna_to_rna(WT_ASD_CORE)
WT_SD_RNA = dna_to_rna(WT_SD)


# ------------------------------
# ΔG 계산
# ------------------------------
def duplex_dG(a, b):
    return RNA.duplexfold(a, b).energy


# ------------------------------
# O-ASD 생성
# ------------------------------
def revcomp_dna(seq):
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(comp[b] for b in reversed(seq))

def make_O_ASD(sd_dna):
    core = revcomp_dna(sd_dna)             # 6bp core
    full10 = "AA" + core + "AA"            # professor-style
    return core, full10


# ------------------------------
# SD 후보 생성 (GC=2~4)
# ------------------------------
def generate_SD():
    bases = ["A","T","C","G"]
    out = []
    for t in product(bases, repeat=6):
        s = "".join(t)
        gc = sum(1 for b in s if b in ("G","C"))
        if 2 <= gc <= 4:
            out.append(s)
    return out  # ~1020개


# ------------------------------
# Spacer 후보 생성 (AT-only)
# ------------------------------
def generate_spacer():
    bases = ["A","T"]
    return ["".join(t) for t in product(bases, repeat=6)]  # 64개


# ------------------------------
# CDS window 로딩 (sliding 6bp)
# ------------------------------
def load_cds_windows(dna, k=6):
    dna = dna.replace("\n","").replace(" ","")
    return [ dna[i:i+k] for i in range(len(dna)-k+1) ]


# ------------------------------
# Step1: evaluate UTR20
# ------------------------------
def eval_UTR20(args):
    sd_dna, sp_dna = args

    oasd_core, oasd_full = make_O_ASD(sd_dna)

    osd_rna = dna_to_rna(sd_dna)
    oasd_core_rna = dna_to_rna(oasd_core)
    oasd_full_rna = dna_to_rna(oasd_full)

    # 1) strong binding
    dg1 = duplex_dG(osd_rna, oasd_full_rna)
    strong_ok = (-10 <= dg1 <= -8.5)

    # 2) WT-ASD_core vs O-SD
    dg2 = duplex_dG(osd_rna, WT_ASD_CORE_RNA)
    ortho1_ok = (dg2 > 0)

    # 3) WT-SD vs O-ASD_core
    dg3 = duplex_dG(WT_SD_RNA, oasd_core_rna)
    ortho2_ok = (dg3 > 0)

    return {
        "SD": sd_dna,
        "Spacer": sp_dna,
        "OASD_core": oasd_core,
        "OASD_full": oasd_full,
        "dg_osd_oasd_full": dg1,
        "dg_wt_asdcore_osd": dg2,
        "dg_wt_sd_oasdcore": dg3,
        "strong_ok": strong_ok,
        "ortho_ok": (strong_ok and ortho1_ok and ortho2_ok),
    }


# ------------------------------
# Step2: count harmful CDS bindings
# ------------------------------
WT_ASD_CORE_RNA = dna_to_rna(WT_ASD_CORE)

def eval_crossbinding(args):
    sd_dna, win = args
    osd_rna = dna_to_rna(sd_dna)
    win_rna = dna_to_rna(win)
    dg = duplex_dG(osd_rna, win_rna)
    return (dg <= -1.0)   # harmful 기준 (교수님 threshold)


# ------------------------------
# MAIN PIPELINE
# ------------------------------
def run_full_pipeline(cds_dna):

    # --- make combinations ---
    SDs = generate_SD()          # ~1020
    SPs = generate_spacer()      # 64
    combos = [(sd, sp) for sd in SDs for sp in SPs]
    print(f"[STEP1] UTR20 combinations = {len(combos)} (expected ~65280)")

    # multiprocessing cores
    NPROC = max(1, cpu_count()-1)
    print(f"Using CPU cores = {NPROC}")

    # --- Step1: evaluate UTR20 ---
    with Pool(NPROC) as pool:
        step1 = list(pool.imap_unordered(eval_UTR20, combos, chunksize=500))

    # filter strict
    strict = [r for r in step1 if r["ortho_ok"]]
    print(f"[STEP1] strict SD candidates = {len(strict)}")

    if len(strict) == 0:
        print("No SD passed strict orthogonality.")
        return

    # -------------------------
    # Step2: CDS windows
    # -------------------------
    cds_windows = load_cds_windows(cds_dna)
    print(f"[STEP2] CDS windows = {len(cds_windows)}")

    # 전체 해밀턴량 = strict SD count × CDS windows
    total_checks = len(strict) * len(cds_windows)
    print(f"[STEP2] total ΔG evaluations = {total_checks:,}")

    # harmful windows 계산
    final_rows = []

    with Pool(NPROC) as pool:
        for r in strict:
            sd = r["SD"]
            # pack jobs
            jobs = [(sd, w) for w in cds_windows]
            harmful = sum(pool.imap_unordered(eval_crossbinding, jobs, chunksize=5000))
            r["harmful_windows"] = harmful
            final_rows.append(r)

    # -------------------------
    # SAVE CSV
    # -------------------------
    os.makedirs("results_320M", exist_ok=True)
    out = "results_320M/final_candidates.csv"

    keys = final_rows[0].keys()
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(final_rows)

    print(f"[DONE] saved final candidates → {out}")


# ------------------------------
# ENTRY POINT
# ------------------------------
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

    CDS = CDS.replace("\n", "").replace(" ", "")
    run_full_pipeline(CDS)
