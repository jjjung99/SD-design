import os
import csv
from multiprocessing import Pool, cpu_count
import RNA
import pandas as pd

# ==================================================
# Fixed sequences (친구 코드와 동일)
# ==================================================
WT_SD_DNA = "AGGAGG"
WT_SPACER_9_DNA = "GAAACAGCT"
WT_SPACER_FOR_SDSP12_DNA = WT_SPACER_9_DNA[:6]
WT_SDSP12_DNA = WT_SD_DNA + WT_SPACER_FOR_SDSP12_DNA

WT_ASD_6_DNA = "CTCCTT"
WT_ASD_12_DNA = "ATCACCTCCTTA"

UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"
AU_RICH_DNA = "TTAATTAA"
UP15_DNA = (UTR_DNA + AU_RICH_DNA)[-15:]

NPROC = max(1, cpu_count() - 1)

INPUT_CSV = "C:/Users/Minnie/Desktop/OSD_design/step1_select_spacer_length_binding_only_SD_random_afterSD15_GFP_out/binding_range_candidates_L9.csv"
OUT_CSV = "ortho_pass_candidates_L9.csv"

# ==================================================
# Utils 
# ==================================================
def clean_dna(seq):
    return "".join(ch for ch in seq.upper() if ch in "ATGC")

def dna_to_rna(seq):
    return clean_dna(seq).replace("T", "U")

def revcomp_dna(seq):
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    seq = clean_dna(seq)
    return "".join(comp[b] for b in reversed(seq))

def duplex_dg(a_dna, b_dna):
    try:
        a_rna = dna_to_rna(a_dna)
        b_rna = dna_to_rna(b_dna)
        return float(RNA.duplexfold(a_rna, b_rna).energy)
    except:
        return 1000

# ==================================================
# Worker 
# ==================================================
def check_orthogonality(row):

    # clean 적용
    sd = clean_dna(row["sd_dna"])
    spacer = clean_dna(row["spacer_dna"])

    sdsp = sd + spacer

    # O-system
    osdsp12 = sdsp[:12]
    oasd12 = revcomp_dna(osdsp12)
    oasd6 = revcomp_dna(sd)

    # =========================
    # 조건 1 
    # =========================
    context = UP15_DNA + sd + spacer
    dg1 = duplex_dg(WT_ASD_12_DNA, context)
    if not (dg1 > 0):
        return None

    # =========================
    # 조건 2
    # =========================
    dg2 = duplex_dg(oasd12, WT_SDSP12_DNA)
    if not (dg2 > 0):
        return None

    # =========================
    # 조건 3
    # =========================
    dg3 = duplex_dg(osdsp12, WT_ASD_12_DNA)
    if not (dg3 > 0):
        return None

    # =========================
    # 조건 4 
    # =========================
    dg4 = duplex_dg(sd, WT_ASD_6_DNA)
    if not (dg4 > 0):
        return None

    # =========================
    # 조건 5
    # =========================
    dg5 = duplex_dg(oasd6, WT_SD_DNA)
    if not (dg5 > 0):
        return None

    return {
        "SD": sd,
        "Spacer": spacer,
        "SD_spacer": sdsp,
        "oASD6": oasd6,
        "oASD12": oasd12,
        "dg_bind": row["dg_bind_oasd12__context21"],
        "pass_all": True
    }

# ==================================================
# Main
# ==================================================
def main():

    print("[INFO] Loading binding candidates...")

    df = pd.read_csv(INPUT_CSV)
    rows = df.to_dict("records")

    print(f"[INFO] Input candidates: {len(rows)}")

    results = []

    with Pool(NPROC) as pool:
        for r in pool.imap_unordered(check_orthogonality, rows, chunksize=1000):
            if r:
                results.append(r)

    print(f"[RESULT] Orthogonality pass: {len(results)}")

    df_out = pd.DataFrame(results)
    df_out.to_csv(OUT_CSV, index=False)

    print(f"[SAVED] {OUT_CSV}")
    print(df_out.head(10))

# ==================================================
# Run
# ==================================================
if __name__ == "__main__":
    main()