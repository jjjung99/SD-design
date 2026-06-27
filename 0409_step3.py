import RNA
import pandas as pd
from multiprocessing import Pool, cpu_count

# ==================================================
# Parameters
# ==================================================
STRUCT_MODE = "tight"   # "tight" (-4) or "loose" (-6)
STRUCT_CUTOFF = -4 if STRUCT_MODE == "tight" else -6

NPROC = max(1, cpu_count()-1)

INPUT_CSV = "ortho_pass_candidates_L9.csv"
OUT_CSV = "step3_pass_candidates_L9.csv"

# ==================================================
# Fixed sequences
# ==================================================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"
AU_RICH_DNA = "TTAATTAA"
AUG_DNA = "ATG"

# ==================================================
# Utils
# ==================================================
def dna_to_rna(seq):
    return seq.replace("T","U")

def fold_dg(seq):
    try:
        return RNA.fold(seq)[1]
    except:
        return 1000

# ==================================================
# Worker 
# ==================================================
def step3_worker(row):

    # 기존 정보 유지
    result = row.copy()

    sd = row["SD"]
    spacer = row["Spacer"]

    # =========================
    # 구조 계산 (CDS 없음)
    # =========================
    seq = dna_to_rna(
        UTR_DNA + AU_RICH_DNA + sd + spacer + AUG_DNA
    )

    dg_struct = fold_dg(seq)

    # =========================
    # 필터
    # =========================
    if dg_struct < STRUCT_CUTOFF:
        return None

    # 결과 추가
    result["dg_struct"] = dg_struct

    return result

# ==================================================
# Main
# ==================================================
def main():

    print("[INFO] Loading orthogonality-passed candidates...")

    df = pd.read_csv(INPUT_CSV)
    rows = df.to_dict("records")

    print(f"[INFO] Input candidates: {len(rows)}")

    results = []

    with Pool(NPROC) as pool:
        for r in pool.imap_unordered(step3_worker, rows, chunksize=100):
            if r:
                results.append(r)

    print(f"[RESULT] Step3 pass: {len(results)}")

    df_out = pd.DataFrame(results)

    # 구조 좋은 순으로 정렬
    df_out = df_out.sort_values("dg_struct", ascending=False)

    df_out.to_csv(OUT_CSV, index=False)

    print(f"[SAVED] {OUT_CSV}")
    print(df_out.head(10))

# ==================================================
# Run
# ==================================================
if __name__ == "__main__":
    main()