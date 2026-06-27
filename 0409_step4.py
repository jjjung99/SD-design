import csv
import random
import RNA
import pandas as pd

# =========================================================
# 1. 파일 경로
# =========================================================
INPUT_CSV = r"C:\Users\Minnie\Desktop\OSD_design\step3_pass_candidates_L9.csv"
CDS_FASTA = r"C:\Users\Minnie\Desktop\OSD_design\data\cds_from_genomic.fna"

OUTPUT_CSV = r"C:\Users\Minnie\Desktop\OSD_design\step4_final_candidates.csv"
TOP20_CSV = r"C:\Users\Minnie\Desktop\OSD_design\top20_candidates.csv"

# =========================================================
# 2. 설정값
# =========================================================
UPSTREAM_LEN = 30
DOWNSTREAM_LEN = 30
N_CDS_SAMPLE = 30   # 속도 vs 정확도 밸런스

ACCESSIBILITY_RATIO_MIN = 0.5

# =========================================================
# 3. 고정 서열
# =========================================================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"
AU_RICH_DNA = "TTAATTAA"
AUG_DNA = "ATG"

# =========================================================
# 4. 유틸
# =========================================================
def dna_to_rna(seq):
    return seq.replace("T", "U")

def calc_mfe(seq):
    structure, mfe = RNA.fold(dna_to_rna(seq))
    return structure, mfe

# =========================================================
# CDS 로드
# =========================================================
def load_cds(path):
    cds_list = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip().upper()
                if len(seq) >= DOWNSTREAM_LEN:
                    cds_list.append(seq)
    return cds_list

# =========================================================
# RBS window 
# ATG 기준 -30 ~ +30
# =========================================================
def build_rbs_window(sd, spacer, cds):

    upstream_full = UTR_DNA + AU_RICH_DNA + sd + spacer

    # ATG 기준 -30
    upstream = upstream_full[-UPSTREAM_LEN:]

    # ATG + CDS 앞 30nt
    downstream = AUG_DNA + cds[:DOWNSTREAM_LEN]

    return upstream + downstream

# =========================================================
# SD 접근성
# =========================================================
def check_sd_accessibility(structure, sd_len):
    sd_start = UPSTREAM_LEN - sd_len
    sd_struct = structure[sd_start:sd_start + sd_len]

    paired = sd_struct.count("(") + sd_struct.count(")")
    return paired < (sd_len / 2)

# =========================================================
# 후보 로드
# =========================================================
def load_candidates(path):
    rows = []
    with open(path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

# =========================================================
# Step4 평가
# =========================================================
def evaluate_candidate(row, cds_pool):

    sd = row.get("sd_dna") or row.get("SD")
    spacer = row.get("spacer_dna") or row.get("Spacer")

    # CDS 샘플링
    sampled_cds = random.sample(cds_pool, N_CDS_SAMPLE)

    mfe_list = []
    accessible_count = 0

    for cds in sampled_cds:

        seq = build_rbs_window(sd, spacer, cds)
        structure, mfe = calc_mfe(seq)

        mfe_list.append(mfe)

        if check_sd_accessibility(structure, len(sd)):
            accessible_count += 1

    dg_avg = sum(mfe_list) / len(mfe_list)
    accessibility_ratio = accessible_count / len(mfe_list)

    # =====================================================
    # Step3 변수 유지 + Step4 추가
    # =====================================================
    result = row.copy()

    dg_bind = float(row["dg_bind"])
    dg_struct = float(row["dg_struct"])

    result["dg_avg_cds"] = round(dg_avg, 3)
    result["dg_total"] = round(dg_struct + dg_avg, 3)
    result["binding_minus_struct"] = round(dg_bind - dg_struct, 3)
    result["Rank"] = round(dg_bind - dg_struct - dg_avg, 3)

    result["accessibility_ratio"] = round(accessibility_ratio, 3)

    # 최종 필터
    result["final_pass"] = accessibility_ratio >= ACCESSIBILITY_RATIO_MIN

    # ranking score
    result["score"] = round(accessibility_ratio - abs(dg_avg), 3)

    return result

# =========================================================
# 실행
# =========================================================
def main():

    print("[1] Step3 결과 로드")
    candidates = load_candidates(INPUT_CSV)
    print(f"후보 수: {len(candidates)}")

    print("[2] CDS 로드")
    cds_pool = load_cds(CDS_FASTA)
    print(f"CDS 수: {len(cds_pool)}")

    print("[3] Step4 평가 시작")

    results = []

    for i, row in enumerate(candidates, 1):

        result = evaluate_candidate(row, cds_pool)
        results.append(result)

        if i % 100 == 0 or i == len(candidates):
            print(f"진행률: {i}/{len(candidates)}")

    print("[4] 결과 저장")

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print("[5] Top20 추출")

    df = pd.DataFrame(results)

    top20 = df.sort_values("score", ascending=False).head(20)
    top20.to_csv(TOP20_CSV, index=False)

    print("\n===== TOP 20 =====")
    print(top20[["SD", "Spacer", "score", "accessibility_ratio"]])

    passed = sum(r["final_pass"] for r in results)

    print("\n===== 완료 =====")
    print(f"전체 후보: {len(results)}")
    print(f"최종 통과: {passed}")
    print(f"Top20 저장: {TOP20_CSV}")

# =========================================================
# Run
# =========================================================
if __name__ == "__main__":
    main()