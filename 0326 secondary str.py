import csv
import RNA

# =========================
# 파일 경로
# =========================
INPUT_CSV = "fixed_spacer6_sdsp_out/candidates_pass_all_random.csv"
OUTPUT_CSV = "fixed_spacer6_sdsp_out/structure_evaluated_random.csv"

# =========================
# 시퀀스
# =========================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"
AU_RICH_DNA = "TTAATTAA"

EMPTY_ORF_PREFIX = (
    "CTGCTGGGTGAGCTTTCTCCGTAAACTTAAAGGAAAAGATTCCGTTGAAAGATT"
    "CAAAGCTATCGTTCAGCGTATACAAGAGACTTCCTCCTGAGACTCGTGTTCCC"
    "GTACCGAACTCT"
)

# =========================
# 유틸
# =========================
def dna_to_rna(seq):
    return seq.replace("T", "U")

# =========================
# RBS window (핵심 수정)
# -30 ~ +30 구현
# =========================
def build_rbs_window(sd, spacer):
    upstream = (UTR_DNA + AU_RICH_DNA)[-30:]

    downstream = (
        "ATG"
        + EMPTY_ORF_PREFIX[:30]   # 핵심 수정 (ATG +30)
    )

    return upstream + sd + spacer + downstream


# =========================
# 구조 계산
# =========================
def calc_structure(seq):
    rna = dna_to_rna(seq)
    structure, mfe = RNA.fold(rna)
    return structure, mfe


# =========================
# SD 접근성 판단
# =========================
def check_sd_accessibility(structure, sd):
    upstream_len = 30  # 정확히 -30 기준
    sd_start = upstream_len
    sd_len = len(sd)

    sd_struct = structure[sd_start:sd_start + sd_len]

    paired = sd_struct.count("(") + sd_struct.count(")")
    accessible = paired < (sd_len / 2)

    return accessible, sd_struct


# =========================
# 메인 실행
# =========================
def main():
    results = []

    mfe_values = []  # 분포 확인용

    with open(INPUT_CSV, "r") as f:
        reader = csv.DictReader(f)

        for row in reader:
            sd = row["sd_dna"]
            spacer = row["spacer_dna"]

            # 1. RBS window
            window = build_rbs_window(sd, spacer)

            # 2. 구조 계산
            structure, mfe = calc_structure(window)
            mfe_values.append(mfe)

            # 3. SD 접근성
            accessible, sd_struct = check_sd_accessibility(structure, sd)

            # 4. ΔG 기준 (논문 기반 optimal)
            pass_mfe = (-12 <= mfe <= -4)

            # 5. 최종 판단
            pass_all = pass_mfe and accessible

            results.append({
                "sd": sd,
                "spacer": spacer,
                "mfe": round(mfe, 3),
                "sd_structure": sd_struct,
                "sd_accessible": accessible,
                "pass_mfe": pass_mfe,
                "pass_all_structure": pass_all
            })

    # =========================
    # CSV 저장
    # =========================
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    # =========================
    # 요약 출력
    # =========================
    print(f"Saved: {OUTPUT_CSV}")
    print(f"Total: {len(results)}")
    print(f"Pass structure: {sum(r['pass_all_structure'] for r in results)}")

    print("\nMFE distribution:")
    print(f"Min: {min(mfe_values):.2f}")
    print(f"Max: {max(mfe_values):.2f}")
    print(f"Mean: {sum(mfe_values)/len(mfe_values):.2f}")


if __name__ == "__main__":
    main()