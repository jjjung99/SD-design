import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple

import RNA  # ViennaRNA python bindings
import matplotlib.pyplot as plt


# ============================================================
# 1) 고정값 (불변)  — docx 기준
# ============================================================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"  # 50bp fixed
AU_RICH_DNA = "TTAATTAA"  # 8bp fixed
AUG_DNA = "ATG"  # start codon fixed

EMPTY_ORF_PREFIX = (
    "CTGCTGGGTGAGCTTTCTCCGTAAACTTAAAGGAAAAGATTCCGTTGAAAGATT"
    "CAAAGCTATCGTTCAGCGTATACAAGAGACTTCCTCCTGAGACTCGTGTTCCC"
    "GTACCGAACTCT"
)

WT_ASD_12_RNA = "AUCACCUCCUUA"  # fixed for orthogonality check
WT_SD_RNA = "AGGAGG"            # WT SD 6bp (reverse crosstalk check)


# ============================================================
# 2) 변하는 값 (탐색 대상) — docx 기준
# ============================================================
SD_LEN = 6
SPACER_LEN = 6

# SD: GC count 2~4
SD_GC_MIN = 2
SD_GC_MAX = 4

# Spacer: DNA A/T only
SPACER_DNA_ALPHABET = ("A", "T")


# ============================================================
# 3) 결합성/직교성 조건 — 사용자가 확정한 (1)~(4) 그대로
# ============================================================
# (1) 결합성: ΔG(O-ASD12 ↔ SD+spacer(12)) ∈ [-10, -8.5]
DG_BIND_MIN = -10.0
DG_BIND_MAX = -8.5

# (2)(3)(4) 직교성: ΔG > 0 (반드시 양수)
DG_ORTHO_MIN = 0.0


# ============================================================
# 4) 저장 경로
# ============================================================
OUT_DIR = "2월 3일 3차시도 결과"  # 네 로컬 경로로 바꾸려면 여기만 수정
os.makedirs(OUT_DIR, exist_ok=True)

ALL_RESULTS_CSV = os.path.join(OUT_DIR, "all_results.csv")
FILTERED_CSV = os.path.join(OUT_DIR, "filtered_candidates.csv")

PLOT_DG_BIND = os.path.join(OUT_DIR, "plot_dg_binding.png")
PLOT_MFE = os.path.join(OUT_DIR, "plot_mfe.png")
PLOT_ORTHO_VIOL = os.path.join(OUT_DIR, "plot_ortho_violations.png")


# ============================================================
# 5) 유틸
# ============================================================
def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")


def gc_count_dna(seq: str) -> int:
    return sum(1 for b in seq if b in ("G", "C"))


def is_valid_sd(sd_dna: str) -> bool:
    g = gc_count_dna(sd_dna)
    return SD_GC_MIN <= g <= SD_GC_MAX


def duplex_dg_rna(a_rna: str, b_rna: str) -> float:
    d = RNA.duplexfold(a_rna, b_rna)
    return float(d.energy)


def fold_mfe_rna(seq_rna: str) -> float:
    _, mfe = RNA.fold(seq_rna)
    return float(mfe)


def make_o_asd12_from_sdsp12(sdsp12_rna: str) -> str:
    """
    O-ASD12 = reverse-complement( SD+spacer 12 ) in RNA space
    """
    comp = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(sdsp12_rna))


# ============================================================
# 6) UP15 정의 (핵심 수정): (UTR + AU_RICH)의 마지막 15bp
# ============================================================
UP15_DNA = (UTR_DNA + AU_RICH_DNA)[-15:]


# ============================================================
# 7) 결과 Row 정의
# ============================================================
@dataclass
class Row:
    sd_dna: str
    spacer_dna: str

    sd_rna: str
    spacer_rna: str
    sdsp12_rna: str
    o_asd12_rna: str

    # (1) 결합성
    dg_oasd12__sdsp12: float

    # (2)(3)(4) 직교성
    dg_wt_asd12__osd6: float
    dg_wt_asd12__up15_sdsp27: float
    dg_oasd12__wt_sd6: float

    # 분포/추가 분석용
    mfe_sdsp12: float
    mfe_up15_sdsp27: float

    # pass flags
    pass_binding: bool
    pass_ortho2: bool
    pass_ortho3: bool
    pass_ortho4: bool
    pass_all: bool


# ============================================================
# 8) 단일 조합 평가 (계산 후 flag 기록, 생성 단계에서 필터링 X)
# ============================================================
def evaluate_one(sd_dna: str, spacer_dna: str) -> Dict:
    sd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)

    sdsp12_rna = sd_rna + spacer_rna
    o_asd12_rna = make_o_asd12_from_sdsp12(sdsp12_rna)

    # (1) 결합성
    dg_bind = duplex_dg_rna(o_asd12_rna, sdsp12_rna)

    # (2) 직교성
    dg_ortho2 = duplex_dg_rna(WT_ASD_12_RNA, sd_rna)

    # (3) 직교성: UP15 + SD + spacer = 27bp
    up15_sdsp27_rna = dna_to_rna(UP15_DNA + sd_dna + spacer_dna)
    dg_ortho3 = duplex_dg_rna(WT_ASD_12_RNA, up15_sdsp27_rna)

    # (4) 직교성 (reverse crosstalk)
    dg_ortho4 = duplex_dg_rna(o_asd12_rna, WT_SD_RNA)

    # hairpin(MFE) 기록 (요구사항: 분포 제공)
    mfe_sdsp12 = fold_mfe_rna(sdsp12_rna)
    mfe_up15_sdsp27 = fold_mfe_rna(up15_sdsp27_rna)

    pass_binding = (DG_BIND_MIN <= dg_bind <= DG_BIND_MAX)
    pass_2 = (dg_ortho2 > DG_ORTHO_MIN)
    pass_3 = (dg_ortho3 > DG_ORTHO_MIN)
    pass_4 = (dg_ortho4 > DG_ORTHO_MIN)
    pass_all = pass_binding and pass_2 and pass_3 and pass_4

    row = Row(
        sd_dna=sd_dna,
        spacer_dna=spacer_dna,
        sd_rna=sd_rna,
        spacer_rna=spacer_rna,
        sdsp12_rna=sdsp12_rna,
        o_asd12_rna=o_asd12_rna,
        dg_oasd12__sdsp12=dg_bind,
        dg_wt_asd12__osd6=dg_ortho2,
        dg_wt_asd12__up15_sdsp27=dg_ortho3,
        dg_oasd12__wt_sd6=dg_ortho4,
        mfe_sdsp12=mfe_sdsp12,
        mfe_up15_sdsp27=mfe_up15_sdsp27,
        pass_binding=pass_binding,
        pass_ortho2=pass_2,
        pass_ortho3=pass_3,
        pass_ortho4=pass_4,
        pass_all=pass_all,
    )
    return asdict(row)


def worker(task: Tuple[str, str]) -> Dict:
    sd_dna, spacer_dna = task
    return evaluate_one(sd_dna, spacer_dna)


# ============================================================
# 9) 전체 조합 생성 (규칙 기반 생성으로 중복 제거/효율화)
#    - SD: 4^6 중 GC 2~4만
#    - spacer: 2^6 (A/T only)
# ============================================================
def generate_tasks() -> List[Tuple[str, str]]:
    sd_alphabet = ("A", "C", "G", "T")
    spacer_alphabet = SPACER_DNA_ALPHABET

    sds = []
    for tup in product(sd_alphabet, repeat=SD_LEN):
        sd = "".join(tup)
        if is_valid_sd(sd):
            sds.append(sd)

    spacers = ["".join(t) for t in product(spacer_alphabet, repeat=SPACER_LEN)]

    return [(sd, sp) for sd in sds for sp in spacers]


# ============================================================
# 10) CSV 저장
# ============================================================
def save_csv(path: str, rows: List[Dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ============================================================
# 11) 시각화 (요구사항: ΔG 분포, hairpin(MFE) 분포, orthogonality 위반 분포)
# ============================================================
def plot_distributions(rows: List[Dict]) -> None:
    dg_bind = [r["dg_oasd12__sdsp12"] for r in rows]
    mfe_12 = [r["mfe_sdsp12"] for r in rows]
    mfe_27 = [r["mfe_up15_sdsp27"] for r in rows]

    # (A) ΔG 결합성 분포
    plt.figure()
    plt.hist(dg_bind, bins=80, alpha=0.85)
    plt.axvline(DG_BIND_MIN)
    plt.axvline(DG_BIND_MAX)
    plt.title("ΔG Binding: O-ASD12 ↔ (SD+spacer)12")
    plt.xlabel("kcal/mol")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(PLOT_DG_BIND, dpi=200)
    plt.close()

    # (B) MFE 분포
    plt.figure()
    plt.hist(mfe_12, bins=80, alpha=0.7, label="MFE(SD+spacer 12)")
    plt.hist(mfe_27, bins=80, alpha=0.7, label="MFE(UP15+SD+spacer 27)")
    plt.title("MFE (hairpin) distributions")
    plt.xlabel("kcal/mol")
    plt.ylabel("count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_MFE, dpi=200)
    plt.close()

    # (C) 직교성 위반 분포(카운트)
    viol_2 = sum(1 for r in rows if not r["pass_ortho2"])
    viol_3 = sum(1 for r in rows if not r["pass_ortho3"])
    viol_4 = sum(1 for r in rows if not r["pass_ortho4"])

    plt.figure()
    labels = [
        "viol(2) WT-ASD12↔O-SD6 (ΔG<=0)",
        "viol(3) WT-ASD12↔UP15+SD+spacer (ΔG<=0)",
        "viol(4) O-ASD12↔WT-SD6 (ΔG<=0)",
    ]
    counts = [viol_2, viol_3, viol_4]
    plt.bar(range(len(counts)), counts)
    plt.xticks(range(len(counts)), labels, rotation=15, ha="right")
    plt.title("Orthogonality violations count")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(PLOT_ORTHO_VIOL, dpi=200)
    plt.close()


# ============================================================
# 12) Main
# ============================================================
def main():
    tasks = generate_tasks()
    print(f"UP15_DNA (len={len(UP15_DNA)}): {UP15_DNA}")
    print(f"Total tasks (SD GC2~4 x spacer AT-only): {len(tasks):,}")

    nproc = max(1, cpu_count() - 1)
    print(f"Using processes: {nproc}")

    rows: List[Dict] = []
    with Pool(processes=nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(worker, tasks, chunksize=500), start=1):
            rows.append(r)
            if i % 20000 == 0:
                print(f"  computed {i:,}/{len(tasks):,}")

    # 1) 전체 계산 결과 저장
    save_csv(ALL_RESULTS_CSV, rows)
    print(f"Saved all results: {ALL_RESULTS_CSV}")

    # 2) 필터 통과 후보 저장
    filtered = [r for r in rows if r["pass_all"]]
    if filtered:
        save_csv(FILTERED_CSV, filtered)
        print(f"Saved filtered candidates: {FILTERED_CSV} (n={len(filtered):,})")
    else:
        print("Filtered candidates: n=0 (no rows passed all filters).")

    # 3) 시각화 저장
    plot_distributions(rows)
    print("Saved plots:")
    print(f" - {PLOT_DG_BIND}")
    print(f" - {PLOT_MFE}")
    print(f" - {PLOT_ORTHO_VIOL}")


if __name__ == "__main__":
    main()
