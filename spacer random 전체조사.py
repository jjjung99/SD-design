import os
import csv
import random
from itertools import product
from multiprocessing import Pool, cpu_count
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple

import RNA
import matplotlib.pyplot as plt


# ============================================================
# 1) 고정값 (불변)
# ============================================================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"  # 50bp
AU_RICH_DNA = "TTAATTAA"  # 8bp
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

WT_ASD_12_RNA = "AUCACCUCCUUA"  # WT-ASD 12bp
WT_SD_RNA = "AGGAGG"            # WT SD 6bp

UP15_DNA = (UTR_DNA + AU_RICH_DNA)[-15:]


# ============================================================
# 2) 변이 규칙
# ============================================================
SD_LEN = 6
SPACER_LEN = 6

SD_GC_MIN = 2
SD_GC_MAX = 4

# spacer는 이제 A/C/G/T 전부 허용
SPACER_DNA_ALPHABET = ("A", "C", "G", "T")

# spacer 생성 방식 선택:
# - "all": 전수(4^6=4096 spacer, 총 13M 조합)
SPACER_MODE = "all"

# ============================================================
# 3) 결합성/직교성 조건
# ============================================================
DG_BIND_MIN = -10.0
DG_BIND_MAX = -8.5
DG_ORTHO_MIN = 0.0  # 반드시 양수


# ============================================================
# 4) 2차구조(헤어핀) 기반 컷
# ============================================================
STRUCT_MFE_CUTOFF = -10.0
UPSTREAM_NT = 30
DOWNSTREAM_NT = 60


# ============================================================
# 5) 저장 경로/파일
# ============================================================
OUT_DIR = "2월 3일 결과_with_CODING_structure_spacerRandom"
os.makedirs(OUT_DIR, exist_ok=True)

ALL_RESULTS_CSV = os.path.join(OUT_DIR, "all_results.csv")
CORE_CANDIDATES_CSV = os.path.join(OUT_DIR, "core_candidates.csv")
FINAL_CANDIDATES_CSV = os.path.join(OUT_DIR, "final_candidates.csv")

PLOT_DG_BIND = os.path.join(OUT_DIR, "plot_dg_binding.png")
PLOT_MFE_SHORT = os.path.join(OUT_DIR, "plot_mfe_short.png")
PLOT_MFE_LOCAL = os.path.join(OUT_DIR, "plot_mfe_local_window.png")
PLOT_ORTHO_VIOL = os.path.join(OUT_DIR, "plot_ortho_violations.png")


# ============================================================
# 6) 유틸
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
    comp = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(sdsp12_rna))


def build_full_mrna_rna(sd_dna: str, spacer_dna: str) -> Tuple[str, int]:
    full_dna = (
        UTR_DNA
        + AU_RICH_DNA
        + sd_dna
        + spacer_dna
        + AUG_DNA
        + EMPTY_ORF_PREFIX
        + CODING_DNA
    )
    full_rna = dna_to_rna(full_dna)
    sd_start = len(UTR_DNA) + len(AU_RICH_DNA)
    return full_rna, sd_start


def local_window(seq: str, sd_start: int, up: int, down: int) -> str:
    start = max(0, sd_start - up)
    end = min(len(seq), sd_start + SD_LEN + SPACER_LEN + down)
    return seq[start:end]


# ============================================================
# 7) Row 정의
# ============================================================
@dataclass
class Row:
    sd_dna: str
    spacer_dna: str

    sd_rna: str
    spacer_rna: str
    sdsp12_rna: str
    o_asd12_rna: str

    dg_oasd12__sdsp12: float
    dg_wt_asd12__osd6: float
    dg_wt_asd12__up15_sdsp27: float
    dg_oasd12__wt_sd6: float

    mfe_sdsp12: float
    mfe_up15_sdsp27: float
    mfe_local_window: float

    pass_binding: bool
    pass_ortho2: bool
    pass_ortho3: bool
    pass_ortho4: bool
    pass_core: bool
    pass_structure: bool
    pass_final: bool


# ============================================================
# 8) 단일 평가
# ============================================================
def evaluate_one(sd_dna: str, spacer_dna: str) -> Dict:
    sd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)

    sdsp12_rna = sd_rna + spacer_rna
    o_asd12_rna = make_o_asd12_from_sdsp12(sdsp12_rna)

    # (1) 결합성
    dg_bind = duplex_dg_rna(o_asd12_rna, sdsp12_rna)

    # (2) 직교성: WT-ASD12 ↔ O-SD6
    dg_ortho2 = duplex_dg_rna(WT_ASD_12_RNA, sd_rna)

    # (3) 직교성: WT-ASD12 ↔ UP15+SD+spacer(27)
    up15_sdsp27_rna = dna_to_rna(UP15_DNA + sd_dna + spacer_dna)
    dg_ortho3 = duplex_dg_rna(WT_ASD_12_RNA, up15_sdsp27_rna)

    # (4) 직교성: O-ASD12 ↔ WT-SD6
    dg_ortho4 = duplex_dg_rna(o_asd12_rna, WT_SD_RNA)

    # 짧은 구간 구조(기록용)
    mfe_sdsp12 = fold_mfe_rna(sdsp12_rna)
    mfe_up15_sdsp27 = fold_mfe_rna(up15_sdsp27_rna)

    # core pass
    pass_binding = (DG_BIND_MIN <= dg_bind <= DG_BIND_MAX)
    pass_2 = (dg_ortho2 > DG_ORTHO_MIN)
    pass_3 = (dg_ortho3 > DG_ORTHO_MIN)
    pass_4 = (dg_ortho4 > DG_ORTHO_MIN)
    pass_core = pass_binding and pass_2 and pass_3 and pass_4

    # local 구조 MFE (core만)
    mfe_local = 0.0
    pass_structure = False
    if pass_core:
        full_rna, sd_start = build_full_mrna_rna(sd_dna, spacer_dna)
        win = local_window(full_rna, sd_start, UPSTREAM_NT, DOWNSTREAM_NT)
        mfe_local = fold_mfe_rna(win)
        pass_structure = (mfe_local > STRUCT_MFE_CUTOFF)

    pass_final = pass_core and pass_structure

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
        mfe_local_window=mfe_local,
        pass_binding=pass_binding,
        pass_ortho2=pass_2,
        pass_ortho3=pass_3,
        pass_ortho4=pass_4,
        pass_core=pass_core,
        pass_structure=pass_structure,
        pass_final=pass_final,
    )
    return asdict(row)


def worker(task: Tuple[str, str]) -> Dict:
    return evaluate_one(task[0], task[1])


# ============================================================
# 9) 전체 조합 생성
# ============================================================
def random_spacer(n: int) -> List[str]:
    # 중복 줄이려고 set으로 뽑고, 모자라면 더 뽑는 방식
    random.seed(RANDOM_SEED)
    s = set()
    while len(s) < n:
        sp = "".join(random.choice(SPACER_DNA_ALPHABET) for _ in range(SPACER_LEN))
        s.add(sp)
    return list(s)


def generate_tasks() -> List[Tuple[str, str]]:
    # SD 후보: 4^6 전부 중에서 GC 조건만 남김
    sd_alphabet = ("A", "C", "G", "T")
    sds = []
    for tup in product(sd_alphabet, repeat=SD_LEN):
        sd = "".join(tup)
        if is_valid_sd(sd):
            sds.append(sd)

    # spacer 후보: 전수 or 랜덤 샘플
    if SPACER_MODE == "all":
        spacers = ["".join(t) for t in product(SPACER_DNA_ALPHABET, repeat=SPACER_LEN)]
    else:
        spacers = random_spacer(N_SPACER_SAMPLES)

    # 전체 태스크
    tasks = [(sd, sp) for sd in sds for sp in spacers]
    return tasks


def save_csv(path: str, rows: List[Dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ============================================================
# 10) Plot
# ============================================================
def plot_distributions(rows: List[Dict]) -> None:
    dg_bind = [r["dg_oasd12__sdsp12"] for r in rows]
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

    mfe_12 = [r["mfe_sdsp12"] for r in rows]
    mfe_27 = [r["mfe_up15_sdsp27"] for r in rows]
    plt.figure()
    plt.hist(mfe_12, bins=80, alpha=0.7, label="MFE(SD+spacer 12)")
    plt.hist(mfe_27, bins=80, alpha=0.7, label="MFE(UP15+SD+spacer 27)")
    plt.title("Short-region MFE distributions")
    plt.xlabel("kcal/mol")
    plt.ylabel("count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_MFE_SHORT, dpi=200)
    plt.close()

    viol_2 = sum(1 for r in rows if not r["pass_ortho2"])
    viol_3 = sum(1 for r in rows if not r["pass_ortho3"])
    viol_4 = sum(1 for r in rows if not r["pass_ortho4"])
    plt.figure()
    labels = [
        "viol(2) WT-ASD12↔O-SD6 (ΔG<=0)",
        "viol(3) WT-ASD12↔UP15+SD+spacer (ΔG<=0)",
        "viol(4) O-ASD12↔WT-SD6 (ΔG<=0)",
    ]
    plt.bar(range(3), [viol_2, viol_3, viol_4])
    plt.xticks(range(3), labels, rotation=15, ha="right")
    plt.title("Orthogonality violations count")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(PLOT_ORTHO_VIOL, dpi=200)
    plt.close()

    local_mfe = [r["mfe_local_window"] for r in rows if r["pass_core"]]
    if local_mfe:
        plt.figure()
        plt.hist(local_mfe, bins=80, alpha=0.85)
        plt.axvline(STRUCT_MFE_CUTOFF)
        plt.title(f"Local-window MFE around SD (with CODING), cutoff={STRUCT_MFE_CUTOFF}")
        plt.xlabel("kcal/mol")
        plt.ylabel("count")
        plt.tight_layout()
        plt.savefig(PLOT_MFE_LOCAL, dpi=200)
        plt.close()


# ============================================================
# 11) Main
# ============================================================
def main():
    print(f"UP15_DNA (len={len(UP15_DNA)}): {UP15_DNA}")
    tasks = generate_tasks()

    print(f"SPACER_MODE={SPACER_MODE}")
    if SPACER_MODE == "random":
        print(f"N_SPACER_SAMPLES={N_SPACER_SAMPLES}, seed={RANDOM_SEED}")

    print(f"Total tasks: {len(tasks):,}")
    print(f"Structure window: upstream={UPSTREAM_NT}, downstream={DOWNSTREAM_NT}")
    print(f"Structure MFE cutoff: mfe_local_window > {STRUCT_MFE_CUTOFF}")

    nproc = max(1, cpu_count() - 1)
    print(f"Using processes: {nproc}")

    rows: List[Dict] = []
    with Pool(processes=nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(worker, tasks, chunksize=500), start=1):
            rows.append(r)
            if i % 20000 == 0:
                print(f"  computed {i:,}/{len(tasks):,}")

    save_csv(ALL_RESULTS_CSV, rows)
    print(f"Saved: {ALL_RESULTS_CSV}")

    core = [r for r in rows if r["pass_core"]]
    if core:
        save_csv(CORE_CANDIDATES_CSV, core)
        print(f"Saved: {CORE_CANDIDATES_CSV} (n={len(core):,})")
    else:
        print("Core candidates: n=0")

    final = [r for r in rows if r["pass_final"]]
    if final:
        save_csv(FINAL_CANDIDATES_CSV, final)
        print(f"Saved: {FINAL_CANDIDATES_CSV} (n={len(final):,})")
    else:
        print("Final candidates: n=0 (consider relaxing STRUCT_MFE_CUTOFF)")

    plot_distributions(rows)
    print("Saved plots:")
    print(f" - {PLOT_DG_BIND}")
    print(f" - {PLOT_MFE_SHORT}")
    print(f" - {PLOT_MFE_LOCAL}")
    print(f" - {PLOT_ORTHO_VIOL}")


if __name__ == "__main__":
    main()
