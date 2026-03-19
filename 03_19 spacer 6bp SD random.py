# fixed_spacer6_sdsp_candidates.py
# 조건:
# - Spacer = 6bp 고정
# - Spacer alphabet = A/T/G/C 전체
# - SD = 6bp
# - SD = A/T/G/C 전체 6bp 전수조사
# - SD + Spacer 조합 후보 생성
# - ViennaRNA duplex energy 기반 필터링

import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple

import RNA


# ============================================================
# 0) 사용자 최종 조건
# ============================================================
SD_LEN = 6
SPACER_LEN = 6

# 기존 GC 조건 변수는 남겨두되, SD 생성에는 더 이상 사용하지 않음
SD_GC_MIN = 2
SD_GC_MAX = 2

SD_ALPHABET = ("A", "T", "G", "C")
SPACER_ALPHABET = ("A", "T", "G", "C")

# binding 조건
DG_BIND_MIN = -9.0
DG_BIND_MAX = -7.5

# orthogonality 조건
DG_ORTHO_MIN = 0.0  # 반드시 > 0


# ============================================================
# 1) 고정 서열
# ============================================================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"  # 50bp
AU_RICH_DNA = "TTAATTAA"  # 8bp
AUG_DNA = "ATG"

WT_SD_DNA = "AGGAGG"
WT_SPACER_6_DNA = "GAAACAGCT"[-6:]   # 기존 WT spacer 9bp의 뒤 6bp 사용

WT_ASD_6_DNA = "CTCCTT"
WT_ASD_12_DNA = "ATCACCTCCTTA"  # 6+6 대응용

UP15_DNA = (UTR_DNA + AU_RICH_DNA)[-15:]


# ============================================================
# 2) 출력 폴더
# ============================================================
OUT_DIR = "fixed_spacer6_sdsp_out"
os.makedirs(OUT_DIR, exist_ok=True)

CANDIDATES_CSV = os.path.join(OUT_DIR, "candidates_pass_all.csv")
ALL_RESULTS_CSV = os.path.join(OUT_DIR, "all_results.csv")
SUMMARY_CSV = os.path.join(OUT_DIR, "summary.csv")


# ============================================================
# 3) 유틸
# ============================================================
def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")


def gc_count_dna(seq: str) -> int:
    return sum(1 for b in seq if b in ("G", "C"))


def is_valid_sd(sd_dna: str) -> bool:
    gc = gc_count_dna(sd_dna)
    return SD_GC_MIN <= gc <= SD_GC_MAX


def duplex_dg(a_rna: str, b_rna: str) -> float:
    d = RNA.duplexfold(a_rna, b_rna)
    return float(d.energy)


def revcomp_rna(seq_rna: str) -> str:
    comp = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq_rna))


def build_context_up15_sd_spacer(sd_dna: str, spacer_dna: str) -> str:
    return dna_to_rna(UP15_DNA + sd_dna + spacer_dna)


def save_csv(path: str, rows: List[Dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_sd_list() -> List[str]:
    # SD를 GC 조건 없이 A/T/G/C 전체 6bp 전수조사
    return ["".join(t) for t in product(SD_ALPHABET, repeat=SD_LEN)]


def make_spacer_list() -> List[str]:
    return ["".join(t) for t in product(SPACER_ALPHABET, repeat=SPACER_LEN)]


# ============================================================
# 4) 평가 row
# ============================================================
@dataclass
class Row:
    sd_dna: str
    spacer_dna: str
    osdsp_dna: str
    oasd_dna: str

    dg_bind_oasd__osdsp: float
    dg_ortho_ctx_wtasd12__up15_sdsp: float
    dg_ortho_oasd12__wt_sdsp12: float
    dg_ortho_osdsp12__wt_asd12: float
    dg_ortho_osd6__wt_asd6: float
    dg_ortho_oasd6__wt_sd6: float

    pass_binding: bool
    pass_ortho_ctx: bool
    pass_ortho_oasd_wt_sdsp: bool
    pass_ortho_osdsp_wt_asd: bool
    pass_ortho_osd6_wt_asd6: bool
    pass_ortho_oasd6_wt_sd6: bool

    pass_all: bool
    first_fail: str


def evaluate_one(task: Tuple[str, str]) -> Dict:
    sd_dna, spacer_dna = task

    sd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)

    osdsp_dna = sd_dna + spacer_dna
    osdsp_rna = sd_rna + spacer_rna

    oasd_rna = revcomp_rna(osdsp_rna)
    oasd_dna = oasd_rna.replace("U", "T")

    wt_asd6_rna = dna_to_rna(WT_ASD_6_DNA)
    wt_asd12_rna = dna_to_rna(WT_ASD_12_DNA)
    wt_sd6_rna = dna_to_rna(WT_SD_DNA)
    wt_sdsp12_rna = dna_to_rna(WT_SD_DNA + WT_SPACER_6_DNA)

    # (1) 기능성 binding
    dg_bind = duplex_dg(oasd_rna, osdsp_rna)
    pass_binding = (DG_BIND_MIN <= dg_bind <= DG_BIND_MAX)

    # (2) WT-ASD(12) vs UP15+SD+Spacer
    ctx_rna = build_context_up15_sd_spacer(sd_dna, spacer_dna)
    dg_ctx = duplex_dg(wt_asd12_rna, ctx_rna)
    pass_ctx = (dg_ctx > DG_ORTHO_MIN)

    # (3) O-ASD(12) vs WT-SDsp(12)
    dg_oasd_wt_sdsp = duplex_dg(oasd_rna, wt_sdsp12_rna)
    pass_oasd_wt_sdsp = (dg_oasd_wt_sdsp > DG_ORTHO_MIN)

    # (4) O-SDsp(12) vs WT-ASD(12)
    dg_osdsp_wt_asd = duplex_dg(osdsp_rna, wt_asd12_rna)
    pass_osdsp_wt_asd = (dg_osdsp_wt_asd > DG_ORTHO_MIN)

    # (5) O-SD(6) vs WT-ASD(6)
    dg_osd6_wt_asd6 = duplex_dg(sd_rna, wt_asd6_rna)
    pass_osd6_wt_asd6 = (dg_osd6_wt_asd6 > DG_ORTHO_MIN)

    # (6) O-ASD(6) vs WT-SD(6)
    oasd6_rna = revcomp_rna(sd_rna)
    dg_oasd6_wt_sd6 = duplex_dg(oasd6_rna, wt_sd6_rna)
    pass_oasd6_wt_sd6 = (dg_oasd6_wt_sd6 > DG_ORTHO_MIN)

    pass_all = (
        pass_binding
        and pass_ctx
        and pass_oasd_wt_sdsp
        and pass_osdsp_wt_asd
        and pass_osd6_wt_asd6
        and pass_oasd6_wt_sd6
    )

    first_fail = "PASS_ALL"
    if not pass_binding:
        first_fail = "binding"
    elif not pass_ctx:
        first_fail = "ortho_ctx"
    elif not pass_oasd_wt_sdsp:
        first_fail = "ortho_oasd_wt_sdsp"
    elif not pass_osdsp_wt_asd:
        first_fail = "ortho_osdsp_wt_asd"
    elif not pass_osd6_wt_asd6:
        first_fail = "ortho_osd6_wt_asd6"
    elif not pass_oasd6_wt_sd6:
        first_fail = "ortho_oasd6_wt_sd6"

    row = Row(
        sd_dna=sd_dna,
        spacer_dna=spacer_dna,
        osdsp_dna=osdsp_dna,
        oasd_dna=oasd_dna,

        dg_bind_oasd__osdsp=dg_bind,
        dg_ortho_ctx_wtasd12__up15_sdsp=dg_ctx,
        dg_ortho_oasd12__wt_sdsp12=dg_oasd_wt_sdsp,
        dg_ortho_osdsp12__wt_asd12=dg_osdsp_wt_asd,
        dg_ortho_osd6__wt_asd6=dg_osd6_wt_asd6,
        dg_ortho_oasd6__wt_sd6=dg_oasd6_wt_sd6,

        pass_binding=pass_binding,
        pass_ortho_ctx=pass_ctx,
        pass_ortho_oasd_wt_sdsp=pass_oasd_wt_sdsp,
        pass_ortho_osdsp_wt_asd=pass_osdsp_wt_asd,
        pass_ortho_osd6_wt_asd6=pass_osd6_wt_asd6,
        pass_ortho_oasd6_wt_sd6=pass_oasd6_wt_sd6,

        pass_all=pass_all,
        first_fail=first_fail,
    )
    return asdict(row)


# ============================================================
# 5) main
# ============================================================
def main():
    print("=== Fixed spacer 6bp / SD+Spacer candidate search ===")
    print(f"SD length: {SD_LEN}")
    print(f"Spacer length: {SPACER_LEN}")
    print("SD condition: all 6bp combinations of A/T/G/C")
    print(f"Binding condition: {DG_BIND_MIN} ~ {DG_BIND_MAX}")
    print(f"Orthogonality condition: DG > {DG_ORTHO_MIN}")
    print(f"WT spacer 6bp used: {WT_SPACER_6_DNA}")

    sds = make_sd_list()
    spacers = make_spacer_list()

    print(f"Valid SD count = {len(sds):,}")
    print(f"Spacer count   = {len(spacers):,}")
    print(f"Total pairs    = {len(sds) * len(spacers):,}")

    tasks = [(sd, sp) for sd in sds for sp in spacers]

    nproc = max(1, cpu_count() - 1)
    print(f"Using processes: {nproc}")

    all_rows = []
    pass_rows = []

    with Pool(processes=nproc) as pool:
        for idx, r in enumerate(pool.imap_unordered(evaluate_one, tasks, chunksize=1000), start=1):
            all_rows.append(r)
            if r["pass_all"]:
                pass_rows.append(r)

            if idx % 50000 == 0 or idx == len(tasks):
                print(f"{idx:,}/{len(tasks):,} done | pass_all = {len(pass_rows):,}")

    save_csv(ALL_RESULTS_CSV, all_rows)
    save_csv(CANDIDATES_CSV, pass_rows)

    unique_sd = sorted(set(r["sd_dna"] for r in pass_rows))
    unique_spacer = sorted(set(r["spacer_dna"] for r in pass_rows))

    summary = [{
        "SD_LEN": SD_LEN,
        "SPACER_LEN": SPACER_LEN,
        "SD_MODE": "all_ATGC_6bp",
        "DG_BIND_MIN": DG_BIND_MIN,
        "DG_BIND_MAX": DG_BIND_MAX,
        "DG_ORTHO_MIN": DG_ORTHO_MIN,
        "WT_SPACER_6_DNA": WT_SPACER_6_DNA,
        "n_valid_sd": len(sds),
        "n_spacers": len(spacers),
        "n_total_pairs": len(tasks),
        "n_pass_all": len(pass_rows),
        "n_unique_sd_in_pass": len(unique_sd),
        "n_unique_spacer_in_pass": len(unique_spacer),
    }]
    save_csv(SUMMARY_CSV, summary)

    print("\n=== DONE ===")
    print(f"Saved all results : {ALL_RESULTS_CSV}")
    print(f"Saved candidates  : {CANDIDATES_CSV}")
    print(f"Saved summary     : {SUMMARY_CSV}")
    print(f"Unique SD in pass : {len(unique_sd):,}")
    print(f"Unique spacer in pass : {len(unique_spacer):,}")


if __name__ == "__main__":
    main()