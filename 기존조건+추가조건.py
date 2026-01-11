import os
import csv
from itertools import product
from multiprocessing import Pool, cpu_count

import RNA
import matplotlib.pyplot as plt


# ===============================
# Fixed (불변) 서열
# ===============================
UTR_DNA = "ATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCC"  # 50bp fixed
AU_RICH_DNA = "TTAATTAA"  # 8bp fixed
AUG_DNA = "ATG"

EMPTY_ORF_PREFIX = (
    "CTGCTGGGTGAGCTTTCTCCGTAAACTTAAAGGAAAAGATTCCGTTGAAAGATT"
    "CAAAGCTATCGTTCAGCGTATACAAGAGACTTCCTCCTGAGACTCGTGTTCCC"
    "GTACCGAACTCT"
)

# ===============================
# WT 모티프
# ===============================
WT_SD_DNA = "AGGAGG"
WT_ASD_CORE_DNA = "CTCCTT"

# (추가조건) WT-ASD 12bp (RNA)
WT_ASD_12_RNA = "AUCACCUCCUUA"


# ===============================
# 설계 변수 규칙
# ===============================
SD_LEN = 6
SPACER_LEN = 6

# SD GC 2~4개 허용
GC_MIN = 2
GC_MAX = 4

# Spacer AT-only
SPACER_BASES = ["A", "T"]

# ===============================
# 조건 (필수)
# ===============================
# Translation-competent 강결합 (O-SD : O-ASDcore6)
DG_STRONG_MIN = -10.0
DG_STRONG_MAX = -8.5

# 직교성: ΔG > 0 기준
ORTHO_POSITIVE_CUTOFF = 0.0

# (추가조건) Upstream 15bp 포함 cross-binding 검증
UPSTREAM_LEN = 15

# 구조 확인용 컨텍스트 길이
FOLD_ORF_PREFIX_N = 60


# ===============================
# Output
# ===============================
OUT_DIR = "orthogonal_sd_results_ALL_only"
os.makedirs(OUT_DIR, exist_ok=True)

ALL_CSV = os.path.join(OUT_DIR, "all_combinations.csv")
FINAL_STRICT_CSV = os.path.join(OUT_DIR, "final_strict.csv")
BAR_VIOL = os.path.join(OUT_DIR, "bar_orthogonality_violations.png")

HIST_STRONG = os.path.join(OUT_DIR, "hist_strong_dG_O_SD_O_ASD_core6.png")
HIST_WT_ASDCORE = os.path.join(OUT_DIR, "hist_cross_WT_ASDcore_CTCCTT__O_SD.png")
HIST_WT_SD = os.path.join(OUT_DIR, "hist_cross_WT_SD_AGGAGG__O_ASDcore6.png")
HIST_WT_ASD12_SDSP12 = os.path.join(OUT_DIR, "hist_cross_WT_ASD12__O_SDSp12.png")
HIST_WT_ASD12_UP15 = os.path.join(OUT_DIR, "hist_cross_WT_ASD12__UP15_O_SDSp12.png")
HIST_MFE = os.path.join(OUT_DIR, "hist_hairpin_MFE_context.png")


# ===============================
# Utilities
# ===============================
def dna_to_rna(seq: str) -> str:
    return seq.replace("T", "U")

def revcomp_dna(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))

def duplex_dG(rna_a: str, rna_b: str) -> float:
    d = RNA.duplexfold(rna_a, rna_b)
    return float(d.energy)

def fold_mfe(rna_seq: str) -> float:
    _, mfe = RNA.fold(rna_seq)
    return float(mfe)

def upstream_15_dna() -> str:
    # SD 바로 앞 15bp: (UTR + AU_RICH) 의 마지막 15bp
    upstream = (UTR_DNA + AU_RICH_DNA)
    return upstream[-UPSTREAM_LEN:]

def build_context_rna(sd_dna: str, spacer_dna: str) -> str:
    ctx_dna = (
        UTR_DNA
        + AU_RICH_DNA
        + sd_dna
        + spacer_dna
        + AUG_DNA
        + EMPTY_ORF_PREFIX[:FOLD_ORF_PREFIX_N]
    )
    return dna_to_rna(ctx_dna)

def save_hist(vals, path, title, xlabel):
    plt.figure()
    plt.hist(vals, bins=60)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


# ===============================
# Precompute WT RNAs
# ===============================
WT_SD_RNA = dna_to_rna(WT_SD_DNA)
WT_ASD_CORE_RNA = dna_to_rna(WT_ASD_CORE_DNA)

UP15_DNA = upstream_15_dna()
UP15_RNA = dna_to_rna(UP15_DNA)


# ===============================
# Candidate generation
# ===============================
def generate_sd_list():
    bases = ["A", "T", "C", "G"]
    sds = []
    for tup in product(bases, repeat=SD_LEN):
        sd = "".join(tup)
        gc = sum(1 for b in sd if b in ("G", "C"))
        if GC_MIN <= gc <= GC_MAX:
            sds.append(sd)
    return sds

def generate_spacer_list():
    return ["".join(tup) for tup in product(SPACER_BASES, repeat=SPACER_LEN)]


# ===============================
# Score one combo (ALL 조건 계산 + strict pass 여부)
# ===============================
def score_one(args):
    sd_dna, spacer_dna = args

    osd_rna = dna_to_rna(sd_dna)
    spacer_rna = dna_to_rna(spacer_dna)

    # SD+Spacer (12bp)
    osdsp12_rna = osd_rna + spacer_rna

    # O-ASD core (6bp) = revcomp(O-SD)
    o_asd_core_dna = revcomp_dna(sd_dna)
    oasd_core_rna = dna_to_rna(o_asd_core_dna)

    # (필수 A) Translation-competent strong binding: O-SD vs O-ASDcore6
    dg_strong = duplex_dG(osd_rna, oasd_core_rna)
    strong_ok = (DG_STRONG_MIN <= dg_strong <= DG_STRONG_MAX)

    # (필수 B1) WT-ASD(core CTCCTT) : O-SD  -> ΔG > 0
    dg_WT_ASDcore__O_SD = duplex_dG(WT_ASD_CORE_RNA, osd_rna)
    ortho_core1_ok = (dg_WT_ASDcore__O_SD > ORTHO_POSITIVE_CUTOFF)

    # (필수 B2) WT-SD(AGGAGG) : O-ASD(core) -> ΔG > 0
    dg_WT_SD__O_ASDcore = duplex_dG(WT_SD_RNA, oasd_core_rna)
    ortho_core2_ok = (dg_WT_SD__O_ASDcore > ORTHO_POSITIVE_CUTOFF)

    # (추가조건 C1) WT-ASD12 : O-(SD+Spacer)12 (계산은 저장)
    dg_WT_ASD12__O_SDSp12 = duplex_dG(WT_ASD_12_RNA, osdsp12_rna)

    # (추가조건 C2, 필수) WT-ASD12 : UP15 + O-(SD+Spacer)12 -> ΔG > 0
    up15_osdsp27_rna = UP15_RNA + osdsp12_rna
    dg_WT_ASD12__UP15_O_SDSp27 = duplex_dG(WT_ASD_12_RNA, up15_osdsp27_rna)
    ortho_up15_ok = (dg_WT_ASD12__UP15_O_SDSp27 > ORTHO_POSITIVE_CUTOFF)

    # (옵션) 구조 (랭킹에서 쓸 예정)
    mfe_ctx = fold_mfe(build_context_rna(sd_dna, spacer_dna))

    # strict 통과 여부 = ALL 필수 조건 통과
    strict_ok = (strong_ok and ortho_core1_ok and ortho_core2_ok and ortho_up15_ok)

    # 위반 카운트(분포용)
    # 여기서는 "필수 직교성/추가조건" 위반만 카운트 (strong_ok는 기능성이라 따로)
    v_core1 = int(not ortho_core1_ok)
    v_core2 = int(not ortho_core2_ok)
    v_up15 = int(not ortho_up15_ok)
    viol_count = v_core1 + v_core2 + v_up15

    return {
        "SD_DNA": sd_dna,
        "Spacer_DNA": spacer_dna,
        "UP15_DNA": UP15_DNA,
        "O_ASD_core_DNA": o_asd_core_dna,

        "dg_O_SD__O_ASD_core6": dg_strong,
        "strong_ok_core6": int(strong_ok),

        "dg_WT_ASDcore__O_SD": dg_WT_ASDcore__O_SD,
        "ortho_core1_ok": int(ortho_core1_ok),

        "dg_WT_SD__O_ASDcore": dg_WT_SD__O_ASDcore,
        "ortho_core2_ok": int(ortho_core2_ok),

        "dg_WT_ASD12__O_SDSp12": dg_WT_ASD12__O_SDSp12,

        "dg_WT_ASD12__UP15_O_SDSp27": dg_WT_ASD12__UP15_O_SDSp27,
        "ortho_up15_ok": int(ortho_up15_ok),

        "hairpin_MFE_context": mfe_ctx,

        "violations_count": viol_count,
        "final_strict_ok": int(strict_ok),
    }


# ===============================
# Main
# ===============================
def main():
    sds = generate_sd_list()
    spacers = generate_spacer_list()
    total = len(sds) * len(spacers)

    print(f"[INFO] SD candidates: {len(sds)} (GC {GC_MIN}~{GC_MAX})")
    print(f"[INFO] Spacer candidates: {len(spacers)} (AT-only)")
    print(f"[INFO] Total combos: {total}")
    print(f"[INFO] WT-ASD12 (RNA): {WT_ASD_12_RNA}")
    print(f"[INFO] WT-ASD core (RNA): {WT_ASD_CORE_RNA}")
    print(f"[INFO] WT-SD (RNA): {WT_SD_RNA}")
    print(f"[INFO] UP15 (DNA): {UP15_DNA}")

    combos = [(sd, sp) for sd in sds for sp in spacers]
    nproc = max(1, cpu_count() - 1)
    print(f"[INFO] Using processes: {nproc}")

    # Columns
    fieldnames = list(score_one((sds[0], spacers[0])).keys())

    # Collect for strict and plots
    strict_rows = []
    strong_vals = []
    wt_asdcore_vals = []
    wt_sd_vals = []
    wt_asd12_sdsp12_vals = []
    wt_asd12_up15_vals = []
    mfe_vals = []
    viol_hist = [0, 0, 0, 0]  # 0~3

    with open(ALL_CSV, "w", newline="") as f_all:
        w_all = csv.DictWriter(f_all, fieldnames=fieldnames)
        w_all.writeheader()

        with Pool(processes=nproc) as pool:
            for row in pool.imap_unordered(score_one, combos, chunksize=800):
                w_all.writerow(row)

                # collect distributions
                strong_vals.append(row["dg_O_SD__O_ASD_core6"])
                wt_asdcore_vals.append(row["dg_WT_ASDcore__O_SD"])
                wt_sd_vals.append(row["dg_WT_SD__O_ASDcore"])
                wt_asd12_sdsp12_vals.append(row["dg_WT_ASD12__O_SDSp12"])
                wt_asd12_up15_vals.append(row["dg_WT_ASD12__UP15_O_SDSp27"])
                mfe_vals.append(row["hairpin_MFE_context"])

                vc = row["violations_count"]
                if vc < 0: vc = 0
                if vc > 3: vc = 3
                viol_hist[vc] += 1

                # strict rows
                if row["final_strict_ok"] == 1:
                    strict_rows.append(row)

    # Save strict file
    with open(FINAL_STRICT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(strict_rows)

    print("\n[SAVED CSV]")
    print(f"- {ALL_CSV}")
    print(f"- {FINAL_STRICT_CSV} (n={len(strict_rows)})")

    # Save plots
    save_hist(strong_vals, HIST_STRONG, "Strong ΔG: O-SD vs O-ASD(core6)", "ΔG (kcal/mol)")
    save_hist(wt_asdcore_vals, HIST_WT_ASDCORE, "Cross ΔG: WT-ASD(core CTCCTT) vs O-SD", "ΔG (kcal/mol)")
    save_hist(wt_sd_vals, HIST_WT_SD, "Cross ΔG: WT-SD(AGGAGG) vs O-ASD(core)", "ΔG (kcal/mol)")
    save_hist(wt_asd12_sdsp12_vals, HIST_WT_ASD12_SDSP12, "Cross ΔG: WT-ASD12 vs O-(SD+Spacer)12", "ΔG (kcal/mol)")
    save_hist(wt_asd12_up15_vals, HIST_WT_ASD12_UP15, "Cross ΔG: WT-ASD12 vs UP15 + O-(SD+Spacer)12", "ΔG (kcal/mol)")
    save_hist(mfe_vals, HIST_MFE, "Hairpin MFE (context folding)", "MFE (kcal/mol)")

    plt.figure()
    plt.bar([0, 1, 2, 3], viol_hist)
    plt.title("Orthogonality violations count (core1/core2/up15) 0~3")
    plt.xlabel("violations_count")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(BAR_VIOL, dpi=200)
    plt.close()

    print("\n[SAVED PLOTS]")
    print(f"- {HIST_STRONG}")
    print(f"- {HIST_WT_ASDCORE}")
    print(f"- {HIST_WT_SD}")
    print(f"- {HIST_WT_ASD12_SDSP12}")
    print(f"- {HIST_WT_ASD12_UP15}")
    print(f"- {HIST_MFE}")
    print(f"- {BAR_VIOL}")


if __name__ == "__main__":
    main()
