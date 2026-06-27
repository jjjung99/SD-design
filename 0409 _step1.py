import RNA
import pandas as pd
import random
from itertools import product
from multiprocessing import Pool, cpu_count

# ==================================================
# Parameters
# ==================================================
SD_LEN = 6
SPACER_LEN = 9
CDS_PREFIX_LEN = 6

BASES = ["A","T","C","G"]

DG_MIN = -9
DG_MAX = -7.5

# ==================================================
# Load CDS (GFP 기반)
# ==================================================
def load_cds(path):
    seqs = []
    for enc in ["utf-8", "cp949", "latin-1"]:
        try:
            with open(path, encoding=enc) as f:
                for line in f:
                    line = line.strip().upper()
                    if not line or line.startswith(">"):
                        continue
                    clean = "".join([b for b in line if b in "ATCG"])
                    if len(clean) >= 20:
                        seqs.append(clean)
            print(f"[INFO] CDS loaded with {enc}")
            return seqs
        except:
            continue
    raise ValueError("CDS load failed")

CODING_SEQS = load_cds("C:/Users/Minnie/Desktop/OSD_design/data/o-SD test gfp.dna")

# ==================================================
# Utils
# ==================================================
def dna_to_rna(seq):
    return seq.replace("T","U")

def revcomp_dna(seq):
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(comp.get(b, "") for b in reversed(seq))

def duplex_dG(a,b):
    try:
        return RNA.duplexfold(a,b).energy
    except:
        return 1000

# ==================================================
# Generate
# ==================================================
def generate_sd():
    for tup in product(BASES, repeat=SD_LEN):
        yield "".join(tup)

def generate_spacer():
    for tup in product(BASES, repeat=SPACER_LEN):
        yield "".join(tup)

# ==================================================
# Binding filter 
# ==================================================
def check_binding(args):
    sd, spacer = args

    # CDS prefix 추가
    cds = random.choice(CODING_SEQS)
    cds_prefix = cds[:CDS_PREFIX_LEN]

    # target (15bp window)
    target = sd + spacer + cds_prefix
    target_rna = dna_to_rna(target)

    o_asd = revcomp_dna(target)
    o_asd_rna = dna_to_rna(o_asd)

    dg = duplex_dG(o_asd_rna, target_rna)

    if DG_MIN <= dg <= DG_MAX:
        return {
            "SD": sd,
            "Spacer": spacer,
            "dg_bind": dg
        }

    return None

# ==================================================
# Main
# ==================================================
def main():

    sds = list(generate_sd())
    spacers = list(generate_spacer())

    combos = [(s, p) for s in sds for p in spacers]

    print(f"[INFO] Total combinations: {len(combos)}")

    results = []

    with Pool(cpu_count()-1) as pool:
        for r in pool.imap_unordered(check_binding, combos, chunksize=1000):
            if r:
                results.append(r)

    df = pd.DataFrame(results)

    print(f"[RESULT] Binding candidates: {len(df)}")

    df.to_csv("binding_range_candidates_L9.csv", index=False)

    print(df.head(10))


if __name__ == "__main__":
    main()