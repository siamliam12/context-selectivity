# C2 step 1: converts every variant's position from hg19 to hg38.
# Reads positions only. No scores.
import os
import pandas as pd
from pyliftover import LiftOver

POP = "decomposition/data/interim/variant_population.csv"
CHAIN = "decomposition/data/raw/genome/hg19ToHg38.over.chain.gz"
PILOT = "decomposition/data/interim/pilot30.vcf"
OUT = "decomposition/data/interim/liftover_hg38.csv"

assert os.path.exists(CHAIN), f"chain file not found: {CHAIN}"
lo = LiftOver(CHAIN)

pop = pd.read_csv(POP)
v = pop[["variant", "SNP", "CHR", "BP_hg19", "A1", "A2"]].drop_duplicates("variant").reset_index(drop=True)
print(f"{len(pop):,} variant-trait rows -> {len(v):,} unique variants")

rows = []
for r in v.itertuples(index=False):
    hits = lo.convert_coordinate(f"chr{r.CHR}", int(r.BP_hg19) - 1)  # pyliftover counts from 0
    if not hits:
        status, chrom, pos, strand = "unmapped", None, None, None
    elif len(hits) > 1:
        status, chrom, pos, strand = "multiple", None, None, None
    else:
        chrom, pos0, strand, _ = hits[0]
        pos = pos0 + 1  # back to counting from 1
        if chrom != f"chr{r.CHR}":
            status = "other_chrom"
        elif strand != "+":
            status = "minus_strand"
        else:
            status = "ok"
    rows.append((r.variant, chrom, pos, strand, status))

lift = pd.DataFrame(rows, columns=["variant", "chr_hg38", "BP_hg38", "strand", "liftover_status"])
v = v.merge(lift, on="variant", validate="one_to_one")

print("\nLiftover result:")
print(v["liftover_status"].value_counts().to_string())
print(f"Share not ok: {(v['liftover_status'] != 'ok').mean():.4%}")
print("\nRows that are not ok (first 50):")
print(v[v["liftover_status"] != "ok"].head(50).to_string(index=False))

# --- Compare with the Part A pilot, which was already lifted to hg38 ---
pil = pd.read_csv(PILOT, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["chrom", "pos", "id", "ref", "alt"])
cmp = pil.merge(v, left_on="id", right_on="SNP", how="left")
found = int(cmp["BP_hg38"].notna().sum())
same = int((cmp["pos"] == cmp["BP_hg38"]).sum())
print(f"\nPilot: {len(pil)} variants, {found} found by rsID, {same} with the same hg38 position")

v.to_csv(OUT, index=False)
print(f"\nsaved {OUT}: {len(v):,} rows")