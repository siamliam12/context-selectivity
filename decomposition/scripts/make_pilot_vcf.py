import random
import pandas as pd
from pyliftover import LiftOver
from pyfaidx import Fasta

S3 = "data/raw/sniff/media-1.xlsx"
CHAIN = "data/raw/genome/hg19ToHg38.over.chain.gz"
GENOME = "data/raw/genome/hg38.fa"
SEED = 20261028

# One row per unique variant. Traits are not needed for timing.
d = pd.read_excel(S3, "S3").drop_duplicates("SNP")
d["is_indel"] = (d["A1"].str.len() > 1) | (d["A2"].str.len() > 1)

lo = LiftOver(CHAIN)
genome = Fasta(GENOME)

def to_vcf_row(r):
    chrom = f"chr{r.CHR}"
    hit = lo.convert_coordinate(chrom, int(r.BP) - 1)   # liftover uses 0-based positions
    if not hit:
        return None                                      # failed to convert
    chrom38, pos0 = hit[0][0], hit[0][1]
    pos = pos0 + 1                                       # back to 1-based for the VCF
    for ref, alt in [(r.A1, r.A2), (r.A2, r.A1)]:
        if str(genome[chrom38][pos0:pos0 + len(ref)]).upper() == ref.upper():
            return f"{chrom38}\t{pos}\t{r.SNP}\t{ref}\t{alt}\t.\t.\t."
    return None                                          # neither allele matches the genome

def pick(pool, n, rng):
    rows = []
    for i in rng.sample(range(len(pool)), len(pool)):    # random order
        row = to_vcf_row(pool.iloc[i])
        if row:
            rows.append(row)
        if len(rows) == n:
            return rows

rng = random.Random(SEED)
snps = d[~d["is_indel"]].reset_index(drop=True)
indels = d[d["is_indel"]].reset_index(drop=True)

pilot30 = pick(snps, 24, rng) + pick(indels, 6, rng)    # 30 variants, 20% indels
pilot10 = pilot30[:8] + pilot30[24:26]                   # 8 SNPs + 2 indels from the same set

header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
for name, rows in [("pilot10", pilot10), ("pilot30", pilot30)]:
    with open(f"data/interim/{name}.vcf", "w") as f:
        f.write(header + "\n".join(rows) + "\n")
    print(name, len(rows), "variants written")