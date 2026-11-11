# C2 step 8a: builds the two frozen C2 tables.
#   variants.csv        one row per unique variant: positions in both genomes, the authors' letters,
#                       the score-file row, and the groups (mutation class, CpG, GC) computed in hg19.
#   variant_traits.csv  one row per variant and trait: the C1 table, frozen as it is.
# Reads the genomes and our own tables only. No scores.
import os
import pandas as pd
import pysam
from coord_check import check_table

CTX = "decomposition/data/interim/variant_context_hg38.csv"
ROWS = "decomposition/data/interim/score_rows.csv"
POP = "decomposition/data/interim/variant_population.csv"
OUT_V = "decomposition/data/frozen/variants.csv"
OUT_T = "decomposition/data/frozen/variant_traits.csv"
GC_FLANK = 50
TRANSITIONS = [{"A", "G"}, {"C", "T"}]

ctx = pd.read_csv(CTX)
rows = pd.read_csv(ROWS)
pop = pd.read_csv(POP)
v = ctx.merge(rows, on="variant", how="left", validate="one_to_one")
assert len(v) == 5704 and v["score_row"].notna().all(), "every variant must have a score-file row"
assert ((v["score_chr"] == "chr" + v["CHR"].astype(str)) & (v["score_pos"] == v["BP_hg19"])).all(), \
    "score-file position must equal the hg19 position"

# --- Groups in hg19, the genome the scores were computed on ---
fa19 = pysam.FastaFile("decomposition/data/raw/genome/hg19.fa")
def context_hg19(r):
    ref, alt = r.s_ref, r.s_alt
    if len(ref) != 1 or len(alt) != 1:
        return None, None, None
    chrom, pos = r.score_chr, int(r.score_pos)
    mclass = "transition" if {ref, alt} in TRANSITIONS else "transversion"
    left = fa19.fetch(chrom, pos - 2, pos - 1).upper()
    right = fa19.fetch(chrom, pos, pos + 1).upper()
    cpg = (ref == "C" and right == "G") or (ref == "G" and left == "C")
    win = fa19.fetch(chrom, pos - 1 - GC_FLANK, pos + GC_FLANK).upper()
    acgt = sum(win.count(b) for b in "ACGT")
    gc = (win.count("G") + win.count("C")) / acgt if acgt else None
    return mclass, cpg, gc

v = v.rename(columns={"mut_class": "mut_class_hg38", "cpg": "cpg_hg38", "gc_101bp": "gc_101bp_hg38"})
res = v.apply(context_hg19, axis=1, result_type="expand")
res.columns = ["mut_class", "cpg", "gc_101bp"]
v = pd.concat([v, res], axis=1)

both = v[v["mut_class"].notna() & v["mut_class_hg38"].notna()]
print(f"SNPs with groups in both genomes: {len(both):,}")
print(f"mutation class same in hg19 and hg38: {(both['mut_class'] == both['mut_class_hg38']).mean():.4f}")
print(f"CpG same in hg19 and hg38: {(both['cpg'].astype(bool) == both['cpg_hg38'].astype(bool)).mean():.4f}")
print(f"GC correlation, hg19 vs hg38: {both['gc_101bp'].corr(both['gc_101bp_hg38']):.4f}")

# --- Permanent position checks, both genomes ---
check_table(v, chrom_col="score_chr", pos_col="score_pos", ref_col="s_ref",
            label="frozen variants (score file, hg19)", build="hg19")
check_table(v, label="frozen variants (our hg38 side)", build="hg38")

# --- Final columns ---
out = pd.DataFrame({
    "variant": v["variant"], "SNP": v["SNP"], "CHR": v["CHR"], "BP_hg19": v["BP_hg19"],
    "A1": v["A1"], "A2": v["A2"], "is_indel": v["is_indel"], "palindromic": v["palindromic"],
    "score_row": v["score_row"].astype("int64"),
    "score_ref": v["s_ref"], "score_alt": v["s_alt"],
    "score_ref_is": v["letters"].map({"same": "A1", "swapped": "A2"}),
    "mut_class": v["mut_class"], "cpg": v["cpg"], "gc_101bp": v["gc_101bp"],
    "chr_hg38": v["chr_hg38"], "BP_hg38": v["BP_hg38"].astype("Int64"), "strand_hg38": v["strand"],
    "liftover_status": v["liftover_status"], "REF_hg38": v["REF"], "ALT_hg38": v["ALT"],
    "ref_check_hg38": v["ref_check"],
    "mut_class_hg38": v["mut_class_hg38"], "cpg_hg38": v["cpg_hg38"], "gc_101bp_hg38": v["gc_101bp_hg38"],
})
assert set(pop["variant"]) == set(out["variant"]), "the two tables must hold the same variants"
print(f"\nvariants: {len(out):,} ({int(out['is_indel'].sum())} indels); variant-trait rows: {len(pop):,}")

def write_once(df, path):
    text = df.to_csv(index=False)
    if os.path.exists(path):
        assert open(path).read() == text, f"{path} exists and differs; frozen files are never overwritten"
        print(f"{path} already exists and matches")
    else:
        with open(path, "w") as f:
            f.write(text)
        print(f"saved {path}: {len(df):,} rows, {df.shape[1]} columns")

write_once(out, OUT_V)
write_once(pop, OUT_T)