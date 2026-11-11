# C2 step 5: records mutation class, CpG status and local GC content for each SNP.
# Fixed now, before any result, as groups to split results by later.
# Reads the genome only. No scores.
import pandas as pd
import pysam

IN = "decomposition/data/interim/ref_alt_hg38.csv"
FASTA = "decomposition/data/raw/genome/hg38.fa"
OUT = "decomposition/data/interim/variant_context_hg38.csv"
GC_FLANK = 50   # bases on each side of the variant (window = 101 bases)

fa = pysam.FastaFile(FASTA)
v = pd.read_csv(IN)
TRANSITIONS = [{"A", "G"}, {"C", "T"}]

def context(r):
    if r.is_indel or pd.isna(r.REF):
        return None, None, None, None
    chrom, pos = r.chr_hg38, int(r.BP_hg38)
    ref, alt = r.REF, r.ALT
    mclass = "transition" if {ref, alt} in TRANSITIONS else "transversion"
    left = fa.fetch(chrom, pos - 2, pos - 1).upper()      # the base before
    right = fa.fetch(chrom, pos, pos + 1).upper()          # the base after
    cpg = (ref == "C" and right == "G") or (ref == "G" and left == "C")
    win = fa.fetch(chrom, max(0, pos - 1 - GC_FLANK), pos + GC_FLANK).upper()
    acgt = sum(win.count(b) for b in "ACGT")
    gc = (win.count("G") + win.count("C")) / acgt if acgt else None
    return mclass, cpg, gc, len(win) - acgt

res = v.apply(context, axis=1, result_type="expand")
res.columns = ["mut_class", "cpg", "gc_101bp", "n_N_in_window"]
v = pd.concat([v, res], axis=1)

snp = v[(~v["is_indel"]) & v["REF"].notna()].copy()
snp["cpg"] = snp["cpg"].astype(bool)
print(f"SNPs with REF and ALT: {len(snp):,}")
print("\nMutation class:")
print(snp["mut_class"].value_counts().to_string())
ts = (snp["mut_class"] == "transition").sum()
tv = (snp["mut_class"] == "transversion").sum()
print(f"Transition / transversion ratio: {ts / tv:.2f}")
print(f"\nCpG sites: {int(snp['cpg'].sum()):,} ({snp['cpg'].mean():.1%})")
print(f"Transitions at CpG sites: {int((snp['cpg'] & (snp['mut_class'] == 'transition')).sum()):,}")
print("\nGC content in the 101-base window:")
print(snp["gc_101bp"].describe().to_string())
print(f"\nWindows containing unknown bases (N): {int((snp['n_N_in_window'] > 0).sum())}")

v.to_csv(OUT, index=False)
print(f"\nsaved {OUT}: {len(v):,} rows")