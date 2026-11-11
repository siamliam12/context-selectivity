# C2 step 2: finds the reference (REF) and alternative (ALT) letter for each variant from the hg38 genome.
# Reads the genome and positions only. No scores.
# Fix (9 Nov): each letter string is checked on its own; if both fit (indels in repeats), A1 is taken
# provisionally and flagged, to be confirmed against the score file in step 4.
import pandas as pd
import pysam

IN = "decomposition/data/interim/liftover_hg38.csv"
FASTA = "decomposition/data/raw/genome/hg38.fa"
PILOT = "decomposition/data/interim/pilot30.vcf"
OUT = "decomposition/data/interim/ref_alt_hg38.csv"

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
fa = pysam.FastaFile(FASTA)
v = pd.read_csv(IN)

def fits(chrom, pos, allele):
    return fa.fetch(chrom, pos - 1, pos - 1 + len(allele)).upper() == allele

def ref_alt(r):
    a1, a2 = str(r.A1).upper(), str(r.A2).upper()
    is_snp = len(a1) == 1 and len(a2) == 1
    if r.liftover_status not in ("ok", "minus_strand"):
        return None, None, "no_hg38_position"
    if r.liftover_status == "minus_strand":
        if not is_snp:
            return None, None, "minus_strand_indel"
        a1, a2 = COMP[a1], COMP[a2]          # read from the other strand
    chrom, pos = r.chr_hg38, int(r.BP_hg38)
    m1, m2 = fits(chrom, pos, a1), fits(chrom, pos, a2)
    if m1 and not m2:
        return a1, a2, "A1_is_ref"
    if m2 and not m1:
        return a2, a1, "A2_is_ref"
    if m1 and m2:
        return a1, a2, "both_match_A1_assumed"   # indel in a repeat: genome can't decide
    return None, None, "neither_matches"

res = v.apply(ref_alt, axis=1, result_type="expand")
res.columns = ["REF", "ALT", "ref_check"]
v = pd.concat([v, res], axis=1)
v["is_indel"] = (v["A1"].str.len() > 1) | (v["A2"].str.len() > 1)
v["palindromic"] = (~v["is_indel"]) & v.apply(
    lambda r: COMP.get(str(r.A1).upper()) == str(r.A2).upper(), axis=1)

print("Liftover status x REF check:")
print(pd.crosstab(v["liftover_status"], v["ref_check"], margins=True).to_string())

ok = v[v["liftover_status"] == "ok"]
print("\nSNPs and indels (liftover ok only):")
print(pd.crosstab(ok["is_indel"].map({False: "SNP", True: "indel"}), ok["ref_check"]).to_string())

print(f"\nPalindromic SNPs (A/T or C/G): {int(v['palindromic'].sum())}, "
      f"on a flipped strand: {int((v['palindromic'] & (v['liftover_status'] == 'minus_strand')).sum())}")

bad = v[v["ref_check"] == "neither_matches"]
print(f"\nNeither letter matches the genome: {len(bad)}")

snp_a2 = v[(~v["is_indel"]) & (v["ref_check"] == "A2_is_ref")]
print(f"\nSNPs where A2 is REF: {len(snp_a2)}")
print(snp_a2[["variant", "SNP", "chr_hg38", "BP_hg38", "A1", "A2", "REF", "ALT"]].to_string(index=False))

pil = pd.read_csv(PILOT, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3, 4],
                  names=["chrom", "pos", "id", "ref", "alt"])
cmp = pil.merge(v, left_on="id", right_on="SNP", how="left")
agree = (cmp["ref"] == cmp["REF"]) & (cmp["alt"] == cmp["ALT"])
print(f"\nPilot: {len(pil)} variants, {int(agree.sum())} with the same REF and ALT")
print(cmp.loc[~agree, ["id", "pos", "ref", "alt", "A1", "A2", "REF", "ALT", "ref_check"]].to_string(index=False))

v.to_csv(OUT, index=False)
print(f"\nsaved {OUT}: {len(v):,} rows")