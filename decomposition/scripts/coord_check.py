"""Permanent position check (C2 step 6; hg19 added in step 4f).

Every script that uses variant positions imports this and runs it on every variant.
It confirms that the REF letters sit at the stated position in the stated genome.
It is never switched off.
"""
import pysam

FASTA = {
    "hg38": "decomposition/data/raw/genome/hg38.fa",
    "hg19": "decomposition/data/raw/genome/hg19.fa",
}
_open = {}

def _genome(build):
    if build not in _open:
        _open[build] = pysam.FastaFile(FASTA[build])
    return _open[build]

def ref_fits(chrom, pos, ref, build="hg38"):
    """True if the genome at chrom:pos (counting from 1) reads ref."""
    return _genome(build).fetch(chrom, pos - 1, pos - 1 + len(ref)).upper() == str(ref).upper()

def check_table(df, chrom_col="chr_hg38", pos_col="BP_hg38", ref_col="REF", label="table", build="hg38"):
    """Checks every row that has a REF. Stops with a clear error if any row fails."""
    rows = df[df[ref_col].notna()]
    bad = [i for i, r in rows.iterrows()
           if not ref_fits(r[chrom_col], int(r[pos_col]), str(r[ref_col]), build)]
    if bad:
        raise AssertionError(f"{label}: {len(bad)} of {len(rows)} variants fail the {build} position check, "
                             f"e.g. rows {bad[:5]}")
    print(f"{label}: {build} position check passed for {len(rows):,} variants")
    return len(rows)