# C2 step 4f: checks the score file's REF letters against the hg19 genome, for our variants
# and for 10,000 random rows of the whole file. Reads letters only, never scores.
import gcsfs, h5py
import numpy as np
import pandas as pd
from coord_check import check_table, ref_fits

SR = "decomposition/data/interim/score_rows.csv"
IDX = "decomposition/data/interim/score_index.parquet"
SRC = "seqnn-share/sniff/borzoi_veps/scores_indel_snp_distmatched.h5"
SEED = 20261110

# 1. Our 5,704 variants: the file's REF, at the file's position, in hg19
sr = pd.read_csv(SR)
check_table(sr, chrom_col="score_chr", pos_col="score_pos", ref_col="s_ref",
            label="our variants, file REF", build="hg19")

# 2. 10,000 random rows from the whole file
idx = pd.read_parquet(IDX)
rng = np.random.default_rng(SEED)
pick = np.sort(rng.choice(len(idx), size=10_000, replace=False))
samp = idx.iloc[pick].copy()
samp["chr"] = samp["chr"].astype(str)

fs = gcsfs.GCSFileSystem(token="anon")
def read_strings(h, name, rows):
    ds = h[name]
    base, size = ds.id.get_offset(), ds.dtype.itemsize
    starts = [int(base + int(r) * size) for r in rows]
    out = []
    for i in range(0, len(starts), 1000):
        s = starts[i:i + 1000]
        got = fs.cat_ranges([SRC] * len(s), s, [x + size for x in s])
        out += [b.rstrip(b"\x00").decode("ascii") for b in got]
    return out

with fs.open(SRC, "rb") as f, h5py.File(f, "r") as h:
    samp["s_ref"] = read_strings(h, "ref_allele", samp["row"])
    samp["s_alt"] = read_strings(h, "alt_allele", samp["row"])

samp["is_indel"] = (samp["s_ref"].str.len() > 1) | (samp["s_alt"].str.len() > 1)
for build in ["hg19", "hg38"]:
    samp[f"fits_{build}"] = [ref_fits(c, int(p), r, build)
                             for c, p, r in zip(samp["chr"], samp["pos"], samp["s_ref"])]
print(f"\nrandom sample of {len(samp):,} rows (seed {SEED}), share whose REF fits each genome:")
print(samp.groupby(samp["is_indel"].map({False: "SNP", True: "indel"}))[["fits_hg19", "fits_hg38"]]
      .mean().to_string())
print(f"\nrows whose REF does not fit hg19: {int((~samp['fits_hg19']).sum())}")
print(samp.loc[~samp["fits_hg19"]].head(10).to_string(index=False))