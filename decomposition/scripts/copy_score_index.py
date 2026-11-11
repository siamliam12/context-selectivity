# C2 step 4a: copies the score file's variant list (rsID, chromosome, position) to a local file,
# with each variant's row number in the score matrix. Never reads the score matrix or the letters.
import time
import gcsfs, h5py
import numpy as np
import pandas as pd

fs = gcsfs.GCSFileSystem(token="anon")
SRC = "seqnn-share/sniff/borzoi_veps/scores_indel_snp_distmatched.h5"
OUT = "decomposition/data/interim/score_index.parquet"

t0 = time.time()
with fs.open(SRC, "rb", block_size=64 * 2**20) as f, h5py.File(f, "r") as h:
    for name in ["snp", "chr", "pos", "ref_allele", "alt_allele", "logD2"]:
        d = h[name]
        print(f"{name}: chunks {d.chunks}, compression {d.compression}", flush=True)
    idx = pd.DataFrame({
        "snp": pd.Series(h["snp"][:]).str.decode("ascii"),
        "chr": pd.Series(h["chr"][:]).str.decode("ascii").astype("category"),
        "pos": h["pos"][:].astype(np.int64),
    })
    print(f"read 3 lists in {time.time() - t0:.0f} s", flush=True)

idx.insert(0, "row", np.arange(len(idx), dtype=np.int64))
print(idx.head().to_string(index=False))
print(f"\n{len(idx):,} rows")
print("Rows per chromosome:")
print(idx["chr"].value_counts().sort_index().to_string())
print(f"Duplicate rsIDs: {int(idx['snp'].duplicated().sum()):,}")
print(f"Duplicate chr+pos: {int(idx.duplicated(['chr', 'pos']).sum()):,}")
sorted_ok = idx.groupby("chr", observed=True)["pos"].apply(lambda s: s.is_monotonic_increasing).all()
print(f"Positions sorted within each chromosome: {sorted_ok}")

idx.to_parquet(OUT, index=False)
print(f"\nsaved {OUT}")