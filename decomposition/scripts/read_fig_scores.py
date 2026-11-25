# D2 piece 1: reads the standardized scores of the 89 Fig. 4C tracks for all 5,704 variants
# (sets A, B and C), by score_row, and saves them in long form for Python and DuckDB.
import gcsfs, h5py
import numpy as np
import pandas as pd
from coord_check import check_table

V = "decomposition/data/frozen/variants.csv"
TM = "decomposition/data/frozen/track_map.csv"
import sys
WHICH = sys.argv[1] if len(sys.argv) > 1 else "dm"
SRC = "seqnn-share/sniff/borzoi_veps/" + ("scores.h5" if WHICH == "raw" else "scores_indel_snp_distmatched.h5")
OUT = f"decomposition/data/interim/fig_scores_long{'' if WHICH == 'dm' else '_' + WHICH}.parquet"

v = pd.read_csv(V)
v["score_chr"] = "chr" + v["CHR"].astype(str)
check_table(v, chrom_col="score_chr", pos_col="BP_hg19", ref_col="score_ref", label="variants", build="hg19")

tm = pd.read_csv(TM).sort_values("col")
cols = tm["col"].to_numpy()
first, last = int(cols.min()), int(cols.max())

fs = gcsfs.GCSFileSystem(token="anon")
with fs.open(SRC, "rb") as f, h5py.File(f, "r") as h:
    ds = h["logD2"]
    base, item, ncol, dt = ds.id.get_offset(), ds.dtype.itemsize, ds.shape[1], ds.dtype
    check = ds[int(v["score_row"].iloc[0]), first:last + 1]        # read by h5py, for the self-test

size = (last - first + 1) * item
starts = [int(base + (int(r) * ncol + first) * item) for r in v["score_row"]]
raw = []
for i in range(0, len(starts), 1000):
    s = starts[i:i + 1000]
    raw += fs.cat_ranges([SRC] * len(s), s, [x + size for x in s])
mat = np.vstack([np.frombuffer(b, dtype=dt) for b in raw])
assert np.array_equal(mat[0], check), "byte arithmetic does not match h5py"
mat = mat[:, cols - first].astype(np.float64)
print(f"scores read: {mat.shape[0]:,} variants x {mat.shape[1]} tracks; missing values: {int(np.isnan(mat).sum())}")
print(f"score range: {mat.min():.3f} to {mat.max():.3f}")

long = (pd.DataFrame(mat, columns=cols).assign(variant=v["variant"])
          .melt(id_vars="variant", var_name="col", value_name="score"))
long["col"] = long["col"].astype(int)
long = long.merge(tm[["col", "group"]], on="col")
long.to_parquet(OUT, index=False)
print(f"saved {OUT}: {len(long):,} rows")