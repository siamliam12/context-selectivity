# D5 step 1: draws 50,000 random variants from the score file (not the figure's), checks their hg19
# letters, converts them to hg38 (C2 method), and reads their 89 standardized scores.
import numpy as np, pandas as pd, gcsfs, h5py
from pyliftover import LiftOver
from coord_check import check_table

SEED, N = 20261116 + 5, 50_000                 # D5's seed, derived from the master seed
IDX = "decomposition/data/interim/score_index.parquet"
CHAIN = "decomposition/data/raw/genome/hg19ToHg38.over.chain.gz"
SRC = "seqnn-share/sniff/borzoi_veps/scores_indel_snp_distmatched.h5"
OUT_V = "decomposition/data/interim/d5_sample_variants.csv"
OUT_S = "decomposition/data/interim/d5_scores_long.parquet"
tm = pd.read_csv("decomposition/data/frozen/track_map.csv").sort_values("col")

# 1. Draw the sample, leaving out the figure's variants
idx = pd.read_parquet(IDX)
fig_rows = set(pd.read_csv("decomposition/data/frozen/variants.csv")["score_row"])
pool = idx[~idx["row"].isin(fig_rows)]
rng = np.random.default_rng(SEED)
s = pool.iloc[np.sort(rng.choice(len(pool), N, replace=False))].copy().reset_index(drop=True)
s["chr"] = s["chr"].astype(str)

# 2. Letters for the sampled rows, then the hg19 position check
fs = gcsfs.GCSFileSystem(token="anon")
def fetch(starts, size):
    out = []
    for i in range(0, len(starts), 1000):
        b = starts[i:i + 1000]
        out += fs.cat_ranges([SRC] * len(b), b, [x + size for x in b])
    return out
with fs.open(SRC, "rb") as f, h5py.File(f, "r") as h:
    info = {k: (h[k].id.get_offset(), h[k].dtype.itemsize) for k in ["ref_allele", "alt_allele"]}
    ds = h["logD2"]
    base, item, ncol, dt = ds.id.get_offset(), ds.dtype.itemsize, ds.shape[1], ds.dtype
    first, last = int(tm["col"].min()), int(tm["col"].max())
    check = ds[int(s["row"].iloc[0]), first:last + 1]
for k, name in [("ref_allele", "ref"), ("alt_allele", "alt")]:
    off, it = info[k]
    s[name] = [x.rstrip(b"\x00").decode() for x in fetch([int(off + int(r) * it) for r in s["row"]], it)]
check_table(s, chrom_col="chr", pos_col="pos", ref_col="ref", label="D5 sample", build="hg19")
s["variant"] = s["chr"].str.replace("chr", "") + ":" + s["pos"].astype(str) + ":" + s["ref"] + ":" + s["alt"]
s["is_indel"] = (s["ref"].str.len() > 1) | (s["alt"].str.len() > 1)

# 3. hg38 positions (failures counted and dropped later)
lo = LiftOver(CHAIN)
def lift(c, p):
    hits = lo.convert_coordinate(c, int(p) - 1)
    if not hits or len(hits) > 1:
        return None, None, "unmapped" if not hits else "multiple"
    c2, p0, strand, _ = hits[0]
    if c2 != c:
        return None, None, "other_chrom"
    return c2, p0 + 1, "ok" if strand == "+" else "minus_strand"
lifted = [lift(c, p) for c, p in zip(s["chr"], s["pos"])]
s["chr_hg38"], s["BP_hg38"], s["liftover_status"] = zip(*lifted)
print("liftover:", s["liftover_status"].value_counts().to_dict())

# 4. The 89 scores (same byte method and self-test as before)
cols = tm["col"].to_numpy()
size = (last - first + 1) * item
mat = np.vstack([np.frombuffer(b, dtype=dt) for b in
                 fetch([int(base + (int(r) * ncol + first) * item) for r in s["row"]], size)])
assert np.array_equal(mat[0], check), "byte arithmetic does not match h5py"
mat = mat[:, cols - first].astype(np.float64)
long = (pd.DataFrame(mat, columns=cols).assign(variant=s["variant"].values)
          .melt(id_vars="variant", var_name="col", value_name="score"))
long["col"] = long["col"].astype(int)
long.merge(tm[["col", "group"]], on="col").to_parquet(OUT_S, index=False)
s.drop(columns=[]).to_csv(OUT_V, index=False)
print(f"sample: {len(s):,} variants ({int(s['is_indel'].sum()):,} indels); saved {OUT_V} and {OUT_S}")