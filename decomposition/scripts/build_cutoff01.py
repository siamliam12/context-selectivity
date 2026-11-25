# Gate 1 try: rebuilds the promoted/demoted lists at cut-off 0.1 from the 660 fine-mapping files,
# finds each variant's score-file row, and reads its 89 standardized scores.
import glob, os
import numpy as np, pandas as pd, gcsfs, h5py

FM = "decomposition/data/raw/sniff/finemap/"
IDX = "decomposition/data/interim/score_index.parquet"
SRC = "seqnn-share/sniff/borzoi_veps/scores_indel_snp_distmatched.h5"
OUT_SETS = "decomposition/data/interim/sets_cutoff01.csv"
OUT_SCORES = "decomposition/data/interim/fig_scores_long_01.parquet"
CUT = 0.1
vt = pd.read_csv("decomposition/data/frozen/variant_traits.csv")
label_of = vt.drop_duplicates("trait").set_index("trait")["label"]
tm = pd.read_csv("decomposition/data/frozen/track_map.csv").sort_values("col")

# 1. Differential PIP for every variant in every trait; keep |dPIP| > 0.1
keep = ["CHR", "SNP", "BP", "A1", "A2", "PIP"]
files = sorted(glob.glob(FM + "sniff/*.txt.gz"))
parts = []
for i, f in enumerate(files):
    s = pd.read_csv(f, sep="\t", usecols=keep, dtype={"A1": str, "A2": str})
    b = pd.read_csv(os.path.join(FM, "baseline", os.path.basename(f)), sep="\t", usecols=keep,
                    dtype={"A1": str, "A2": str})
    m = s.merge(b, on=["CHR", "SNP", "BP", "A1", "A2"], suffixes=("_sniff", "_base"))
    m["dPIP"] = m["PIP_sniff"] - m["PIP_base"]
    m = m[m["dPIP"].abs() > CUT].copy()
    m["trait"] = os.path.basename(f).rsplit(".", 3)[0]
    parts.append(m)
    if i % 30 == 0:
        print(f"  {i}/{len(files)} files read", flush=True)
d = pd.concat(parts, ignore_index=True)
d["label"] = d["trait"].map(label_of)
assert d["label"].notna().all(), "a file name did not match a known trait"
d["variant"] = d["CHR"].astype(str) + ":" + d["BP"].astype(str) + ":" + d["A1"] + ":" + d["A2"]
d["direction"] = np.where(d["dPIP"] > 0, "promoted", "demoted")
d["is_indel"] = (d["A1"].str.len() > 1) | (d["A2"].str.len() > 1)
for cut in [0.1, 0.25]:
    p = d.loc[d["dPIP"] > cut, "variant"].nunique()
    q = d.loc[d["dPIP"] < -cut, "variant"].nunique()
    print(f"cut-off {cut}: {p:,} promoted, {q:,} demoted (unique variants)")
print("paper:      0.1 -> 11,644 promoted, 10,055 demoted; 0.25 -> 3,160 and 2,579")

# 2. Find each variant's row in the score file (hg19 position + letters)
fs = gcsfs.GCSFileSystem(token="anon")
def fetch(starts, size):
    out = []
    for i in range(0, len(starts), 1000):
        s = starts[i:i + 1000]
        out += fs.cat_ranges([SRC] * len(s), s, [x + size for x in s])
    return out
with fs.open(SRC, "rb") as f, h5py.File(f, "r") as h:
    info = {k: (h[k].id.get_offset(), h[k].dtype.itemsize) for k in ["ref_allele", "alt_allele"]}
    ds = h["logD2"]
    base, item, ncol, dt = ds.id.get_offset(), ds.dtype.itemsize, ds.shape[1], ds.dtype
    first, last = int(tm["col"].min()), int(tm["col"].max())
    check_row = None

v = d.drop_duplicates("variant")[["variant", "CHR", "BP", "A1", "A2"]].copy()
idx = pd.read_parquet(IDX)
v["chr"] = ("chr" + v["CHR"].astype(str)).astype(idx["chr"].dtype)
cand = v.merge(idx, left_on=["chr", "BP"], right_on=["chr", "pos"])
def letters(k, rows):
    off, it = info[k]
    return [x.rstrip(b"\x00").decode() for x in fetch([int(off + int(r) * it) for r in rows], it)]
cand["s_ref"], cand["s_alt"] = letters("ref_allele", cand["row"]), letters("alt_allele", cand["row"])
ok = (((cand["s_ref"] == cand["A1"]) & (cand["s_alt"] == cand["A2"]))
      | ((cand["s_ref"] == cand["A2"]) & (cand["s_alt"] == cand["A1"])))
rows = cand[ok].drop_duplicates("variant")[["variant", "row"]].reset_index(drop=True)
print(f"score rows found for {len(rows):,} of {len(v):,} variants")

# 3. Read the 89 standardized scores for those rows (same byte method and self-test as before)
cols = tm["col"].to_numpy()
size = (last - first + 1) * item
raw = fetch([int(base + (int(r) * ncol + first) * item) for r in rows["row"]], size)
mat = np.vstack([np.frombuffer(b, dtype=dt) for b in raw])
with fs.open(SRC, "rb") as f, h5py.File(f, "r") as h:
    assert np.array_equal(mat[0], h["logD2"][int(rows["row"][0]), first:last + 1]), "byte arithmetic wrong"
mat = mat[:, cols - first].astype(np.float64)
long = (pd.DataFrame(mat, columns=cols).assign(variant=rows["variant"].values)
          .melt(id_vars="variant", var_name="col", value_name="score"))
long["col"] = long["col"].astype(int)
long.merge(tm[["col", "group"]], on="col").to_parquet(OUT_SCORES, index=False)
d[d["variant"].isin(rows["variant"])][["variant", "trait", "label", "dPIP", "direction", "is_indel"]] \
    .to_csv(OUT_SETS, index=False)
print(f"saved {OUT_SETS} and {OUT_SCORES}")