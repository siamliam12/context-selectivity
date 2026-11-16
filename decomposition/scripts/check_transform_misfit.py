# C3 step 3b: looks closer at the SNP conversion for tracks below R2 0.999.
# Re-reads the same 20,000 random rows (same seed as step 3), saves them for reuse, and checks
# whether the misfit comes from a few extreme variants or from a curve. Not our variants.
import os
import gcsfs, h5py
import numpy as np
import pandas as pd

MAP = "decomposition/data/frozen/track_map.csv"
BASE = "seqnn-share/sniff/borzoi_veps/"
RAW, DM = BASE + "scores.h5", BASE + "scores_indel_snp_distmatched.h5"
SAMPLE = "decomposition/data/interim/transform_sample.npz"
OUT = "decomposition/results/c3/transform_misfit.csv"
SEED, N, N_ROWS, DROP = 20261111, 20_000, 19_534_182, 20

fs = gcsfs.GCSFileSystem(token="anon")
tm = pd.read_csv(MAP).sort_values("col").reset_index(drop=True)
cols = tm["col"].to_numpy()
first, last = int(cols.min()), int(cols.max())
rows = np.sort(np.random.default_rng(SEED).choice(N_ROWS, size=N, replace=False))

def fetch(path, starts, size):
    out = []
    for i in range(0, len(starts), 1000):
        s = starts[i:i + 1000]
        out += fs.cat_ranges([path] * len(s), s, [x + size for x in s])
    return out

def scores(path):
    with fs.open(path, "rb") as f, h5py.File(f, "r") as h:
        ds = h["logD2"]
        base, item, ncol, dt = ds.id.get_offset(), ds.dtype.itemsize, ds.shape[1], ds.dtype
    size = (last - first + 1) * item
    mat = np.vstack([np.frombuffer(b, dtype=dt) for b in
                     fetch(path, [int(base + (int(r) * ncol + first) * item) for r in rows], size)])
    return mat[:, cols - first].astype(np.float32)

def letters(name):
    with fs.open(DM, "rb") as f, h5py.File(f, "r") as h:
        ds = h[name]
        base, item = ds.id.get_offset(), ds.dtype.itemsize
    return [b.rstrip(b"\x00") for b in fetch(DM, [int(base + int(r) * item) for r in rows], item)]

if os.path.exists(SAMPLE):
    z = np.load(SAMPLE)
    X, Y, is_indel = z["X"], z["Y"], z["is_indel"]
    assert np.array_equal(z["rows"], rows), "saved sample has different rows"
    print("loaded saved sample")
else:
    refs, alts = letters("ref_allele"), letters("alt_allele")
    is_indel = np.array([len(a) != 1 or len(b) != 1 for a, b in zip(refs, alts)])
    X, Y = scores(RAW), scores(DM)
    np.savez(SAMPLE, rows=rows, X=X, Y=Y, is_indel=is_indel)
    print(f"read and saved {SAMPLE}")

def r2(y, yhat):
    return 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)

snp = ~is_indel
recs = []
for j, row in tm.iterrows():
    x, y = X[snp, j].astype(np.float64), Y[snp, j].astype(np.float64)
    a, b = np.polyfit(x, y, 1)
    res = y - (a * x + b)
    keep = np.argsort(np.abs(res))[:-DROP]                      # drop the 20 worst-fitting SNPs
    a2, b2 = np.polyfit(x[keep], y[keep], 1)
    q = np.polyfit(x, y, 2)                                     # a curve instead of a line
    worst = np.argmax(np.abs(res))
    recs.append({"col": row["col"], "group": row["group"],
                 "r2_line": r2(y, a * x + b),
                 "r2_line_without_20": r2(y[keep], a2 * x[keep] + b2),
                 "r2_curve": r2(y, np.polyval(q, x)),
                 "resid_sd": res.std(), "max_abs_resid": np.abs(res).max(),
                 "worst_raw": x[worst], "worst_raw_percentile": (x < x[worst]).mean() * 100})
d = pd.DataFrame(recs)

print(f"\nSNP tracks with R2 >= 0.999: straight line {int((d['r2_line'] >= 0.999).sum())}, "
      f"line without the 20 worst SNPs {int((d['r2_line_without_20'] >= 0.999).sum())}, "
      f"curve {int((d['r2_curve'] >= 0.999).sum())} (of 89)")
print("\nthe 5 worst tracks:")
print(d.nsmallest(5, "r2_line").to_string(index=False, float_format=lambda v: f"{v:.5f}"))

d.to_csv(OUT, index=False)
print(f"\nsaved {OUT}")