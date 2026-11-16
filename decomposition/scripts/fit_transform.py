# C3 step 3: checks how raw scores (scores.h5) turn into standardized scores (distmatched file)
# for the 89 Fig. 4C tracks, on random rows of the whole file (not our variants).
# Fits standardized = a x raw + b per track, separately for SNPs and indels.
import os
import gcsfs, h5py
import numpy as np
import pandas as pd

MAP = "decomposition/data/frozen/track_map.csv"
BASE = "seqnn-share/sniff/borzoi_veps/"
RAW, DM = BASE + "scores.h5", BASE + "scores_indel_snp_distmatched.h5"
OUT = "decomposition/results/c3/transform_fit.csv"
SEED, N, N_ROWS = 20261111, 20_000, 19_534_182

fs = gcsfs.GCSFileSystem(token="anon")
tm = pd.read_csv(MAP).sort_values("col").reset_index(drop=True)
cols = tm["col"].to_numpy()
first, last = int(cols.min()), int(cols.max())
print(f"the 89 tracks are score columns {first} to {last}; side by side: {last - first + 1 == len(cols)}")

rng = np.random.default_rng(SEED)
rows = np.sort(rng.choice(N_ROWS, size=N, replace=False))

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
        check = ds[0:3, first:last + 1]                      # read by h5py, for the self-test
    size = (last - first + 1) * item
    start = lambda r: int(base + (int(r) * ncol + first) * item)
    test = np.vstack([np.frombuffer(b, dtype=dt) for b in fetch(path, [start(r) for r in range(3)], size)])
    assert np.array_equal(test, check), "byte arithmetic does not match h5py"
    mat = np.vstack([np.frombuffer(b, dtype=dt) for b in fetch(path, [start(r) for r in rows], size)])
    return mat[:, cols - first].astype(np.float64)

def letters(name):
    with fs.open(DM, "rb") as f, h5py.File(f, "r") as h:
        ds = h[name]
        base, item = ds.id.get_offset(), ds.dtype.itemsize
    return [b.rstrip(b"\x00") for b in fetch(DM, [int(base + int(r) * item) for r in rows], item)]

refs, alts = letters("ref_allele"), letters("alt_allele")
is_indel = np.array([len(a) != 1 or len(b) != 1 for a, b in zip(refs, alts)])
print(f"random rows: {N:,} (seed {SEED}); SNPs {int((~is_indel).sum()):,}, indels {int(is_indel.sum()):,}")

X, Y = scores(RAW), scores(DM)
print("self-tests passed; scores read for both files")

recs = []
for j, row in tm.iterrows():
    rec = {"col": row["col"], "group": row["group"]}
    for kind, m in [("snp", ~is_indel), ("indel", is_indel)]:
        x, y = X[m, j], Y[m, j]
        a, b = np.polyfit(x, y, 1)
        r2 = 1 - np.sum((y - (a * x + b)) ** 2) / np.sum((y - y.mean()) ** 2)
        rho = pd.Series(x).corr(pd.Series(y), method="spearman")
        rec.update({f"a_{kind}": a, f"b_{kind}": b, f"r2_{kind}": r2, f"spearman_{kind}": rho,
                    f"n_{kind}": int(m.sum())})
    rec["dm_mean"], rec["dm_sd"] = Y[:, j].mean(), Y[:, j].std()
    recs.append(rec)
fit = pd.DataFrame(recs)

for kind in ["snp", "indel"]:
    r2 = fit[f"r2_{kind}"]
    print(f"\n{kind}: R2 min {r2.min():.5f}, median {r2.median():.5f}; tracks with R2 >= 0.999: "
          f"{int((r2 >= 0.999).sum())} of 89; Spearman min {fit[f'spearman_{kind}'].min():.5f}")
print(f"\nslope a, SNP vs indel (median): {fit['a_snp'].median():.4f} vs {fit['a_indel'].median():.4f}")
print(f"intercept b, SNP vs indel (median): {fit['b_snp'].median():.4f} vs {fit['b_indel'].median():.4f}")
print(f"standardized scores in this sample: mean {fit['dm_mean'].median():.3f}, sd {fit['dm_sd'].median():.3f} (median over tracks)")
print("\nworst 5 tracks by SNP R2:")
print(fit.nsmallest(5, "r2_snp")[["col", "group", "r2_snp", "r2_indel", "a_snp", "a_indel"]].to_string(index=False))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fit.to_csv(OUT, index=False)
print(f"\nsaved {OUT}")