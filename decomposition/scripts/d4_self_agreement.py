# D4: self-agreement of the published-style table (reconstructed method, Amendment 1) at half size,
# with trait x LD-block units. Saves the 1,000 splits so Part E reuses exactly the same ones.
import os
import numpy as np
import pandas as pd
from table_fn import group_means, variant_z

SEED, N_SPLITS = 20261116 + 4, 1000          # D4's seed, derived from the master seed
OUT = "decomposition/results/d4/"
os.makedirs(OUT, exist_ok=True)

long = pd.read_parquet("decomposition/data/interim/fig_scores_long.parquet")
rec = pd.read_csv("decomposition/data/frozen/fig4c_reconstructed_set.csv")
bl = pd.read_csv("decomposition/data/raw/ldblocks/EUR_fourier_ls-all.bed", sep=r"\s+")

# 1. Each variant's LD block, from its hg19 position (blocks are start <= position < stop)
rec["CHR"] = rec["variant"].str.split(":").str[0].astype(int)
rec["BP"] = rec["variant"].str.split(":").str[1].astype(int)
rec["block"] = None
for c, b in bl.groupby("chr"):
    m = rec["CHR"].eq(int(c.replace("chr", "")))
    pos, starts, stops = rec.loc[m, "BP"].to_numpy(), b["start"].to_numpy(), b["stop"].to_numpy()
    i = np.clip(np.searchsorted(starts, pos, side="right") - 1, 0, None)
    ok = (pos >= starts[i]) & (pos < stops[i])
    rec.loc[m, "block"] = np.where(ok, c + ":" + starts[i].astype(str), None)
print(f"variants outside every block: {int(rec['block'].isna().sum())} (each becomes its own unit)")
rec["block"] = rec["block"].fillna("none:" + rec["variant"])
rec["unit"] = rec["label"] + "|" + rec["block"]

# 2. Per-variant Z once (it does not depend on the split)
z, dropped = variant_z(group_means(long[long["variant"].isin(rec["variant"])]), ddof=0)
rec = rec[rec["variant"].isin(z.index)].reset_index(drop=True)
Z = z.loc[rec["variant"]].to_numpy()
traits = sorted(rec["label"].unique())
t_idx = pd.Categorical(rec["label"], categories=traits).codes
units = rec.groupby("label")["unit"].unique()
print(f"variant-trait rows: {len(rec):,}; zero-spread variants dropped: {dropped}")
print("units per trait:", {t: len(units[t]) for t in traits})

def trait_means(mask):
    sums, cnt = np.zeros((len(traits), Z.shape[1])), np.zeros(len(traits))
    np.add.at(sums, t_idx[mask], Z[mask])
    np.add.at(cnt, t_idx[mask], 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return sums / cnt[:, None]

def overlap(a, b):
    if np.isnan(a).any() or np.isnan(b).any():
        return np.nan
    A, B = a > 0.5, b > 0.5
    either = (A | B).sum()
    return 1.0 if either == 0 else (A & B).sum() / either

# 3. The 1,000 splits: each trait's units shuffled and cut in half
rng = np.random.default_rng(SEED)
rows = []
for s in range(N_SPLITS):
    for t in traits:
        perm = rng.permutation(units[t])
        k = len(perm) // 2
        rows += [(s, u, 1) for u in perm[:k]] + [(s, u, 2) for u in perm[k:]]
splits = pd.DataFrame(rows, columns=["split", "unit", "half"])
splits.to_parquet(OUT + "splits.parquet", index=False)

# 4. Self-agreement for every split
unit_row = rec["unit"].to_numpy()
res = []
for s, sp in splits.groupby("split"):
    h = dict(zip(sp["unit"], sp["half"]))
    half = np.array([h[u] for u in unit_row])
    T1, T2 = trait_means(half == 1), trait_means(half == 2)
    for i, t in enumerate(traits):
        r = np.corrcoef(T1[i], T2[i])[0, 1] if not (np.isnan(T1[i]).any() or np.isnan(T2[i]).any()) else np.nan
        res.append((s, t, overlap(T1[i], T2[i]), r, 2 * r / (1 + r)))
res = pd.DataFrame(res, columns=["split", "trait", "overlap", "r", "r_sb"])
res.to_csv(OUT + "self_agreement_splits.csv", index=False)

q = lambda x: f"{x.median():.2f} [{x.quantile(0.025):.2f}, {x.quantile(0.975):.2f}]"
summ = res.groupby("trait").agg(overlap=("overlap", q), profile_r_SB=("r_sb", q))
summ.insert(0, "units", [len(units[t]) for t in summ.index])
print("\nself-agreement at half size: median [95% interval]")
print(summ.to_string())
summ.to_csv(OUT + "self_agreement_summary.csv")