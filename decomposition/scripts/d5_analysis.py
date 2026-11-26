# D5 step 3: the genome-wide check. Does each variant's 31-group score profile follow its own locus's
# expression profile beyond chance? Main: window genes, GTEx v8. Checks: nearest gene; GTEx v10.
import os
import numpy as np, pandas as pd
from scipy.stats import rankdata
from table_fn import group_means

SEED, ROUNDS, BOOT = 20261116 + 52, 1000, 1000           # D5 analysis seed, from the master seed
OUT = "decomposition/results/d5/"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(SEED)

d5 = pd.read_csv("decomposition/data/interim/d5_sample_variants.csv").drop_duplicates("variant")
loc = pd.read_parquet("decomposition/data/interim/locus_expression.parquet")
loc = loc[loc["source"] == "d5"]
sc = group_means(pd.read_parquet("decomposition/data/interim/d5_scores_long.parquet"))
groups = list(sc.columns)

# LD block of each variant (hg19), for the bootstrap
bl = pd.read_csv("decomposition/data/raw/ldblocks/EUR_fourier_ls-all.bed", sep=r"\s+")
d5["block"] = None
for c, b in bl.groupby("chr"):
    m = d5["chr"].eq(c)
    pos, starts, stops = d5.loc[m, "pos"].to_numpy(), b["start"].to_numpy(), b["stop"].to_numpy()
    i = np.clip(np.searchsorted(starts, pos, side="right") - 1, 0, None)
    d5.loc[m, "block"] = np.where((pos >= starts[i]) & (pos < stops[i]), c + ":" + starts[i].astype(str), None)
d5["block"] = d5["block"].fillna("none:" + d5["variant"])
d5 = d5.set_index("variant")

def std_ranks(M):
    R = rankdata(M, axis=1)
    R = R - R.mean(axis=1, keepdims=True)
    return R / np.sqrt((R ** 2).sum(axis=1, keepdims=True))

def run(version, genes):
    e = loc[(loc["version"] == version) & (loc["genes"] == genes)].drop_duplicates("variant").set_index("variant")[groups]
    raw_spread = e.std(axis=1)
    le = np.log1p(e)
    level = le.mean(axis=1)                                   # overall expression level, for matching the null
    le = (le - le.mean()) / le.std(ddof=0)                    # each group standardized across the D5 sample
    common = d5.index.intersection(le.index).intersection(sc.index)
    S, X = sc.loc[common].to_numpy(), le.loc[common].to_numpy()
    keep = (S.std(axis=1) > 0) & (raw_spread.loc[common].to_numpy() > 0)
    info = d5.loc[common][keep]
    S, X = S[keep], X[keep]
    Rs, Rx = std_ranks(S), std_ranks(X)
    real = (Rs * Rx).sum(axis=1)                              # Spearman correlation, one per variant
    decile = pd.qcut(level.loc[common][keep].rank(method="first"), 10, labels=False).to_numpy()
    chrom = info["chr"].to_numpy()

    def null_perm():
        perm = np.arange(len(real))
        for d in range(10):
            idx = np.where(decile == d)[0]
            p = rng.permutation(idx)
            for _ in range(20):                               # no pairing within the same chromosome
                bad = chrom[idx] == chrom[p]
                if not bad.any():
                    break
                p[bad] = rng.permutation(p[bad])
            perm[idx] = p
        return perm

    perms = [null_perm() for _ in range(ROUNDS)]
    top_s, top_x = S.argmax(axis=1), X.argmax(axis=1)
    ub, inv = np.unique(info["block"].to_numpy(), return_inverse=True)
    members = [np.where(inv == k)[0] for k in range(len(ub))]

    def summary(mask, name):
        nulls = [np.median((Rs[mask] * Rx[p][mask]).sum(axis=1)) for p in perms]
        diff = np.median(real[mask]) - np.mean(nulls)
        boots = []
        for b in range(BOOT):
            idx = np.concatenate([members[k] for k in rng.integers(0, len(ub), len(ub))])
            idx = idx[mask[idx]]
            p = perms[b % ROUNDS]
            boots.append(np.median(real[idx]) - np.median((Rs[idx] * Rx[p][idx]).sum(axis=1)))
        lo, hi = np.percentile(boots, [2.5, 97.5])
        top_real = (top_s[mask] == top_x[mask]).mean()
        top_null = np.mean([(top_s[mask] == top_x[p][mask]).mean() for p in perms[:200]])
        return {"version": version, "genes": genes, "subset": name, "variants": int(mask.sum()),
                "real_median": round(np.median(real[mask]), 4), "chance_median": round(np.mean(nulls), 4),
                "difference": round(diff, 4), "ci_low": round(lo, 4), "ci_high": round(hi, 4),
                "top_match": round(top_real, 4), "top_match_chance": round(top_null, 4)}

    rows = [summary(np.ones(len(real), bool), "all")]
    if (version, genes) == ("v8", "window"):
        maxabs = np.abs(S).max(axis=1)
        top10 = maxabs >= np.quantile(maxabs, 0.9)
        inside = loc[(loc["version"] == "v8") & (loc["genes"] == "window")].drop_duplicates("variant") \
                    .set_index("variant").loc[info.index, "inside_gene"].to_numpy()
        indel = info["is_indel"].astype(bool).to_numpy()
        for mask, name in [(top10, "top 10% scores"), (~top10, "other 90%"), (inside, "inside a gene"),
                           (~inside, "outside genes"), (~indel, "SNPs"), (indel, "indels")]:
            rows.append(summary(mask, name))
    print(f"{version}/{genes}: {len(real):,} variants used ({int((~keep).sum()):,} left out for zero spread)", flush=True)
    return rows

res = pd.DataFrame(run("v8", "window") + run("v8", "nearest") + run("v10", "window"))
res.to_csv(OUT + "d5_summary.csv", index=False)
pd.set_option("display.width", 250)
print(res.to_string(index=False))