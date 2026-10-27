from itertools import permutations
from pathlib import Path
import duckdb
import matplotlib
matplotlib.use("Agg")  
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/ceiling"
GRIDS = ["primary", "pearson", "spearman_z", "leadvar", "unselected"]

grids = {g: pd.read_csv(RES / f"grid_{g}.tsv", sep="\t", index_col=0) for g in GRIDS}
rows, cols = list(grids["primary"].index), list(grids["primary"].columns)
for g, df in grids.items():
    assert list(df.index) == rows and list(df.columns) == cols, f"{g}: labels differ from primary"
    assert df.shape == (6, 6) and not df.isna().any().any(), f"{g}: not a complete 6x6 grid"
labels = [r.split("_", 1)[1] for r in rows]          # "GTEx_Brain" -> "Brain"
 
print("Diagonal pairs (row i with column i) - each line must be the same tissue:")
for r, c in zip(rows, cols):
    print(f"  {r:18s} <-> {c}")

def parts(A):
    """Diagonal mean, off-diagonal mean, and G for a square grid."""
    n = A.shape[0]
    diag = np.trace(A) / n
    off = (A.sum() - np.trace(A)) / (n * n - n)
    return diag, off, diag - off
 
summary = pd.DataFrame(
    [dict(zip(["grid", "diagonal_mean", "off_diagonal_mean", "G"], [g, *parts(df.to_numpy())]))
     for g, df in grids.items()]
)
main3 = summary.set_index("grid").loc[["primary", "pearson", "spearman_z"], "G"]
signs_agree = bool((np.sign(main3) == np.sign(main3.iloc[0])).all())

A = grids["primary"].to_numpy()
per = pd.DataFrame({
    "tissue": labels,
    "diagonal": np.diag(A),
    "row_off_mean": [(A[i].sum() - A[i, i]) / 5 for i in range(6)],
    "col_off_mean": [(A[:, i].sum() - A[i, i]) / 5 for i in range(6)],
})
per["row_gap"] = per["diagonal"] - per["row_off_mean"]
G_primary = parts(A)[2]
assert abs(per["row_gap"].mean() - G_primary) < 1e-12, "row gaps do not average to G - arithmetic bug"

def shuffle_test(M):
    """G for every way of reassigning the column labels. The real labelling is included."""
    perms = list(permutations(range(M.shape[0])))
    null = np.array([parts(M[:, list(p)])[2] for p in perms])
    G_real = parts(M)[2]
    n_at_least = int((null >= G_real).sum())
    return G_real, perms, null, n_at_least
 
fib = next(i for i, lab in enumerate(labels) if lab.lower().startswith("fibro"))
keep = [i for i in range(6) if i != fib]
A5 = A[np.ix_(keep, keep)]
 
G6, perms6, null6, k6 = shuffle_test(A)
G5, perms5, null5, k5 = shuffle_test(A5)
 
shuf = pd.DataFrame([
    {"test": "6x6 primary", "G": G6, "shuffles": len(null6), "shuffles_G_at_least_real": k6,
     "p_fraction": f"{k6}/{len(null6)}", "p": k6 / len(null6),
     "null_mean": null6.mean(), "null_sd": null6.std(), "z": (G6 - null6.mean()) / null6.std()},
    {"test": "5x5 no fibroblast", "G": G5, "shuffles": len(null5), "shuffles_G_at_least_real": k5,
     "p_fraction": f"{k5}/{len(null5)}", "p": k5 / len(null5),
     "null_mean": null5.mean(), "null_sd": null5.std(), "z": (G5 - null5.mean()) / null5.std()},
])

summary.to_csv(RES / "g_summary.tsv", sep="\t", index=False, float_format="%.6f")
per.to_csv(RES / "per_tissue_primary.tsv", sep="\t", index=False, float_format="%.6f")
shuf.to_csv(RES / "shuffle_test.tsv", sep="\t", index=False, float_format="%.6f")
pd.DataFrame({"column_order": [",".join(map(str, p)) for p in perms6], "G": null6}).to_csv(
    RES / "permutation_null_6x6.tsv", sep="\t", index=False, float_format="%.6f")
pd.DataFrame({"column_order": [",".join(map(str, p)) for p in perms5], "G": null5}).to_csv(
    RES / "permutation_null_5x5.tsv", sep="\t", index=False, float_format="%.6f")
 
 
# ---------- 7. print as boxed tables ----------
def show(title, df, decimals=4):
    print(f"\n{title}")
    duckdb.from_df(df.round(decimals)).show()
 
 
show("G FOR EVERY GRID", summary)
print(f"Sign check (primary, pearson, spearman_z): {'all agree' if signs_agree else 'DISAGREE - no headline number'}")
show("PER-TISSUE TABLE (primary grid)", per)
show("SHUFFLE TEST", shuf)
 
# ---------- 8. histogram of the shuffles with the real G marked ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
panels = [(axes[0], null6, G6, "6x6 grid: 720 shuffles"),
          (axes[1], null5, G5, "5x5 without fibroblast: 120 shuffles")]
for ax, null, G, title in panels:
    ax.hist(null, bins=40, color="grey")
    ax.axvline(G, color="red", linewidth=2, label=f"real G = {G:.4f}")
    ax.set_title(title)
    ax.set_xlabel("G")
    ax.set_ylabel("number of shuffles")
    ax.legend()
fig.tight_layout()
fig.savefig(RES / "permutation_hist.png", dpi=150)
print(f"\nsaved tables and histogram to {RES}")