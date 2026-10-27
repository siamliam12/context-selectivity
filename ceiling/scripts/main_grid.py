from pathlib import Path
 
import duckdb
import matplotlib
matplotlib.use("Agg")  # draw straight to a file (the server has no screen)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import pearsonr, rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]         
ALIGNED = str(ROOT / "data/frozen/aligned.parquet")
OUT = ROOT / "results/ceiling"
P_THRESH = 1e-5

TISSUES = [
    ("Brain",      "GTEx_Brain",      "ROSMAP_Brain"),
    ("Blood",      "GTEx_Blood",      "Lepik_Blood"),
    ("Muscle",     "GTEx_Muscle",     "FUSION_Muscle"),
    ("Fat",        "GTEx_Fat",        "TWINSUK_Fat"),
    ("Skin",       "GTEx_Skin",       "TWINSUK_Skin"),
    ("Fibroblast", "GTEx_FibroBlast", "GENCORD_FibroBlast"),
]

LABELS = [t[0] for t in TISSUES]
ROWS = [t[1] for t in TISSUES]
COLS = [t[2] for t in TISSUES]
GRIDS = ["primary", "pearson", "spearman_z", "leadvar", "unselected"]

names = set(pq.read_schema(ALIGNED).names)
needed = [f"{s}_{x}" for s in ROWS + COLS for x in ("beta", "se")] + [f"{s}_pvalue" for s in ROWS]
missing = [c for c in needed if c not in names]
assert not missing, f"not in aligned.parquet (check spelling and capitals): {missing}"
OUT.mkdir(parents=True, exist_ok=True)
assert not (OUT / "grid_primary.tsv").exists(), "grids already exist; delete them on purpose to rebuild"

grids = {g: pd.DataFrame(np.nan, index=ROWS, columns=COLS) for g in GRIDS}
counts = {}

for r in ROWS:
    wanted = ", ".join(f"{s}_beta, {s}_se" for s in [r] + COLS)
    sel = duckdb.sql(f"""
        SELECT variant, gene, {r}_pvalue AS p, {wanted}
        FROM '{ALIGNED}'
        WHERE {r}_pvalue < {P_THRESH}
    """).df()

    lead = sel.sort_values(["gene", "p", "variant"]).drop_duplicates("gene")
 
    x = sel[f"{r}_beta"]
    xz = sel[f"{r}_beta"] / sel[f"{r}_se"]
    for c in COLS:
        y = sel[f"{c}_beta"]
        grids["primary"].loc[r, c] = spearmanr(x, y)[0]
        grids["pearson"].loc[r, c] = pearsonr(x, y)[0]
        grids["spearman_z"].loc[r, c] = spearmanr(xz, y / sel[f"{c}_se"])[0]
        grids["leadvar"].loc[r, c] = spearmanr(lead[f"{r}_beta"], lead[f"{c}_beta"])[0]
    
    counts[r] = {"filtered_pairs": len(sel), "lead_pairs": len(lead)}
    print(f"{r:16s} {len(sel):>9,} filtered pairs   {len(lead):>6,} lead pairs")

def ranks(column):
    """Read one column in file order and replace every value by its rank (ties share the average rank)."""
    values = pq.read_table(ALIGNED, columns=[column]).column(0).to_numpy()
    return rankdata(values)

print("\nunselected grid: ranking 50 million values per column (a few minutes)...")
col_ranks = {c: ranks(f"{c}_beta") for c in COLS}
for r in ROWS:
    row_ranks = ranks(f"{r}_beta")
    for c in COLS:
        # Spearman = ordinary (Pearson) correlation computed on ranks
        grids["unselected"].loc[r, c] = np.corrcoef(row_ranks, col_ranks[c])[0, 1]
    counts[r]["unselected_pairs"] = len(row_ranks)
    print(f"{r:16s} done")

for g in GRIDS:
    grids[g].to_csv(OUT / f"grid_{g}.tsv", sep="\t", float_format="%.6f")
n_pairs = pd.DataFrame.from_dict(counts, orient="index")
n_pairs.to_csv(OUT / "grid_n_pairs.tsv", sep="\t")
 
 
# ---------- 6. print every grid as a boxed table ----------
def show(title, df):
    t = df.astype(float).round(3)
    t.columns = LABELS                      # short names so the table fits the screen
    t.insert(0, "GTEx row", LABELS)
    print(f"\n{title}")
    duckdb.from_df(t.reset_index(drop=True)).show()
 
 
print("\nRows = GTEx tissue, columns = independent study. Diagonal = same tissue.")
for g in GRIDS:
    show(f"GRID: {g}", grids[g])
 
nt = n_pairs.copy()
nt.insert(0, "GTEx row", LABELS)
print("\nPAIRS USED PER ROW")
duckdb.from_df(nt.reset_index(drop=True)).show()
 
# ---------- 7. heatmap of the primary grid ----------
m = grids["primary"].to_numpy(dtype=float)
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(m, cmap="viridis")
ax.set_xticks(range(6), LABELS, rotation=45, ha="right")
ax.set_yticks(range(6), LABELS)
ax.set_xlabel("Independent study (column)")
ax.set_ylabel("GTEx (row)")
mid = (m.min() + m.max()) / 2
for i in range(6):
    for j in range(6):
        ax.text(j, i, f"{m[i, j]:.3f}", ha="center", va="center",
                color="black" if m[i, j] > mid else "white", fontsize=9)
    ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="red", linewidth=2))
fig.colorbar(im, ax=ax, label="Spearman correlation")
ax.set_title("Primary grid (same-tissue cells outlined in red)")
fig.tight_layout()
fig.savefig(OUT / "heatmap_primary.png", dpi=150)
print(f"\nsaved grids, pair counts and heatmap to {OUT}")