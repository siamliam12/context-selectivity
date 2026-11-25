# Gate 1 debugging: shows the gap between our two best tries and the figure, dot by dot and trait by trait.
import pandas as pd
import gate1_try as g   # reuses the table and sets (the try script's own runs print first)

best = {"mean A 0.25": g.table("dm", g.SETS["A"], agg="mean"),
        "mean A 0.1 SNPs": g.table("dm01", g.SETS["A01_snp"], agg="mean")}
rows = []
for r in g.fig.itertuples():
    row = {"trait": r.label, "tissue": r.group, "figure": r.approx_z}
    for k, w in best.items():
        row[k] = round(w.loc[r.label, r.group], 2)
        row[k + " gap"] = round(w.loc[r.label, r.group] - r.approx_z, 2)
    rows.append(row)
print(pd.DataFrame(rows).to_string(index=False))

print("\nspread of each trait's unlabelled tissues (min to max), to compare with the grey dots in the figure:")
lab = g.fig.groupby("label")["group"].apply(set)
sp = []
for t in g.labels:
    row = {"trait": t}
    for k, w in best.items():
        vals = w.loc[t].drop(list(lab.get(t, set())))
        row[k] = f"{vals.min():.2f} to {vals.max():.2f}"
    sp.append(row)
print(pd.DataFrame(sp).to_string(index=False))