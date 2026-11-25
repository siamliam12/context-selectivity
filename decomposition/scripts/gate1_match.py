# Gate 1 debugging: for each figure column, which of our 15 trait rows matches it best?
import numpy as np
import pandas as pd
import gate1_try as g   # reuses the sets and the table (its own tries print first)

w = g.table("dm", g.SETS["B_cs"], agg="mean")
fig = g.fig

def distance(fig_trait, our_trait):
    lab = fig[fig.label == fig_trait].set_index("group")["approx_z"]
    terms = []
    for grp in w.columns:
        ours = w.loc[our_trait, grp]
        if grp in lab.index:
            terms.append((ours - lab[grp]) ** 2)          # labelled dot: how far off
        else:
            terms.append(max(0.0, ours - 0.5) ** 2)       # grey in figure: penalty only if ours is above 0.5
    return np.sqrt(np.mean(terms))

rows = []
for f in g.labels:
    d = pd.Series({o: distance(f, o) for o in g.labels}).sort_values()
    rows.append({"figure column": f, "own trait": round(d[f], 3),
                 "best match": d.index[0], "best": round(d.iloc[0], 3),
                 "2nd match": d.index[1], "2nd": round(d.iloc[1], 3)})
print("\nwhich of our traits best matches each figure column (smaller = closer):")
print(pd.DataFrame(rows).to_string(index=False))