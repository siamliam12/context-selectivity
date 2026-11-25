# D2 + D3: rebuilds Fig. 4C six ways (sets A, B, C x SD by n and by n-1) and runs the frozen Gate 1 check.
import hashlib, os
import pandas as pd
from table_fn import build_table

FZ = "decomposition/data/frozen/"
import sys
TAG = sys.argv[1] if len(sys.argv) > 1 else "dm"
OUT2, OUT3 = f"decomposition/results/d2/{TAG}/", f"decomposition/results/d3/{TAG}/"
os.makedirs(OUT2, exist_ok=True)
os.makedirs(OUT3, exist_ok=True)

# 0. The frozen files must be unchanged since the freeze
known = dict(line.split()[::-1] for line in open(FZ + "frozen.sha256") if line.strip())
for f in ["variant_traits.csv", "set_c_variants.csv", "fig4c_labels.csv", "track_map.csv", "trait_map.csv"]:
    assert known[f] == hashlib.sha256(open(FZ + f, "rb").read()).hexdigest(), f"{f} changed since the freeze"
print("frozen files match their fingerprints\n")

long = pd.read_parquet(f"decomposition/data/interim/fig_scores_long{'' if TAG == 'dm' else '_' + TAG}.parquet")
vt = pd.read_csv(FZ + "variant_traits.csv")
sc = pd.read_csv(FZ + "set_c_variants.csv")
fig = pd.read_csv(FZ + "fig4c_labels.csv")
labels = sorted(pd.read_csv(FZ + "trait_map.csv")["label"])
SETS = {"A": vt[vt["direction"] == "promoted"][["variant", "label"]],
        "B": vt[["variant", "label"]],
        "C": sc[["variant", "label"]]}
VALUES = [("Cholesterol", "liver", 2.42), ("MCH", "bone_marrow", 1.08),
          ("MC", "blood", 1.48), ("Diastolic BP", "heart", 0.59)]

def check(w):
    v_hit = [abs(w.loc[t, g] - x) <= 0.02 for t, g, x in VALUES]
    d_hit = [abs(w.loc[r.label, r.group] - r.approx_z) <= 0.05 for r in fig.itertuples()]
    bad = []
    for t in labels:
        row = w.loc[t]
        ignore = set(row[(row - 0.5).abs() <= 0.02].index) | set(fig[(fig.label == t) & fig.on_line].group)
        ours = set(row[row > 0.5].index) - ignore
        theirs = set(fig[fig.label == t].group) - ignore
        if ours != theirs:
            bad.append((t, sorted(ours - theirs), sorted(theirs - ours)))
    return v_hit, d_hit, bad

def verdict(nv, dots_missed, traits_missed):
    if nv == 4 and dots_missed <= 2 and traits_missed <= 3:
        return "Exact" if dots_missed == 0 and traits_missed == 0 else "Partial"
    if nv == 3 and dots_missed == 0 and traits_missed == 0:
        return "Partial"
    return "Unresolved"

rows, details = [], {}
for s, mem in SETS.items():
    for ddof in [0, 1]:
        tbl, dropped = build_table(long, mem, ddof=ddof)
        tbl.to_csv(f"{OUT2}table_set{s}_ddof{ddof}.csv", index=False)
        w = tbl.pivot(index="label", columns="group", values="median_z")
        v_hit, d_hit, bad = check(w)
        rows.append({"set": s, "sd": "n" if ddof == 0 else "n-1", "variants": mem["variant"].nunique(),
                     "dropped": dropped, "values_hit": f"{sum(v_hit)}/4", "dots_hit": f"{sum(d_hit)}/17",
                     "traits_match": f"{15 - len(bad)}/15",
                     "verdict": verdict(sum(v_hit), 17 - sum(d_hit), len(bad))})
        details[(s, ddof)] = (w, v_hit, d_hit, bad)

summ = pd.DataFrame(rows)
print(summ.to_string(index=False))
summ.to_csv(OUT3 + "gate1_summary.csv", index=False)

rank = summ["verdict"].map({"Exact": 0, "Partial": 1, "Unresolved": 2})
best = summ.assign(rank=rank).sort_values(["rank", "set", "sd"], kind="stable").iloc[0]
print(f"\nGATE 1: {best['verdict']} (set {best['set']}, SD by {best['sd']})")

w, v_hit, d_hit, bad = details[(best["set"], 0 if best["sd"] == "n" else 1)]
print("\nfour printed values (ours vs printed):")
for (t, g, x), ok in zip(VALUES, v_hit):
    print(f"  {t:13s} {g:12s} {w.loc[t, g]:6.3f} vs {x:.2f}  {'ok' if ok else 'MISS'}")
print("\n17 figure dots (ours vs figure reading):")
for r, ok in zip(fig.itertuples(), d_hit):
    print(f"  {r.label:16s} {r.group:12s} {w.loc[r.label, r.group]:6.3f} vs {r.approx_z:.2f}  {'ok' if ok else 'MISS'}")
print("\ntraits whose labelled tissues differ:", bad if bad else "none")
print("\ntext vs figure (text: HLSRC-bone marrow 0.78, MC-spleen 0.50):")
for t, g in [("HLSRC", "bone_marrow"), ("MPV", "bone_marrow"), ("MC", "spleen"), ("MPV", "spleen"), ("MC", "lung")]:
    print(f"  {t:6s} {g:12s} {w.loc[t, g]:6.3f}")

summ.assign(attempt=TAG).to_csv("decomposition/results/d3/attempts.csv", mode="a", index=False, header=not os.path.exists("decomposition/results/d3/attempts.csv"))