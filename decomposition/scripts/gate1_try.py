# Gate 1 debugging: tries candidate variant sets and methods; prints one line per try and logs them all.
import os
import numpy as np
import pandas as pd

FZ, INT = "decomposition/data/frozen/", "decomposition/data/interim/"
LOG = "decomposition/results/d3/attempts_detail.csv"
fig = pd.read_csv(FZ + "fig4c_labels.csv")
labels = sorted(pd.read_csv(FZ + "trait_map.csv")["label"])
vt = pd.read_csv(FZ + "variant_traits.csv")
sc = pd.read_csv(FZ + "set_c_variants.csv")
SCORES = {"dm": pd.read_parquet(INT + "fig_scores_long.parquet"),
          "raw": pd.read_parquet(INT + "fig_scores_long_raw.parquet")}
VALUES = [("Cholesterol", "liver", 2.42), ("MCH", "bone_marrow", 1.08),
          ("MC", "blood", 1.48), ("Diastolic BP", "heart", 0.59)]
SETS = {"A": vt[vt["direction"] == "promoted"], "B": vt, "C": sc}
SCORES["dm01"] = pd.read_parquet(INT + "fig_scores_long_01.parquet")
s01 = pd.read_csv(INT + "sets_cutoff01.csv")
SETS.update({"A01": s01[s01["direction"] == "promoted"], "B01": s01})
SETS.update({k + "_snp": v[~v["is_indel"].astype(bool)] for k, v in list(SETS.items())})

def val(w, t, g):
    return w.loc[t, g] if (t in w.index and g in w.columns) else np.nan

def table(score, members, ddof=0, agg="median", drop=(), order="group_then_z", within="mean",weight=None,unit=None):
    s = SCORES[score]
    s = s[s["variant"].isin(members["variant"]) & ~s["group"].isin(drop)]
    if order == "group_then_z":              # average tracks per group, then Z across groups
        gm = s.groupby(["variant", "group"])["score"].agg(within).unstack("group")
        z = gm.sub(gm.mean(axis=1), axis=0).div(gm.std(axis=1, ddof=ddof), axis=0)
    else:                                    # Z across the 89 tracks, then average Z per group
        tr = s.pivot_table(index="variant", columns="col", values="score")
        zt = tr.sub(tr.mean(axis=1), axis=0).div(tr.std(axis=1, ddof=ddof), axis=0)
        col2g = s.drop_duplicates("col").set_index("col")["group"]
        z = zt.T.groupby(col2g).mean().T
    zl = z.stack().rename("z").reset_index()
    keep_cols = ["variant", "label"] + ([weight] if weight else []) + ([unit] if unit else [])
    m = members[keep_cols].merge(zl, on="variant")
    if unit:
        m = m.groupby(["label", unit, "group"], as_index=False)["z"].mean()
    if weight:
        m["w"] = m[weight].abs()
        m["wz"] = m["z"] * m["w"]
        g = m.groupby(["label", "group"])
        return (g["wz"].sum() / g["w"].sum()).unstack("group")
    return m.groupby(["label", "group"])["z"].agg(agg).unstack("group")

def check(w):
    v = sum(abs(val(w, t, g) - x) <= 0.02 for t, g, x in VALUES)
    d = sum(abs(val(w, r.label, r.group) - r.approx_z) <= 0.05 for r in fig.itertuples())
    bad = 0
    for t in labels:
        row = w.loc[t]
        ignore = set(row[(row - 0.5).abs() <= 0.02].index) | set(fig[(fig.label == t) & fig.on_line].group)
        bad += (set(row[row > 0.5].index) - ignore) != (set(fig[fig.label == t].group) - ignore)
    return v, d, 15 - bad
TABLES = {}
def run(name, set_name, score="dm", **opts):
    members = SETS[set_name]
    w = table(score, members, **opts)
    TABLES[(name, set_name, score)] = w
    v, d, t = check(w)
    rmse = np.sqrt(np.nanmean([(val(w, r.label, r.group) - r.approx_z) ** 2 for r in fig.itertuples()]))
    row = {"try": name, "set": set_name, "scores": score, "values": f"{v}/4", "dots": f"{d}/17",
           "traits": f"{t}/15", "dot_dist": round(rmse, 3),**{f"{a}-{g}": round(val(w, a, g), 2) for a, g, _ in VALUES}}
    pd.DataFrame([row]).to_csv(LOG, mode="a", index=False, header=not os.path.exists(LOG))
    return row

smp = np.load(INT + "transform_sample.npz")
track_cols = pd.read_csv(FZ + "track_map.csv").sort_values("col")["col"].to_numpy()
raw_long = SCORES["raw"]
def restandardize(f, snp_only):
    X = smp["X"].astype(float)
    X = f(X[~smp["is_indel"]] if snp_only else X)
    mu = pd.Series(X.mean(axis=0), index=track_cols)
    sd = pd.Series(X.std(axis=0), index=track_cols)
    return raw_long.assign(score=(f(raw_long["score"].to_numpy()) - raw_long["col"].map(mu)) / raw_long["col"].map(sd))
SCORES["std_all"] = restandardize(lambda x: x, False)
SCORES["std_snp"] = restandardize(lambda x: x, True)
SCORES["logstd_all"] = restandardize(np.log1p, False)
SCORES["log_raw"] = raw_long.assign(score=np.log1p(raw_long["score"]))

c01 = pd.read_csv(INT + "sets_cutoff01_cs.csv")
c01["in_cs_sniff"] = c01["in_cs_sniff"].astype(bool)
tries = []
for t in [0.1, 0.15, 0.2, 0.25]:
    big = c01[c01["dPIP"].abs() > t]
    SETS[f"B{t}_cs"] = big[big["in_cs_sniff"]]
    SETS[f"A{t}_cs"] = big[big["in_cs_sniff"] & big["direction"].eq("promoted")]
    for s in [f"B{t}_cs", f"A{t}_cs"]:
        tries.append(run("mean", s, "dm01", agg="mean"))
        tries.append(run("mean over credible sets", s, "dm01", agg="mean", unit="cs_sniff"))
out = pd.DataFrame(tries).sort_values("dot_dist")
print("printed values: Cholesterol-liver 2.42, MCH-bone_marrow 1.08, MC-blood 1.48, Diastolic BP-heart 0.59")
print(out.to_string(index=False))

top = out.iloc[0]
w = TABLES[(top["try"], top["set"], top["scores"])]
print(f"\nbest try: {top['try']}, set {top['set']}")
print("\nthe 17 dots")
for r in fig.itertuples():
    print(f"  {r.label:16s} {r.group:12s} {w.loc[r.label, r.group]:6.2f} vs {r.approx_z:.2f}  gap {w.loc[r.label, r.group] - r.approx_z:+.2f}")
print("\ntraits whose labels differ (ours above 0.5 / figure's):")
for t in labels:
    ours, theirs = set(w.columns[w.loc[t] > 0.5]), set(fig[fig.label == t].group)
    if ours != theirs:
        print(f"  {t}: ours {sorted(ours)} | figure {sorted(theirs)}")