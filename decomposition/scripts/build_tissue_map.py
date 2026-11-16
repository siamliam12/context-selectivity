# C3 step 4d: freezes the sub-tissue map (which GTEx sub-tissues make up each of the 31 groups,
# with weights matching the tracks), and stages the group expression tables for v8 and v10.
# Gene expression only. No variant data.
import os, re
import pandas as pd

RAW = "decomposition/data/raw/gtex/"
TS = "decomposition/data/interim/track_samples_gtex.csv"
MAP_OUT = "decomposition/data/frozen/tissue_map.csv"
EXPR_OUT = "decomposition/data/interim/expr_groups_{rel}.csv.gz"
FILES = {"v8": RAW + "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
         "v10": RAW + "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"}
K562 = "Cells - Leukemia cell line (CML)"

def v10_name(s):
    return re.sub(r"\s+", "_", s.replace(" - ", "_").replace("(", "").replace(")", "").strip())

# --- 1. The sub-tissue map ---
ts = pd.read_csv(TS)
tmap = (ts.groupby(["group", "SMTSD_v8"]).size().rename("n_tracks").reset_index()
          .rename(columns={"SMTSD_v8": "subtissue_v8"}))
tmap["weight"] = tmap["n_tracks"] / tmap.groupby("group")["n_tracks"].transform("sum")
tmap["subtissue_v10"] = tmap["subtissue_v8"].map(v10_name)
tmap["source"] = ["ENCODE (step 4e)" if s == K562 else "GTEx" for s in tmap["subtissue_v8"]]
assert tmap["group"].nunique() == 31, "expected 31 groups"
assert (tmap.groupby("group")["weight"].sum() - 1).abs().max() < 1e-12, "weights must add up to 1"
print(f"sub-tissue map: {len(tmap)} rows for 31 groups")
print(tmap.groupby("group").apply(lambda d: "; ".join(f"{s} x {w:.2f}" for s, w in
      zip(d["subtissue_v8"], d["weight"])), include_groups=False).to_string())

# --- 2. Staged expression tables (plain weighted averages of median TPM) ---
gt = tmap[tmap["source"] == "GTEx"]
tables = {}
for rel, path in FILES.items():
    m = pd.read_csv(path, sep="\t", skiprows=2)
    col = "subtissue_v8" if rel == "v8" else "subtissue_v10"
    missing = sorted(set(gt[col]) - set(m.columns))
    assert not missing, f"{rel}: sub-tissues not found: {missing}"
    out = m[["Name", "Description"]].rename(columns={"Name": "gene_id", "Description": "gene_name"})
    out.insert(1, "gene_id_nov", out["gene_id"].str.replace(r"\.\d+$", "", regex=True))
    for g, sub in gt.groupby("group"):
        out[g] = sum(w * m[s] for s, w in zip(sub[col], sub["weight"]))
    out.to_csv(EXPR_OUT.format(rel=rel), index=False)
    tables[rel] = out
    print(f"\n{rel}: {len(out):,} genes x {out.shape[1] - 3} groups -> {EXPR_OUT.format(rel=rel)}")

# --- Check: do v8 and v10 agree? ---
a, b = tables["v8"], tables["v10"]
both = a.merge(b, on="gene_id_nov", suffixes=("_v8", "_v10"))
groups = sorted(gt["group"].unique())
rho = {g: both[f"{g}_v8"].corr(both[f"{g}_v10"], method="spearman") for g in groups}
print(f"\ngenes in both releases: {len(both):,}")
print(f"v8 vs v10, rank agreement across genes per group: min {min(rho.values()):.3f} "
      f"({min(rho, key=rho.get)}), median {pd.Series(rho).median():.3f}")

# --- Freeze the map ---
text = tmap.to_csv(index=False)
if os.path.exists(MAP_OUT):
    assert open(MAP_OUT).read() == text, f"{MAP_OUT} exists and differs; frozen files are never overwritten"
    print(f"\n{MAP_OUT} already exists and matches")
else:
    with open(MAP_OUT, "w") as f:
        f.write(text)
    print(f"\nsaved {MAP_OUT}")