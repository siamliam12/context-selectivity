# C3 step 4g: builds the K562 (bone_marrow) expression column from the 16 ENCODE files
# (8 experiments x 2 replicates), checks it, and adds it to the staged GTEx tables.
# Gene expression only. No variant data.
import os
import numpy as np
import pandas as pd

FILES = "decomposition/data/interim/encode_k562_files.csv"
RAW = "decomposition/data/raw/encode/"
STAGED = "decomposition/data/interim/expr_groups_{rel}.csv.gz"
OUT = "decomposition/data/interim/expr_groups_{rel}_31.csv.gz"
K562_OUT = "decomposition/data/interim/k562_encode_tpm.csv.gz"

# --- Read the 16 files ---
ch = pd.read_csv(FILES)
cols = {}
for r in ch.itertuples():
    t = pd.read_csv(os.path.join(RAW, f"{r.experiment}_{r.file}.tsv"), sep="\t", usecols=["gene_id", "TPM"])
    t = t[t["gene_id"].str.startswith("ENSG")]                    # drop spike-in controls
    t["gene_id_nov"] = t["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    s = t.groupby("gene_id_nov")["TPM"].sum()
    cols[(r.experiment, r.bio_reps)] = s / s.sum() * 1e6         # re-scale to 1 million over genes only
mat = pd.DataFrame(cols)
print(f"{mat.shape[1]} files; {len(mat):,} genes")

# --- Replicates within each experiment, then experiments with each other ---
exps = sorted(set(e for e, _ in mat.columns))
rep_rho = {e: mat[e].corr(method="spearman").iloc[0, 1] for e in exps}
print("replicate agreement within each experiment (rank correlation):")
print("  " + ", ".join(f"{e} {v:.3f}" for e, v in rep_rho.items()))
expm = pd.DataFrame({e: mat[e].mean(axis=1) for e in exps})
ex_rho = expm.corr(method="spearman")
off = ex_rho.where(~np.eye(len(exps), dtype=bool)).stack()
print(f"agreement between experiments: min {off.min():.3f}, median {off.median():.3f}")

k562 = expm.mean(axis=1).rename("bone_marrow")
k562.to_frame().to_csv(K562_OUT)
print(f"saved {K562_OUT}")

# --- Does it look like K562? ---
names = pd.read_csv(STAGED.format(rel="v8"), usecols=["gene_id_nov", "gene_name"]).set_index("gene_id_nov")["gene_name"]
top = k562.sort_values(ascending=False).head(10)
print("\ntop 10 K562 genes (TPM):", ", ".join(f"{names.get(g, g)} {v:.0f}" for g, v in top.items()))
for gname in ["HBG1", "HBG2", "GATA1", "TAL1", "ALB", "PTPRC"]:
    vals = [round(k562.get(i, np.nan), 1) for i in names[names == gname].index]
    print(f"  {gname}: {vals}")

# --- Add it to the staged tables ---
for rel in ["v8", "v10"]:
    st = pd.read_csv(STAGED.format(rel=rel))
    st["bone_marrow"] = st["gene_id_nov"].map(k562)
    ok = st["bone_marrow"].notna()
    gtex_groups = [c for c in st.columns if c not in ("gene_id", "gene_id_nov", "gene_name", "bone_marrow")]
    rho = {g: st.loc[ok, g].corr(st.loc[ok, "bone_marrow"], method="spearman") for g in gtex_groups}
    best = sorted(rho.items(), key=lambda kv: -kv[1])[:3]
    print(f"\n{rel}: K562 value found for {ok.mean():.1%} of genes; most similar GTEx groups: "
          + ", ".join(f"{g} {v:.3f}" for g, v in best))
    st.to_csv(OUT.format(rel=rel), index=False)
    print(f"saved {OUT.format(rel=rel)}: {st.shape[1] - 3} groups")