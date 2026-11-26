# D5 step 2 (also E1's input): for every variant, the genes in the score's window (+-98,304 bp, hg38) and
# their summed expression per tissue group. Also the nearest gene, and whether the variant sits inside a gene.
import re
import numpy as np, pandas as pd

HALF = 98_304
GTF = "decomposition/data/raw/gtex/gencode.v26.GRCh38.genes.gtf"
EXPR = {"v8": "decomposition/data/interim/expr_groups_v8_31.csv.gz",
        "v10": "decomposition/data/interim/expr_groups_v10_31.csv.gz"}
OUT = "decomposition/data/interim/locus_expression.parquet"

# Genes: position, strand, start site
g = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                names=["chrom", "feature", "start", "end", "strand", "attr"])
g = g[g["feature"] == "gene"].copy()
g["gene_id"] = g["attr"].str.extract(r'gene_id "([^"]+)"')[0]
g["gene_id_nov"] = g["gene_id"].str.replace(r"\.\d+$", "", regex=True)
g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"])
g = g.reset_index(drop=True)

# Variants: the D5 sample and the Fig. 4C variants, hg38 positions
d5 = pd.read_csv("decomposition/data/interim/d5_sample_variants.csv")
fv = pd.read_csv("decomposition/data/frozen/variants.csv")
loci = pd.concat([d5[["variant", "chr_hg38", "BP_hg38"]].assign(source="d5"),
                  fv[["variant", "chr_hg38", "BP_hg38"]].assign(source="fig")])
print("without an hg38 position (dropped):", loci.groupby("source")["BP_hg38"].apply(lambda s: int(s.isna().sum())).to_dict())
loci = loci.dropna(subset=["BP_hg38"]).reset_index(drop=True)
loci["BP_hg38"] = loci["BP_hg38"].astype(int)

# For each variant: window genes, nearest gene (by start site), inside a gene body?
win, near, inside = [None] * len(loci), [None] * len(loci), np.zeros(len(loci), bool)
for c, lc in loci.groupby("chr_hg38"):
    gc = g[g["chrom"] == c]
    if gc.empty:
        continue
    st, en, ix = gc["start"].to_numpy(), gc["end"].to_numpy(), gc.index.to_numpy()
    order = np.argsort(gc["tss"].to_numpy())
    tss_sorted, ix_by_tss = gc["tss"].to_numpy()[order], ix[order]
    for i, p in zip(lc.index, lc["BP_hg38"].to_numpy()):
        hit = (st <= p + HALF) & (en >= p - HALF)
        win[i] = ix[hit]
        inside[i] = bool(((st <= p) & (en >= p)).any())
        k = np.searchsorted(tss_sorted, p)
        cand = [j for j in (k - 1, k) if 0 <= j < len(tss_sorted)]
        near[i] = ix_by_tss[min(cand, key=lambda j: abs(tss_sorted[j] - p))]
loci["n_genes"] = [len(w) if w is not None else 0 for w in win]
loci["inside_gene"] = inside
print("genes in the window: median", int(loci["n_genes"].median()), "| none:", int((loci["n_genes"] == 0).sum()))

# Summed expression per group (TPM), window and nearest-gene versions, for v8 and v10
rows = []
for ver, path in EXPR.items():
    e = pd.read_csv(path)
    groups = [c for c in e.columns if c not in ("gene_id", "gene_id_nov", "gene_name")]
    E = g[["gene_id_nov"]].merge(e.drop_duplicates("gene_id_nov"), on="gene_id_nov", how="left")[groups].fillna(0).to_numpy()
    for mode, idxs in [("window", win), ("nearest", [np.array([n]) if n is not None else np.array([], int) for n in near])]:
        S = np.vstack([E[w].sum(axis=0) if w is not None and len(w) else np.zeros(len(groups)) for w in idxs])
        df = pd.DataFrame(S, columns=groups)
        df.insert(0, "genes", mode)
        df.insert(0, "version", ver)
        rows.append(pd.concat([loci[["variant", "source", "n_genes", "inside_gene"]], df], axis=1))
out = pd.concat(rows, ignore_index=True)
out.to_parquet(OUT, index=False)
print(f"saved {OUT}: {len(out):,} rows ({len(loci):,} variants x 2 versions x 2 gene modes)")