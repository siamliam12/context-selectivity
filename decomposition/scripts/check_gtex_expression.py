# C3 step 4b: checks the downloaded GTEx files: genes and tissues in each median-expression table,
# whether every sub-tissue behind our tracks is there, and the gene list. No variant data.
import pandas as pd

RAW = "decomposition/data/raw/gtex/"
TS = pd.read_csv("decomposition/data/interim/track_samples_gtex.csv")
needed = sorted(set(TS["SMTSD_v8"]))
print(f"sub-tissues behind our 89 tracks: {len(needed)}")

for rel, fn in [("v8", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"),
                ("v10", "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")]:
    m = pd.read_csv(RAW + fn, sep="\t", skiprows=2)
    tissues = [c for c in m.columns if c not in ("Name", "Description")]
    missing = [t for t in needed if t not in tissues]
    print(f"\n{rel}: {len(m):,} genes, {len(tissues)} tissues")
    print(f"  our sub-tissues missing from this table: {missing}")
    print(f"  first gene IDs: {list(m['Name'].head(3))}")

gtf = pd.read_csv(RAW + "gencode.v26.GRCh38.genes.gtf", sep="\t", comment="#", header=None,
                  usecols=[0, 2, 3, 4, 6, 8], names=["chrom", "feature", "start", "end", "strand", "attr"])
genes = gtf[gtf["feature"] == "gene"]
print(f"\ngene list: {len(genes):,} genes on {genes['chrom'].nunique()} sequences, "
      f"e.g. {sorted(genes['chrom'].unique())[:4]}")