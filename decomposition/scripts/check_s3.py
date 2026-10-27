import pandas as pd

d = pd.read_excel("data/raw/sniff/media-1.xlsx", "S3")

print("rows:", len(d))
print("unique SNPs:", d["SNP"].nunique())
print("traits:", d["trait"].nunique())
print("missing values:", int(d.isna().sum().sum()))
print("all rows have PIP_borzoi - PIP_base > 0.25:", bool((d["PIP_borzoi"] - d["PIP_base"] > 0.25).all()))
print("all A1 and A2 single letters:", bool(d["A1"].str.len().eq(1).all() and d["A2"].str.len().eq(1).all()))

f = d[(d["PIP_susie"] < 0.1) & (d["PIP_base"] < 0.1)]
print("unique SNPs after filter (PIP_susie < 0.1 and PIP_base < 0.1):", f["SNP"].nunique())

print()
print("rows per trait:")
print(d["trait"].value_counts().to_string())
