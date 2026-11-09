# Checks which fine-mapping file each S2 trait label really comes from,
# by matching rows on the exact variant AND the exact Sniff PIP.
import pandas as pd

d = pd.read_csv("decomposition/data/interim/differential_variants.csv")
s2 = pd.read_excel("decomposition/data/raw/sniff/media-1.xlsx", sheet_name="S2")

def vkey(df):
    return (df["CHR"].astype(int).astype(str) + ":" + df["BP"].astype(int).astype(str)
            + ":" + df["A1"].astype(str) + ":" + df["A2"].astype(str))
s2["variant"] = vkey(s2)

m = d.merge(s2, on="variant", how="inner")
m = m[(m["PIP_sniff"] - m["PIP"]).abs() < 1e-9]
print(f"{len(m):,} rows match on variant and exact Sniff PIP\n")

print("File code -> S2 label it matches (number of rows):")
pairs = (m.groupby(["trait", "Trait"]).size().reset_index(name="rows")
           .sort_values(["trait", "rows"], ascending=[True, False]))
print(pairs.to_string(index=False))

print("\nShare of matched rows where N agrees:", (m["N_sniff"] == m["N"]).mean())

print("\nPCSK9 check, rs7552841 in S2:")
print(s2[s2["SNP"] == "rs7552841"][["SNP", "PIP", "N", "Trait"]].to_string(index=False))