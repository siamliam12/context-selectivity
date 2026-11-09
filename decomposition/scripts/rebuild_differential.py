# Rebuilds Sniff-promoted and Sniff-demoted variants (cut-off 0.25) from the fine-mapping files.
# Checks them against the paper's counts and against S3. Reads PIPs only, no Borzoi scores.
import os, glob
import pandas as pd

# --- Settings ---
RAW = "decomposition/data/raw/sniff/finemap"
XLSX = "decomposition/data/raw/sniff/media-1.xlsx"
OUT = "decomposition/data/interim/differential_variants.csv"
CUT = 0.25
KEY = ["CHR", "BP", "A1", "A2"]
COLS = ["CHR", "SNP", "BP", "A1", "A2", "MAF", "N", "PIP", "CREDIBLE_SET"]

def vkey(df):
    return (df["CHR"].astype(int).astype(str) + ":" + df["BP"].astype(int).astype(str)
            + ":" + df["A1"].astype(str) + ":" + df["A2"].astype(str))

# --- Find the traits ---
traits = sorted({os.path.basename(p).split(".")[0] for p in glob.glob(f"{RAW}/sniff/*.txt.gz")})
print(f"{len(traits)} traits found")

# --- Main loop: one trait, one chromosome at a time ---
parts = []
for trait in traits:
    n_total = 0
    for chrom in range(1, 23):
        b = pd.read_csv(f"{RAW}/baseline/{trait}.{chrom}.txt.gz", sep=r"\s+", usecols=COLS)
        s = pd.read_csv(f"{RAW}/sniff/{trait}.{chrom}.txt.gz", sep=r"\s+", usecols=COLS)
        assert len(b) == len(s), f"{trait} chr{chrom}: row counts differ"
        assert not b.duplicated(KEY).any() and not s.duplicated(KEY).any(), f"{trait} chr{chrom}: duplicates"
        m = b.merge(s, on=KEY, suffixes=("_base", "_sniff"), validate="one_to_one")
        assert len(m) == len(b), f"{trait} chr{chrom}: {len(b) - len(m)} variants did not match"
        assert (m["SNP_base"] == m["SNP_sniff"]).all(), f"{trait} chr{chrom}: SNP IDs differ"
        n_total += len(m)
        m["dPIP"] = m["PIP_sniff"] - m["PIP_base"]
        keep = m[m["dPIP"].abs() > CUT].copy()
        keep["direction"] = ["promoted" if x > 0 else "demoted" for x in keep["dPIP"]]
        keep["trait"] = trait
        parts.append(keep)
    print(f"{trait}: {n_total:,} variants, N = {int(b['N'].iloc[0]):,}", flush=True)

# --- What does credible-set index 0 mean? ---
share0 = s["CREDIBLE_SET"].astype(str).str.endswith(":0").mean()
print(f"\nShare of rows labelled ':0' in the last file read: {share0:.4f}")

# --- Combine and count ---
d = pd.concat(parts, ignore_index=True)
d["variant"] = vkey(d)
print(f"\n== counts at cut-off {CUT} ==")
for direction, expected in [("promoted", 3160), ("demoted", 2579)]:
    sub = d[d["direction"] == direction]
    print(f"{direction}: {len(sub):,} rows, {sub['variant'].nunique():,} unique (paper: {expected:,})")
print("\nRows per trait:")
print(d.groupby(["trait", "direction"]).size().unstack(fill_value=0).to_string())

# --- Compare the promoted list with S3 ---
s3 = pd.read_excel(XLSX, sheet_name="S3")
s3["variant"] = vkey(s3)
prom = d[d["direction"] == "promoted"]
c = prom.merge(s3, on=["trait", "variant"], how="outer", indicator=True)
print("\n== promoted vs S3 (variant-trait rows) ==")
print(c["_merge"].value_counts().to_string())
both = c[c["_merge"] == "both"]
print("max |PIP_base difference|: ", (both["PIP_base_x"] - both["PIP_base_y"]).abs().max())
print("max |Sniff PIP difference|:", (both["PIP_sniff"] - both["PIP_borzoi"]).abs().max())
print("credible-set labels equal: ", (both["CREDIBLE_SET_sniff"].astype(str) == both["CS_borzoi"].astype(str)).mean())

# --- Extras ---
print("\nS3 rows with cred_id 0:", int((pd.to_numeric(s3["cred_id"]) == 0).sum()))
print("S3 credible sets per trait (trait + CS_borzoi):")
print(s3.groupby("trait")["CS_borzoi"].nunique().to_string())
s2 = pd.read_excel(XLSX, sheet_name="S2")
print("\nS2 sample size per figure label:")
print(s2.groupby("Trait")["N"].unique().to_string())

# --- Save ---
os.makedirs(os.path.dirname(OUT), exist_ok=True)
d.to_csv(OUT, index=False)
print(f"\nsaved {OUT}: {len(d):,} rows")
