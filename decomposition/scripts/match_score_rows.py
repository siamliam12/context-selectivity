# C2 step 4b: finds each of our variants in the score file and compares positions and letters.
# Reads the variant list and, for matched rows only, the REF and ALT letters. Never reads scores.
import gcsfs, h5py
import numpy as np
import pandas as pd

V = "decomposition/data/interim/variant_context_hg38.csv"
IDX = "decomposition/data/interim/score_index.parquet"
SRC = "seqnn-share/sniff/borzoi_veps/scores_indel_snp_distmatched.h5"
OUT = "decomposition/data/interim/score_rows.csv"

v = pd.read_csv(V)
idx = pd.read_parquet(IDX)
chr_type = idx["chr"].dtype
print(f"our variants: {len(v):,}; score file rows: {len(idx):,}")
print("most common IDs in the score file:")
print(idx["snp"].value_counts().head(3).to_string())

# --- 1. Which genome version does the score file use? ---
q19 = pd.DataFrame({"variant": v["variant"], "chr": ("chr" + v["CHR"].astype(str)).astype(chr_type),
                    "pos": v["BP_hg19"].astype("int64")})
has38 = v["BP_hg38"].notna() & v["chr_hg38"].isin(chr_type.categories)
q38 = pd.DataFrame({"variant": v.loc[has38, "variant"], "chr": v.loc[has38, "chr_hg38"].astype(chr_type),
                    "pos": v.loc[has38, "BP_hg38"].astype("int64")})
m19 = q19.merge(idx, on=["chr", "pos"])
m38 = q38.merge(idx, on=["chr", "pos"])
print(f"\nposition found using hg19: {m19['variant'].nunique():,} of {len(q19):,}")
print(f"position found using hg38: {m38['variant'].nunique():,} of {len(q38):,}")
build = "hg38" if m38["variant"].nunique() >= m19["variant"].nunique() else "hg19"
cand = (m38 if build == "hg38" else m19).copy()
print(f"score file positions are: {build}")

# --- 2. Read REF and ALT for the candidate rows only ---
fs = gcsfs.GCSFileSystem(token="anon")

def read_strings(h, name, rows):
    ds = h[name]
    base, size = ds.id.get_offset(), ds.dtype.itemsize   # data sits in one plain block
    starts = [int(base + int(r) * size) for r in rows]
    out = []
    for i in range(0, len(starts), 1000):
        s = starts[i:i + 1000]
        got = fs.cat_ranges([SRC] * len(s), s, [x + size for x in s])
        out += [b.rstrip(b"\x00").decode("ascii") for b in got]
    return out

with fs.open(SRC, "rb") as f, h5py.File(f, "r") as h:
    # self-test: the byte-position method must give the same first 5 letters as h5py itself
    assert read_strings(h, "ref_allele", range(5)) == [x.decode() for x in h["ref_allele"][:5]]
    cand["s_ref"] = read_strings(h, "ref_allele", cand["row"])
    cand["s_alt"] = read_strings(h, "alt_allele", cand["row"])
print("letters read for", f"{len(cand):,}", "candidate rows (self-test passed)")

# --- 3. Compare letters ---
cand = cand.merge(v[["variant", "SNP", "A1", "A2", "REF", "ALT", "ref_check", "is_indel"]], on="variant")
if build == "hg19":   # the file uses hg19 positions and hg19-side letters: compare with A1/A2
    ours_ref, ours_alt = cand["A1"], cand["A2"]
else:
    ours_ref, ours_alt = cand["REF"].fillna(cand["A1"]), cand["ALT"].fillna(cand["A2"])
cand["letters"] = np.select(
    [(cand["s_ref"] == ours_ref) & (cand["s_alt"] == ours_alt),
     (cand["s_ref"] == ours_alt) & (cand["s_alt"] == ours_ref)],
    ["same", "swapped"], default="other")
cand["rank"] = cand["letters"].map({"same": 0, "swapped": 1, "other": 2})
best = cand.sort_values(["variant", "rank"]).drop_duplicates("variant")
print(f"\nvariants with more than one row at their position: "
      f"{int((cand.groupby('variant').size() > 1).sum()):,}")
print("\nletters, by our REF check:")
print(pd.crosstab(best["ref_check"].fillna("none"), best["letters"], margins=True).to_string())
print(f"\nrsID agrees (where letters match): {(best.loc[best['rank'] < 2, 'SNP'] == best.loc[best['rank'] < 2, 'snp']).mean():.4f}")

# --- 4. Variants not found at their position: try their rsID ---
found = set(best.loc[best["rank"] < 2, "variant"])
missing = v[~v["variant"].isin(found)]
print(f"\nnot found at their position with matching letters: {len(missing)}")
by_id = missing[["variant", "SNP", "A1", "A2", "liftover_status", "ref_check"]].merge(
    idx, left_on="SNP", right_on="snp")
print(f"of these, rsID present in the score file: {by_id['variant'].nunique()}")
print(by_id.head(40).to_string(index=False))

best = best[best["rank"] < 2][["variant", "row", "chr", "pos", "s_ref", "s_alt", "letters"]]
best = best.rename(columns={"row": "score_row", "chr": "score_chr", "pos": "score_pos"})
best.to_csv(OUT, index=False)
print(f"\nsaved {OUT}: {len(best):,} variants with a score-file row")