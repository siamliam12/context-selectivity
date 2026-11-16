# C3 step 4a: finds the GTEx tissue and sub-tissue of each of the 89 track samples,
# using GTEx's public sample lists for v8 and v10. Reads sample labels only.
import pandas as pd
pd.set_option("display.max_colwidth", None)
pd.set_option("display.width", 250)

MAP = "decomposition/data/frozen/track_map.csv"
RAW = "decomposition/data/raw/gtex/"
OUT = "decomposition/data/interim/track_samples_gtex.csv"
FILES = {"v8": RAW + "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
         "v10": RAW + "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"}

tm = pd.read_csv(MAP)
tm["sampid"] = tm["target_id"].str.replace(r"\.\d+$", "", regex=True)   # drop the trailing ".1"

found = {}
for rel, path in FILES.items():
    sa = pd.read_csv(path, sep="\t", usecols=["SAMPID", "SMTS", "SMTSD", "SMAFRZE"], low_memory=False)
    m = tm.merge(sa, left_on="sampid", right_on="SAMPID", how="left")
    print(f"GTEx {rel}: {len(sa):,} samples; our 89 found: {int(m['SAMPID'].notna().sum())}; "
          f"analysis set: {m['SMAFRZE'].value_counts(dropna=False).to_dict()}")
    found[rel] = m

both = found["v8"][["group", "col", "target_id", "SMTS", "SMTSD", "SMAFRZE"]].merge(
    found["v10"][["col", "SMTSD", "SMAFRZE"]], on="col", suffixes=("_v8", "_v10"))
print("\nsub-tissue identical in v8 and v10 for all 89:", bool((both["SMTSD_v8"] == both["SMTSD_v10"]).all()))

print("\nsub-tissues per group (number of track samples):")
print(both.groupby("group")["SMTSD_v8"].apply(
    lambda s: "; ".join(f"{k} ({n})" for k, n in s.value_counts().sort_index().items())).to_string())

print("\ntrack samples outside GTEx's RNA-seq analysis set:")
out_set = both[(both["SMAFRZE_v8"] != "RNASEQ") | (both["SMAFRZE_v10"] != "RNASEQ")]
print(out_set[["group", "target_id", "SMTSD_v8", "SMAFRZE_v8", "SMAFRZE_v10"]].to_string(index=False))

both.to_csv(OUT, index=False)
print(f"\nsaved {OUT}")