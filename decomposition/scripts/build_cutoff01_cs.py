# Gate 1 try: the cut-off 0.1 lists again, now keeping each variant's credible-set labels (Sniff and baseline).
import glob, os
import numpy as np, pandas as pd

FM = "decomposition/data/raw/sniff/finemap/"
OUT = "decomposition/data/interim/sets_cutoff01_cs.csv"
vt = pd.read_csv("decomposition/data/frozen/variant_traits.csv")
label_of = vt.drop_duplicates("trait").set_index("trait")["label"]
keep = ["CHR", "BP", "A1", "A2", "PIP", "CREDIBLE_SET"]
types = {"A1": str, "A2": str, "CREDIBLE_SET": str}

parts = []
files = sorted(glob.glob(FM + "sniff/*.txt.gz"))
for i, f in enumerate(files):
    s = pd.read_csv(f, sep="\t", usecols=keep, dtype=types)
    b = pd.read_csv(os.path.join(FM, "baseline", os.path.basename(f)), sep="\t", usecols=keep, dtype=types)
    m = s.merge(b, on=["CHR", "BP", "A1", "A2"], suffixes=("_sniff", "_base"))
    m["dPIP"] = m["PIP_sniff"] - m["PIP_base"]
    m = m[m["dPIP"].abs() > 0.1].copy()
    m["trait"] = os.path.basename(f).rsplit(".", 3)[0]
    parts.append(m)
    if i % 30 == 0:
        print(f"  {i}/{len(files)} files read", flush=True)

d = pd.concat(parts, ignore_index=True)
d["label"] = d["trait"].map(label_of)
d["variant"] = d["CHR"].astype(str) + ":" + d["BP"].astype(str) + ":" + d["A1"] + ":" + d["A2"]
d["direction"] = np.where(d["dPIP"] > 0, "promoted", "demoted")
d["is_indel"] = (d["A1"].str.len() > 1) | (d["A2"].str.len() > 1)
d = d.rename(columns={"CREDIBLE_SET_sniff": "cs_sniff", "CREDIBLE_SET_base": "cs_base"})
d["in_cs_sniff"] = ~d["cs_sniff"].str.endswith(":0")
d["in_cs_base"] = ~d["cs_base"].str.endswith(":0")

# Check: at 0.25, our credible-set flags must agree with the frozen C1 table
chk = d.merge(vt[["variant", "trait", "in_cs_sniff", "in_cs_base"]], on=["variant", "trait"], suffixes=("", "_c1"))
print(f"rows also in the C1 table: {len(chk):,} (should be 5,830)")
print(f"Sniff flag agrees: {(chk['in_cs_sniff'] == chk['in_cs_sniff_c1'].astype(bool)).mean():.4f}, "
      f"baseline flag agrees: {(chk['in_cs_base'] == chk['in_cs_base_c1'].astype(bool)).mean():.4f}")
d[["variant", "trait", "label", "dPIP", "direction", "is_indel", "cs_sniff", "cs_base",
   "in_cs_sniff", "in_cs_base"]].to_csv(OUT, index=False)
print(f"saved {OUT}: {len(d):,} rows")