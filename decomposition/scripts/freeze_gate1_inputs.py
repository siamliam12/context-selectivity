# C5: freezes two inputs for Gate 1, before any table is built. No scores involved.
#   fig4c_labels.csv    the 17 labelled (coloured) dots of Fig. 4C, read by eye on 16 Nov 2026
#   set_c_variants.csv  set C, the "Sniff-prioritized" variants (Methods definition), built from S3
import os
import pandas as pd

FROZEN = "decomposition/data/frozen/"
S3 = "decomposition/data/raw/sniff/media-1.xlsx"

def write_once(df, path):
    text = df.to_csv(index=False)
    if os.path.exists(path):
        assert open(path).read() == text, f"{path} exists and differs; frozen files are never overwritten"
        print(f"{path} already exists and matches")
    else:
        with open(path, "w") as f:
            f.write(text)
        print(f"saved {path}: {len(df)} rows")

tm = pd.read_csv(FROZEN + "trait_map.csv")
canon = {str(l).lower(): l for l in tm["label"]}
groups = set(pd.read_csv(FROZEN + "track_map.csv")["group"])

# --- 1. Fig. 4C readings (coloured dots = tissues the authors labelled) ---
lab = pd.DataFrame([
    ("HLSRC", "bone_marrow", 0.82, False),
    ("MPV", "blood", 0.89, False), ("MPV", "bone_marrow", 0.78, False), ("MPV", "spleen", 0.50, True),
    ("Cholesterol", "liver", 2.42, False),
    ("MCH", "bone_marrow", 1.08, False),
    ("MC", "blood", 1.48, False), ("MC", "spleen", 0.83, False), ("MC", "lung", 0.50, True),
    ("Balding", "blood", 0.71, False), ("Balding", "muscle", 0.63, False),
    ("Chronotype", "brain", 0.79, False),
    ("Hair pigment", "blood", 0.72, False),
    ("Eczema", "blood", 1.29, False), ("Eczema", "spleen", 0.60, False),
    ("Age at Menarche", "brain", 0.70, False),
    ("Diastolic BP", "heart", 0.59, False),
], columns=["label", "group", "approx_z", "on_line"])
lab["label"] = lab["label"].str.lower().map(canon)
assert lab["label"].notna().all(), "a trait label does not match trait_map.csv"
assert set(lab["group"]) <= groups, f"unknown group: {set(lab['group']) - groups}"
print(f"{len(lab)} labelled dots in {lab['label'].nunique()} traits; "
      f"traits with none: {sorted(set(tm['label']) - set(lab['label']))}")
write_once(lab, FROZEN + "fig4c_labels.csv")

# --- 2. Set C from S3: PIP < 0.1 in SuSiE and baseline, Sniff raises it by >= 0.25 ---
s3 = pd.read_excel(S3, sheet_name="S3")
for k in ["A1", "A2"]:
    s3[k] = s3[k].astype(str)
s3["variant"] = s3["CHR"].astype(str) + ":" + s3["BP"].astype(str) + ":" + s3["A1"] + ":" + s3["A2"]
s3["dpip"] = s3["PIP_borzoi"] - s3["PIP_base"]
c = s3[(s3["PIP_susie"] < 0.1) & (s3["PIP_base"] < 0.1) & (s3["dpip"] >= 0.25)].copy()
c["is_indel"] = (c["A1"].str.len() > 1) | (c["A2"].str.len() > 1)
print(f"\nset C: {len(c)} variant-trait rows, {c['variant'].nunique()} unique variants, "
      f"{c.loc[~c['is_indel'], 'variant'].nunique()} unique SNPs (the paper reports 1,158)")

match = [k for k in tm.columns if set(c["trait"].astype(str)) <= set(tm[k].astype(str))]
assert match, f"S3 trait names not found in trait_map.csv: {sorted(set(c['trait']))}"
c["label"] = c["trait"].astype(str).map(dict(zip(tm[match[0]].astype(str), tm["label"])))
known = set(pd.read_csv(FROZEN + "variants.csv")["variant"])
assert set(c["variant"]) <= known, "set C has variants missing from variants.csv"
print("rows per trait:", c["label"].value_counts().to_dict())

out = c[["variant", "label", "PIP_susie", "PIP_base", "PIP_borzoi", "dpip", "is_indel", "CS_borzoi"]]
out = out.sort_values(["label", "variant"]).reset_index(drop=True)
write_once(out, FROZEN + "set_c_variants.csv")