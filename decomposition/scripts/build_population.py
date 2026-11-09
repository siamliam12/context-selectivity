# Builds the tidy C1 variant table and the trait map.
# Input: the rebuilt list from rebuild_differential.py. Reads no Borzoi scores.
import os
import pandas as pd

IN = "decomposition/data/interim/differential_variants.csv"
XLSX = "decomposition/data/raw/sniff/media-1.xlsx"
OUT = "decomposition/data/interim/variant_population.csv"
MAP_OUT = "decomposition/data/frozen/trait_map.csv"

# File code -> figure label, matched by name (S2's labels are shifted, so S2 is not used).
CODE_TO_LABEL = {
    "biochemistry_Cholesterol": "Cholesterol",
    "blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT": "HLSRC",
    "blood_MEAN_CORPUSCULAR_HEMOGLOBIN": "MCH",
    "blood_MEAN_PLATELET_VOL": "MPV",
    "blood_MONOCYTE_COUNT": "MC",
    "bmd_HEEL_TSCOREz": "BMD",
    "body_BALDING1": "Balding",
    "body_HEIGHTz": "Height",
    "bp_DIASTOLICadjMEDz": "Diastolic BP",
    "disease_ALLERGY_ECZEMA_DIAGNOSED": "Eczema",
    "lung_FEV1FVCzSMOKE": "FEV:FVC",
    "lung_FVCzSMOKE": "FVC",
    "other_MORNINGPERSON": "Chronotype",
    "pigment_HAIR": "Hair pigment",
    "repro_MENARCHE_AGE": "Age at Menarche",
}

d = pd.read_csv(IN)
print(f"read {len(d):,} rows")

# --- Trait map, checked against S1 ---
s1 = pd.read_excel(XLSX, sheet_name="S1")
tm = pd.DataFrame(list(CODE_TO_LABEL.items()), columns=["trait", "label"])
tm = tm.merge(s1, left_on="label", right_on="Abbreviation", how="left", validate="one_to_one")
assert tm["Trait"].notna().all(), "a label is missing from S1"
assert set(tm["trait"]) == set(d["trait"]), "trait codes in the data and the map differ"
tm = tm.rename(columns={"Trait": "full_name"})[["trait", "label", "full_name"]]
print(tm.to_string(index=False))

# --- Tidy table ---
t = pd.DataFrame({
    "trait": d["trait"],
    "label": d["trait"].map(CODE_TO_LABEL),
    "variant": d["variant"],
    "SNP": d["SNP_base"],
    "CHR": d["CHR"],
    "BP_hg19": d["BP"],
    "A1": d["A1"],
    "A2": d["A2"],
    "MAF": d["MAF_base"],
    "PIP_base": d["PIP_base"],
    "PIP_sniff": d["PIP_sniff"],
    "dPIP": d["dPIP"].round(10),
    "direction": d["direction"],
    "sets": d["direction"].map({"promoted": "A+B", "demoted": "B"}),
    "cs_sniff": d["CREDIBLE_SET_sniff"],
    "cs_base": d["CREDIBLE_SET_base"],
})
t["in_cs_sniff"] = ~t["cs_sniff"].astype(str).str.endswith(":0")
t["in_cs_base"] = ~t["cs_base"].astype(str).str.endswith(":0")
t["is_indel"] = (t["A1"].str.len() > 1) | (t["A2"].str.len() > 1)

# --- Checks ---
assert t["label"].notna().all(), "a trait has no label"
assert not t.duplicated(["trait", "variant"]).any(), "a variant appears twice in one trait"

print("\nRows per trait and direction:")
print(t.groupby(["label", "direction"]).size().unstack(fill_value=0).to_string())

for name, sub in [("Set A (promoted only)", t[t["sets"] == "A+B"]), ("Set B (promoted + demoted)", t)]:
    print(f"\n{name}: {len(sub):,} rows, {sub['variant'].nunique():,} unique variants, "
          f"{int(sub['is_indel'].sum()):,} indel rows, "
          f"{sub['in_cs_sniff'].mean():.3f} in a Sniff credible set")

dem = t[t["direction"] == "demoted"]
print(f"\nDemoted rows in a baseline credible set: {dem['in_cs_base'].mean():.3f}")
both = set(t.loc[t["direction"] == "promoted", "variant"]) & set(dem["variant"])
print(f"Variants promoted in one trait and demoted in another: {len(both)}")

# --- Save ---
t.to_csv(OUT, index=False)
print(f"\nsaved {OUT}: {len(t):,} rows")
if os.path.exists(MAP_OUT):
    old = pd.read_csv(MAP_OUT)
    assert old.equals(tm), f"{MAP_OUT} exists and differs; frozen files are never overwritten"
    print(f"{MAP_OUT} already exists and matches")
else:
    os.makedirs(os.path.dirname(MAP_OUT), exist_ok=True)
    tm.to_csv(MAP_OUT, index=False)
    print(f"saved {MAP_OUT}")