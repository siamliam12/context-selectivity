# C4 (option B): freezes the authors' own trait-tissue claims from the Sniff preprint as the relevance
# standard. These are the six printed numbers also used at Gate 1. No results involved.
import os
import pandas as pd

OUT = "decomposition/data/frozen/author_pairs.csv"
pairs = pd.DataFrame([
    ("Cholesterol",  "liver",       2.42, "cholesterol: most likely liver"),
    ("MCH",          "bone_marrow", 1.08, "hemoglobin: most likely bone marrow"),
    ("HLSRC",        "bone_marrow", 0.78, "reticulocyte count: most likely bone marrow"),
    ("MC",           "blood",       1.48, "monocyte count: most likely blood"),
    ("MC",           "spleen",      0.50, "monocyte count: also spleen"),
    ("Diastolic BP", "heart",       0.59, "diastolic BP: heart"),
], columns=["label", "group", "printed_z", "author_claim"])

traits = set(pd.read_csv("decomposition/data/frozen/trait_map.csv")["label"])
groups = set(pd.read_csv("decomposition/data/frozen/track_map.csv")["group"])
assert set(pairs["label"]) <= traits, f"unknown trait: {set(pairs['label']) - traits}"
assert set(pairs["group"]) <= groups, f"unknown group: {set(pairs['group']) - groups}"
print(pairs.to_string(index=False))

text = pairs.to_csv(index=False)
if os.path.exists(OUT):
    assert open(OUT).read() == text, f"{OUT} exists and differs; frozen files are never overwritten"
    print(f"{OUT} already exists and matches")
else:
    with open(OUT, "w") as f:
        f.write(text)
    print(f"saved {OUT}: {len(pairs)} pairs")