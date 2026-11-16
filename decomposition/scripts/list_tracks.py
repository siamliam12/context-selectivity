# C3 step 1: lists the 6,385 track names in both score files and finds the GTEx and K562 RNA tracks.
# Reads track names only. No scores.
import gcsfs, h5py
import pandas as pd

fs = gcsfs.GCSFileSystem(token="anon")
BASE = "seqnn-share/sniff/borzoi_veps/"
OUT = "decomposition/data/interim/score_tracks.csv"

def tracks(name):
    with fs.open(BASE + name, "rb") as f, h5py.File(f, "r") as h:
        ids = [x.decode() for x in h["target_ids"][:]]
        labels = [x.decode() for x in h["target_labels"][:]]
    return pd.DataFrame({"col": range(len(ids)), "target_id": ids, "target_label": labels})

dm = tracks("scores_indel_snp_distmatched.h5")
raw = tracks("scores.h5")
print(f"distmatched: {len(dm):,} tracks; raw: {len(raw):,} tracks")
print("same ids and labels, same order, in both files:", dm.equals(raw))

dm["kind"] = dm["target_label"].str.split(":").str[0]
print("\ntracks by kind (text before the first colon in the label):")
print(dm["kind"].value_counts().to_string())

gtex = dm[dm["target_id"].str.startswith("GTEX")]
print(f"\nGTEx tracks (id starts with GTEX): {len(gtex)}")
print("labels and how many tracks each:")
print(gtex["target_label"].value_counts().sort_index().to_string())

rna = dm[dm["kind"] == "RNA"]
k562 = rna[rna["target_id"].str.upper().str.startswith("K-562")
           | rna["target_label"].str.contains("K562|K-562|leukemia|CML", case=False)]
print(f"\nK562-like RNA tracks: {len(k562)}")
print(k562.to_string(index=False))

dm.to_csv(OUT, index=False)
print(f"\nsaved {OUT}")