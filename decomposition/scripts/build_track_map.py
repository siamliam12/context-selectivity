# C3 step 2: builds the frozen map for the 89 tracks behind Fig. 4C: score-file column, tissue group,
# and the model's own output index from targets_human.txt (needed in Part E). Reads names only. No scores.
import os
import pandas as pd

TRACKS = "decomposition/data/interim/score_tracks.csv"
TARGETS = "/home/fsh33/research/main/decomposition/data/raw/borzoi/targets_human.txt"
OUT = "decomposition/data/frozen/track_map.csv"

t = pd.read_csv(TRACKS)
sel = t[t["target_id"].str.startswith("GTEX") | t["target_id"].str.startswith("K-562")].copy()
sel["group"] = sel["target_label"].str.replace("RNA:", "", regex=False)
sel["source"] = sel["target_id"].str.startswith("GTEX").map({True: "GTEx", False: "GTEx K562"})

sizes = sel.groupby("group").size()
print(f"tracks: {len(sel)}, groups: {sel['group'].nunique()}")
print("groups by number of tracks:", sizes.value_counts().sort_index().to_dict())
print("groups with 2 tracks:", sorted(sizes[sizes == 2].index))
assert len(sel) == 89 and sel["group"].nunique() == 31, "expected 89 tracks in 31 groups"

# --- Link each track to the model's output index, by name ---
tg = pd.read_csv(TARGETS, sep="\t", index_col=0)
print(f"\ntargets_human.txt: {len(tg):,} rows; columns: {list(tg.columns)}")
strip = lambda s: s.str.replace(r"[+-]$", "", regex=True)
tg["id_core"] = strip(tg["identifier"].astype(str))
sel["id_core"] = strip(sel["target_id"])
links = sel[["col", "id_core"]].merge(tg.reset_index().rename(columns={tg.index.name or "index": "model_index"}),
                                      on="id_core", how="left")
per = links.groupby("col").agg(model_index=("model_index", lambda x: ";".join(str(int(i)) for i in x.dropna())),
                               n_model_tracks=("model_index", "count"),
                               description=("description", "first"))
sel = sel.merge(per, left_on="col", right_index=True)
print("tracks by number of model outputs found:", sel["n_model_tracks"].value_counts().to_dict())
print(sel[["col", "target_id", "group", "model_index", "description"]].head(6).to_string(index=False))
print(sel[sel["source"] == "GTEx K562"][["col", "target_id", "group", "model_index", "description"]].to_string(index=False))
assert (sel["n_model_tracks"] > 0).all(), "every track must link to at least one model output"

out = sel[["col", "target_id", "target_label", "group", "source", "model_index", "n_model_tracks", "description"]]
out = out.sort_values(["group", "col"]).reset_index(drop=True)
print("\nthe 31 groups:", ", ".join(sorted(out["group"].unique())))

text = out.to_csv(index=False)
if os.path.exists(OUT):
    assert open(OUT).read() == text, f"{OUT} exists and differs; frozen files are never overwritten"
    print(f"{OUT} already exists and matches")
else:
    with open(OUT, "w") as f:
        f.write(text)
    print(f"saved {OUT}: {len(out)} tracks")