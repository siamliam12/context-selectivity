# C3 step 4f: for the K562 experiments that pass our rule, lists their GRCh38 gene-expression files
# and downloads ENCODE's default ones into data/raw/encode/. Gene expression only.
import json, os, urllib.request, urllib.parse
import pandas as pd
pd.set_option("display.width", 250)

RAW = "decomposition/data/raw/encode/"
EXPS = ["ENCSR000CPH", "ENCSR000EYO", "ENCSR637VLS", "ENCSR000AEM",
        "ENCSR000AEO", "ENCSR000AEP", "ENCSR000AEQ", "ENCSR545DKY"]

params = [("type", "File"), ("output_type", "gene quantifications"), ("file_format", "tsv"),
          ("assembly", "GRCh38"), ("status", "released"), ("format", "json"), ("limit", "all")]
params += [("dataset", f"/experiments/{e}/") for e in EXPS]
for f in ["accession", "dataset", "biological_replicates", "genome_annotation",
          "file_size", "href", "preferred_default", "date_created"]:
    params.append(("field", f))
url = "https://www.encodeproject.org/search/?" + urllib.parse.urlencode(params)
req = urllib.request.Request(url, headers={"Accept": "application/json"})
graph = json.load(urllib.request.urlopen(req))["@graph"]

files = pd.DataFrame([{
    "experiment": g["dataset"].split("/")[-2], "file": g["accession"],
    "bio_reps": g.get("biological_replicates"), "annotation": g.get("genome_annotation"),
    "default": bool(g.get("preferred_default", False)),
    "MB": round(g.get("file_size", 0) / 1e6, 1), "created": g.get("date_created", "")[:10],
    "href": g.get("href")} for g in graph])
files = files.sort_values(["experiment", "default", "annotation"], ascending=[True, False, False])
print(f"{len(files)} GRCh38 gene-quantification files across "
      f"{files['experiment'].nunique()} of {len(EXPS)} experiments")
print(files.drop(columns="href").to_string(index=False))

# The rule (revised 16 Nov, before any expression value was read): GENCODE V29 files, which every
# experiment has, with all biological replicates. ENCODE's "default" flag marks only one replicate.
chosen = files[files["annotation"] == "V29"].copy()
reps = chosen.groupby("experiment")["bio_reps"].apply(lambda s: sorted(r[0] for r in s))
print("\nreplicates per experiment (V29):")
print(reps.to_string())
assert (chosen.groupby("experiment").size() >= 2).all(), "every experiment should have 2 replicates in V29"
# chosen = pd.concat(chosen)
missing = sorted(set(EXPS) - set(chosen["experiment"]))
print(f"\nchosen: {len(chosen)} files from {chosen['experiment'].nunique()} experiments; "
      f"experiments with no GRCh38 file: {missing}")

os.makedirs(RAW, exist_ok=True)
for r in chosen.itertuples():
    out = os.path.join(RAW, f"{r.experiment}_{r.file}.tsv")
    if not os.path.exists(out):
        urllib.request.urlretrieve("https://www.encodeproject.org" + r.href, out)
    print("  have", out)
chosen.to_csv("decomposition/data/interim/encode_k562_files.csv", index=False)
print("\nsaved decomposition/data/interim/encode_k562_files.csv")   