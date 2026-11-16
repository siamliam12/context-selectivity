# C3 step 4e: lists ENCODE's released K562 polyA-plus RNA-seq experiments, so the K562 source
# for the bone_marrow group can be chosen by the rule fixed beforehand. Reads ENCODE's catalogue only.
import json, urllib.request, urllib.parse, urllib.error
import pandas as pd
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 80)

params = [("type", "Experiment"), ("assay_title", "polyA plus RNA-seq"),
          ("biosample_ontology.term_name", "K562"), ("status", "released"),
          ("format", "json"), ("limit", "all")]
for f in ["accession", "assay_title", "biosample_summary", "lab.title", "award.rfa",
          "date_released", "description"]:
    params.append(("field", f))
url = "https://www.encodeproject.org/search/?" + urllib.parse.urlencode(params)
try:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    graph = json.load(urllib.request.urlopen(req))["@graph"]
except urllib.error.HTTPError as e:
    raise SystemExit(f"ENCODE search failed: {e.code} (a 404 here means no experiments matched)")

rows = [{"experiment": g.get("accession"), "summary": g.get("biosample_summary"),
         "lab": (g.get("lab") or {}).get("title"), "phase": (g.get("award") or {}).get("rfa"),
         "released": g.get("date_released"), "description": g.get("description")} for g in graph]
df = pd.DataFrame(rows).sort_values(["phase", "experiment"])
print(f"{len(df)} released K562 polyA-plus RNA-seq experiments")
print(df.to_string(index=False))
df.to_csv("decomposition/data/interim/encode_k562_experiments.csv", index=False)
print("\nsaved decomposition/data/interim/encode_k562_experiments.csv")