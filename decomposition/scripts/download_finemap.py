# Downloads the Sniff and baseline fine-mapping files. Skips files already downloaded.
import gcsfs, os

fs = gcsfs.GCSFileSystem(token="anon")
src = "seqnn-share/sniff/sniff_results_ukb_15/"
dst = "decomposition/data/raw/sniff/finemap/"

for path in fs.find(src):
    if not path.endswith(".txt.gz"):
        continue
    rel = path[len(src):]
    out = os.path.join(dst, rel)
    if os.path.exists(out):
        continue
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fs.get(path, out)
    print(rel, flush=True)
