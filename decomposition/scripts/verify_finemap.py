# Compares each downloaded file's size with the bucket. Reports missing or wrong-size files.
import gcsfs, os

fs = gcsfs.GCSFileSystem(token="anon")
src = "seqnn-share/sniff/sniff_results_ukb_15/"
dst = "decomposition/data/raw/sniff/finemap/"

n, bad = 0, 0
for path, info in fs.find(src, detail=True).items():
    if not path.endswith(".txt.gz"):
        continue
    n += 1
    out = os.path.join(dst, path[len(src):])
    local = os.path.getsize(out) if os.path.exists(out) else -1
    if local != info["size"]:
        bad += 1
        print("MISMATCH", path[len(src):], local, info["size"])
print(f"{n} files in bucket, {bad} missing or wrong size")
