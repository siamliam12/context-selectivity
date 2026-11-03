# Prints the first 5 lines of one small file from each folder.
# Streams only the start of each file; downloads nothing else.
import gcsfs, gzip

fs = gcsfs.GCSFileSystem(token="anon")
for folder in ["baseline", "sniff"]:
    path = f"seqnn-share/sniff/sniff_results_ukb_15/{folder}/biochemistry_Cholesterol.21.txt.gz"
    print(f"== {folder} ==")
    with fs.open(path, "rb") as f, gzip.open(f, "rt") as g:
        for i, line in enumerate(g):
            print(line.rstrip())
            if i >= 4:
                break
