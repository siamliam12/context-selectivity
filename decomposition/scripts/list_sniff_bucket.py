# Lists file names and sizes in the Sniff bucket. Opens no files.
import gcsfs

fs = gcsfs.GCSFileSystem(token="anon")

print("== top folder ==")
for p in fs.ls("seqnn-share/sniff/"):
    print(p)

target = "seqnn-share/sniff/sniff_results_ukb_15/"
if fs.exists(target):
    print("\n== sniff_results_ukb_15 (all files) ==")
    for path, info in fs.find(target, detail=True).items():
        print(f"{info['size']:>15,}  {path}")
else:
    print("\nsniff_results_ukb_15/ not found")
