# C2 step 3: lists the released Borzoi score files and shows what is inside the standardized one.
# Prints names, shapes and types. For one-dimensional lists (variant names, positions, letters,
# track names) it prints the first 5 entries. It never reads the 2-D score matrix.
import gcsfs, h5py

fs = gcsfs.GCSFileSystem(token="anon")
folder = "seqnn-share/sniff/borzoi_veps/"

print("== borzoi_veps folder ==")
paths = fs.find(folder, detail=True)
for path, info in paths.items():
    print(f"{info['size']:>18,}  {path}")

matches = [p for p in paths if p.endswith("distmatched.h5")]
assert matches, "no distmatched .h5 file found in the folder"
target = matches[0]

print(f"\n== inside {target} ==")
with fs.open(target, "rb") as f, h5py.File(f, "r") as h:
    for k, val in h.attrs.items():
        print(f"file attribute {k}: {val}")
    def show(name, obj):
        if isinstance(obj, h5py.Dataset):
            line = f"{name}: shape {obj.shape}, type {obj.dtype}"
            if obj.ndim == 1:
                line += f", first 5: {obj[:5].tolist()}"
            print(line)
        else:
            print(f"{name}/ (group)")
    h.visititems(show)