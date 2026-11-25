# Gate 1 lead: lists the authors' processed_borzoi_parquets folder and shows one file's columns.
import gcsfs
import pyarrow.parquet as pq

fs = gcsfs.GCSFileSystem(token="anon")
base = "seqnn-share/sniff/processed_borzoi_parquets/"
files = fs.find(base, detail=True)
print(f"{len(files)} files, {sum(i['size'] for i in files.values()) / 1e9:.1f} GB in total")
for p, i in list(files.items())[:15]:
    print(f"{i['size']:>15,}  {p}")
first = next(p for p in files if p.endswith(".parquet"))
with fs.open(first, "rb") as f:
    pf = pq.ParquetFile(f)
    print(f"\n{first}: {pf.metadata.num_rows:,} rows")
    print(pf.schema_arrow)