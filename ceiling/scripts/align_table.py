from datetime import datetime
from pathlib import Path
 
import duckdb
import numpy as np
import pandas as pd

FROZEN = Path("data/frozen/frozen_variant_set.tsv")
INTERIM = Path("data/interim")
OUT = Path("data/frozen/aligned.parquet")
LOG = Path("logs/c3_aligned.log")

FILES = [
    "GTEx_Blood",
    "GTEx_Brain",
    "GTEx_Fat", 
    "GTEx_FibroBlast",
    "GTEx_Muscle", 
    "GTEx_Skin", 
    "Lepik_Blood",  
    "ROSMAP_Brain",  
    "TWINSUK_Fat",  
    "GENCORD_FibroBlast",  
    "FUSION_Muscle",  
    "TWINSUK_Skin",
]
STATS = ["beta", "se", "pvalue"]

def log(msg):
    print(msg)
    with LOG.open("a") as fh:
        fh.write(msg + "\n")
 
def read_frozen():
    # dtype=str + keep_default_na=False: IDs stay exact text; nothing becomes NaN.
    df = pd.read_csv(FROZEN, sep="\t", dtype=str, keep_default_na=False)
    assert {"variant", "gene"} <= set(df.columns), f"frozen columns are {list(df.columns)}"
    return df[["variant", "gene"]].reset_index(drop=True)

# ---------- 0. guards ----------
assert len(FILES) == 12 and len(set(FILES)) == 12, "FILES must list 12 distinct IDs"
assert not OUT.exists(), f"{OUT} already exists; delete it on purpose if a rebuild is intended"
LOG.parent.mkdir(parents=True, exist_ok=True)
log(f"--- C3 run {datetime.now().isoformat(timespec='seconds')} ---")

# ---------- 1. frozen set + seat numbers ----------
frozen = read_frozen()
n = len(frozen)
assert n >= 50_000, f"frozen set has only {n} pairs"
assert not frozen.duplicated().any(), "frozen set contains duplicate (variant, gene) pairs"
assert (frozen != "").all().all(), "frozen set contains empty IDs"
log(f"frozen pairs: {n}")
 
keys = frozen.copy()
keys["idx"] = np.arange(n, dtype=np.int64)  # seat number = position in the frozen file
 
con = duckdb.connect()  # in-memory; nothing here needs to persist
con.register("frozen_keys", keys)

# ---------- 2. pull beta / se / pvalue from each file, back in frozen order ----------
aligned = frozen.copy()
for qtd in FILES:
    path = INTERIM / f"{qtd}.parquet"
    assert path.exists(), f"missing {path}"
    part = con.execute(f"""
        SELECT k.idx, p.beta, p.se, p.pvalue
        FROM frozen_keys k
        JOIN read_parquet('{path.as_posix()}') p
          ON p.variant = k.variant AND p.gene = k.gene
        ORDER BY k.idx
    """).df()
 
    # every frozen pair matched exactly once: right count AND seats exactly 0..n-1
    assert len(part) == n, f"{qtd}: {len(part)} rows returned, expected {n}"
    assert np.array_equal(part["idx"].to_numpy(), np.arange(n)), f"{qtd}: pairs missing or duplicated"

    # values are real numbers in legal ranges (single-file checks)
    for s in STATS:
        assert part[s].dtype == np.float64, f"{qtd}.{s} has dtype {part[s].dtype}"
    assert np.isfinite(part[STATS].to_numpy()).all(), f"{qtd}: missing or non-finite values"
    assert (part["se"] > 0).all(), f"{qtd}: se <= 0 present"
    assert part["pvalue"].between(0, 1).all(), f"{qtd}: pvalue outside [0, 1]"
 
    for s in STATS:
        aligned[f"{qtd}_{s}"] = part[s].to_numpy()  # positional; safe because seats == 0..n-1
    log(f"{qtd}: {len(part)} rows matched, checks passed")

# ---------- 3. whole-table checks ----------
expected_cols = ["variant", "gene"] + [f"{q}_{s}" for q in FILES for s in STATS]
assert list(aligned.columns) == expected_cols, "column order is wrong"
assert aligned.shape == (n, 38), f"shape is {aligned.shape}"
assert not aligned.isna().any().any(), "missing values in aligned table"

# ---------- 4. save, then re-check what is actually on disk ----------
aligned.to_parquet(OUT, index=False)
 
back = pd.read_parquet(OUT)
again = read_frozen()  # fresh read of the frozen file, independent of everything above
assert back.shape == (n, 38), f"reloaded shape is {back.shape}"
assert list(back.columns) == expected_cols, "reloaded columns differ"
assert not back.isna().any().any(), "reloaded table has missing values"
assert back["variant"].tolist() == again["variant"].tolist(), "variant order differs from frozen file"
assert back["gene"].tolist() == again["gene"].tolist(), "gene order differs from frozen file"
 
log(f"wrote {OUT}: {back.shape[0]} rows x {back.shape[1]} columns")
log("columns: " + ", ".join(expected_cols))
log("C3 checks passed")