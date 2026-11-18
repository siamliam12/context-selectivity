# D1 check: the Python table function vs the DuckDB table, set A, SD by n. Compares numbers only.
import pandas as pd
from table_fn import build_table

long = pd.read_parquet("decomposition/data/interim/fig_scores_long.parquet")
vt = pd.read_csv("decomposition/data/frozen/variant_traits.csv")
py, dropped = build_table(long, vt[vt["direction"] == "promoted"], ddof=0)
sql = pd.read_csv("decomposition/results/d1/table_setA_sql.csv")

m = py.merge(sql, on=["label", "group"], suffixes=("_py", "_sql"), how="outer", indicator=True)
print(f"cells: python {len(py)}, sql {len(sql)}, in both {int((m['_merge'] == 'both').sum())}")
print(f"variants dropped for zero spread: {dropped}")
print(f"largest difference in median Z: {(m['median_z_py'] - m['median_z_sql']).abs().max():.2e}")
print(f"variant counts identical: {bool((m['n_variants_py'] == m['n_variants_sql']).all())}")
py.to_csv("decomposition/results/d1/table_setA_py.csv", index=False)