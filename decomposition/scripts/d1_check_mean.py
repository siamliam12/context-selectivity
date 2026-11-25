# D1 re-check after Amendment 1: Python vs DuckDB, mean over the reconstructed set; and the old median check still holds.
import pandas as pd
from table_fn import build_table

long = pd.read_parquet("decomposition/data/interim/fig_scores_long.parquet")
rec = pd.read_csv("decomposition/data/frozen/fig4c_reconstructed_set.csv")
vt = pd.read_csv("decomposition/data/frozen/variant_traits.csv")

def compare(py, sql_path, label):
    sql = pd.read_csv(sql_path).rename(columns={"median_z": "trait_z"})
    m = py.merge(sql, on=["label", "group"], suffixes=("_py", "_sql"), how="outer", indicator=True)
    print(f"{label}: cells in both {int((m['_merge'] == 'both').sum())} of {len(py)}; "
          f"largest difference {(m['trait_z_py'] - m['trait_z_sql']).abs().max():.1e}; "
          f"counts identical {bool((m['n_variants_py'] == m['n_variants_sql']).all())}")

py_mean, d1 = build_table(long, rec, ddof=0, agg="mean")
compare(py_mean, "decomposition/results/d1/table_recon_mean_sql.csv", "mean, reconstructed set")
py_med, d2 = build_table(long, vt[vt["direction"] == "promoted"], ddof=0, agg="median")
compare(py_med, "decomposition/results/d1/table_setA_sql.csv", "median, set A (old check)")
py_mean.to_csv("decomposition/results/d1/table_recon_mean_py.csv", index=False)