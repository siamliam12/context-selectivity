import duckdb
from scipy.stats import spearmanr

row = duckdb.sql("""
    SELECT GTEx_Brain_beta,
           ROSMAP_Brain_beta, Lepik_Blood_beta, FUSION_Muscle_beta,
           TWINSUK_Fat_beta, TWINSUK_Skin_beta, GENCORD_FibroBlast_beta
    FROM 'data/frozen/aligned.parquet'
    WHERE GTEx_Brain_pvalue < 1e-5
""").df()
print(len(row))

cols = ["ROSMAP_Brain", "Lepik_Blood", "FUSION_Muscle",
        "TWINSUK_Fat", "TWINSUK_Skin", "GENCORD_FibroBlast"]
rhos = {}
for c in cols:
    rhos[c] = spearmanr(row["GTEx_Brain_beta"], row[c + "_beta"])[0]
    print(f"{c:20s} {rhos[c]:.3f}")

diag = rhos["ROSMAP_Brain"]
off = sum(v for k, v in rhos.items() if k != "ROSMAP_Brain") / 5
print(f"brain row gap = {diag:.3f} - {off:.3f} = {diag - off:.3f}")