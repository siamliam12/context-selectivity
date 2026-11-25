CREATE TABLE g AS
SELECT variant, "group", avg(score) AS m
FROM 'decomposition/data/interim/fig_scores_long.parquet'
GROUP BY variant, "group";

CREATE TABLE z AS
SELECT variant, "group", (m - mu) / sd AS z
FROM (SELECT variant, "group", m,
             avg(m) OVER (PARTITION BY variant) AS mu,
             stddev_pop(m) OVER (PARTITION BY variant) AS sd
      FROM g)
WHERE sd > 0;

COPY (
  SELECT t.label, z."group", avg(z.z) AS trait_z, count(*) AS n_variants
  FROM z JOIN read_csv_auto('decomposition/data/frozen/fig4c_reconstructed_set.csv') AS t USING (variant)
  GROUP BY t.label, z."group"
  ORDER BY t.label, z."group"
) TO 'decomposition/results/d1/table_recon_mean_sql.csv' (HEADER);