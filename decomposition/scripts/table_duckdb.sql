CREATE TABLE g AS
SELECT variant,"group",avg(score) AS m
FROM 'decomposition/data/interim/fig_scores_long.parquet'
GROUP BY variant,"group";

CREATE TABLE z AS
SELECT variant,"group",(m-mu)/sd AS z
FROM(
    SELECT variant,"group",m,
            avg(m) OVER (PARTITION BY variant) as mu,
            stddev_pop(m) OVER (PARTITION BY variant) as sd
    FROM g 
)
WHERE sd>0;

COPY (
    SELECT t.label,z."group",median(z.z) AS median_z,count(*) AS n_variants
    FROM z
    JOIN read_Csv_auto('decomposition/data/frozen/variant_traits.csv') AS t USING (variant)
    WHERE t.direction = 'promoted'
    GROUP BY t.label, z."group"
    ORDER BY t.label, z."group"
) TO 'decomposition/results/d1/table_setA_sql.csv' (HEADER);