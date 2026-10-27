-- STEP 0: Base
CREATE TABLE step0 AS
SELECT variant, gene FROM 'GTEx_Blood.parquet';

-- Step 1: Join with Brain
CREATE TABLE step1 AS
SELECT a.variant, a.gene
FROM step0 a
INNER JOIN 'GTEx_Brain.parquet' b
	ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step1;
