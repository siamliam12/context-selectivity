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

-- Step 2: Join with fat
CREATE TABLE step2 AS
SELECT a.variant,a.gene
FROM step1 a
INNER JOIN 'GTEx_Fat.parquet' b
	ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step2;

-- Step 3: Join with GTEx_FibroBlast
CREATE TABLE step3 AS
SELECT a.variant,a.gene
FROM step2 a
INNER JOIN 'GTEx_FibroBlast.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step3;

-- Step 4: Join with GTEx_Muscle
CREATE TABLE step4 AS
SELECT a.variant,a.gene
FROM step3 a
INNER JOIN 'GTEx_Muscle.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step4;


-- Step 5: Join with GTEx_Skin
CREATE TABLE step5 AS
SELECT a.variant,a.gene
FROM step4 a
INNER JOIN 'GTEx_Skin.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step5;


-- Step 6: Join with Lepik_Blood
CREATE TABLE step6 AS
SELECT a.variant,a.gene
FROM step5 a
INNER JOIN 'Lepik_Blood.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step6;


-- Step 7: Join with ROSMAP_Brain
CREATE TABLE step7 AS
SELECT a.variant,a.gene
FROM step6 a
INNER JOIN 'ROSMAP_Brain.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step7;


-- Step 8: Join with TWINSUK_Fat
CREATE TABLE step8 AS
SELECT a.variant,a.gene
FROM step7 a
INNER JOIN 'TWINSUK_Fat.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step8;


-- Step 9: Join with GENCORD_FibroBlast
CREATE TABLE step9 AS
SELECT a.variant,a.gene
FROM step8 a
INNER JOIN 'GENCORD_FibroBlast.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step9;


-- Step 10: Join with FUSION_Muscle
CREATE TABLE step10 AS
SELECT a.variant,a.gene
FROM step9 a
INNER JOIN 'FUSION_Muscle.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step10;


-- Step 11: Join with TWINSUK_Skin
CREATE TABLE step11 AS
SELECT a.variant,a.gene
FROM step10 a
INNER JOIN 'TWINSUK_Skin.parquet' b
        ON a.variant = b.variant AND a.gene = b.gene;
SELECT count(*) FROM step11;


COPY step11 TO '../frozen/frozen_variant_set.tsv' (HEADER, DELIMITER '\t');

















