import duckdb
from pathlib import Path 

con = duckdb.connect("data/interim/ceiling.duckdb")

sumstats_dir = Path("data/raw/sumstats")
tsv_files = list(sumstats_dir.glob("*/*.tsv.gz"))

for file_path in tsv_files:
    base_name = file_path.name.replace(".tsv.gz","")
    output_parquet = f"data/interim/{base_name}.parquet"

    if Path(output_parquet).exists():
        print(f"Skipping {file_path} (Already processed)...")
        continue

    print(f"Processing {file_path} -> {output_parquet}")

    query = f"""
    COPY(
        SELECT DISTINCT variant,
            regexp_replace(molecular_trait_id,'\\.[0-9]+$','') AS gene,
            chromosome,position,ref,alt,beta,se,pvalue,an
        FROM read_csv('{file_path}',delim='\t',header=true,nullstr='NA',types={{'chromosome':'VARCHAR'}})
    ) TO '{output_parquet}' (FORMAT parquet)
    """
    con.sql(query)

print("All files processed successfully!")
