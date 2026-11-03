# Prints S3's sheets, columns, first rows, and any column with few distinct values (like traits).
import pandas as pd

s3 = pd.read_excel("decomposition/data/raw/sniff/media-1.xlsx", sheet_name=None)
for name, df in s3.items():
    print("== sheet:", name, df.shape)
    print(list(df.columns))
    print(df.head(3).to_string())
    for col in df.columns:
        if df[col].nunique() <= 20:
            print(col, sorted(df[col].dropna().astype(str).unique()))
