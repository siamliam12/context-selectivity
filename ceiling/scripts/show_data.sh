#!/bin/bash


if [ -z "$1" ]; then
    echo "Error: No File specified."
    echo "Usage: $0 <filename.parquet>"
    exit 1
fi 

FILE_PATH="$HOME/research/ceiling/data/interim/$1"
if [ ! -f "$FILE_PATH" ]; then
    echo "Error:File '$FILE_PATH' does not exist."
    exit 1
fi

echo "Viewing top 10 rows of $1..."
~/tools/duckdb -c "SELECT * FROM '$FILE_PATH' LIMIT 10;"