#!/bin/bash

mkdir -p dataset/gtex_data
cd dataset/gtex_data || exit

declare -A gtex_downloads=(
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000176/QTD000176.all.tsv.gz"]="GTEx_Brain.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz"]="GTEx_Blood.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000281/QTD000281.all.tsv.gz"]="GTEx_Muscle.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000116/QTD000116.all.tsv.gz"]="GTEx_Fat.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000316/QTD000316.all.tsv.gz"]="GTEx_Skin.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000216/QTD000216.all.tsv.gz"]="GTEx_FibroBlast.tsv.gz"
)

for url in "${!gtex_downloads[@]}";do
	original_file=$(basename "$url")
	new_file="${gtex_downloads[$url]}"
	
	echo "Starting download for: $original_file"
	if wget -q --show-progress "$url"; then
		echo "[SUCCESS] Download Complete: $original_file"
		mv "$original_file" "$new_file"
		echo "[RENAMED] $original_file Saved as -> $new_file"
		echo "###################################"
	else
		echo "[FAILED] Could not download: $original_file"
		echo "##################################"
	fi
done

cd ../..
mkdir -p dataset/independent
cd dataset/independent || exit

declare -A independent_downloads=(
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000025/QTD000434/QTD000434.all.tsv.gz"]="ROSMAP_Brain.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000019/QTD000373/QTD000373.all.tsv.gz"]="Lepik_Blood.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000011/QTD000085/QTD000085.all.tsv.gz"]="FUSION_Muscle.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000029/QTD000534/QTD000534.all.tsv.gz"]="TWINSUK_Fat.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000029/QTD000544/QTD000544.all.tsv.gz"]="TWINSUK_Skin.tsv.gz"
["https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000012/QTD000100/QTD000100.all.tsv.gz"]="GENCORD_FibroBlast.tsv.gz"
)

for url in "${!independent_downloads[@]}"; do
    original_file=$(basename "$url")
    new_file="${independent_downloads[$url]}"
    
    echo "Starting download for: $original_file"
    
    if wget -q --show-progress "$url"; then
        echo "[SUCCESS] Download complete: $original_file"
        mv "$original_file" "$new_file"
        echo "[RENAMED] Saved as -> $new_file"
        echo "----------------------------------------"
    else
        echo "[FAILED] Could not download: $original_file"
        echo "----------------------------------------"
    fi
done

cd ..
mkdir -p dataset/metadata
cd dataset/metadata || exit

wget -q --show-progress "https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata_r8.tsv"
wget -q --show-progress "https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md"
wget -q --show-progress "https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/tabix_ftp_paths.tsv"

