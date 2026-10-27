#!/usr/bin/env bash
set -euo pipefail

base=/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1
dataset_dir="$base/downloads/GSE294109"
files_dir="$dataset_dir/files"
metadata_dir="$base/metadata"
log_dir="$base/logs"

mkdir -p "$files_dir" "$metadata_dir" "$log_dir"

filelist_url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294109/suppl/filelist.txt
curl -fsSL --retry 10 --retry-delay 2 "$filelist_url" -o "$metadata_dir/GSE294109_filelist.txt"

awk -F '\t' '$1 == "File" {print $2}' "$metadata_dir/GSE294109_filelist.txt" | while IFS= read -r filename; do
    gsm=${filename%%_*}
    gsm_group="${gsm:0:${#gsm}-3}nnn"
    printf 'https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s\n' "$gsm_group" "$gsm" "$filename"
done > "$metadata_dir/GSE294109_individual_urls.txt"

cd "$files_dir"
xargs -n 1 -P 12 wget -c --tries=20 --timeout=90 --no-verbose < "$metadata_dir/GSE294109_individual_urls.txt"

find "$files_dir" -maxdepth 1 -type f -printf '%f\t%s\n' | sort > "$metadata_dir/GSE294109_downloaded_files_sizes.tsv"

