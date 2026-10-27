#!/usr/bin/env bash
set -euo pipefail

base=/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1
meta="$base/metadata"
pre="$base/precheck"
g264="$base/downloads/GSE264586"
g294="$base/downloads/GSE294109"

mkdir -p "$pre" "$g264/extracted" "$g294/files"

report="$pre/precheck_report_v1.txt"
: > "$report"

printf 'Project 26210 immunotherapy scRNA-seq download precheck\n' >> "$report"
date -Is >> "$report"
printf '\n[GSE264586 archive]\n' >> "$report"
stat -c 'path=%n\nsize_bytes=%s' "$g264/GSE264586_RAW.tar" >> "$report"
tar -tf "$g264/GSE264586_RAW.tar" > "$pre/GSE264586_tar_members.txt"
printf 'tar_member_count=%s\n' "$(wc -l < "$pre/GSE264586_tar_members.txt")" >> "$report"
printf 'barcodes_count=%s\n' "$(grep -c '\.barcodes\.tsv\.gz$' "$pre/GSE264586_tar_members.txt")" >> "$report"
printf 'features_count=%s\n' "$(grep -c '\.features\.tsv\.gz$' "$pre/GSE264586_tar_members.txt")" >> "$report"
printf 'matrix_count=%s\n' "$(grep -c '\.matrix\.mtx\.gz$' "$pre/GSE264586_tar_members.txt")" >> "$report"

if [ ! -f "$g264/extracted/.extract_complete" ]; then
    tar -xf "$g264/GSE264586_RAW.tar" -C "$g264/extracted"
    touch "$g264/extracted/.extract_complete"
fi

find "$g264/extracted" -type f -name '*.gz' -print0 | xargs -0 -n 1 -P 8 gzip -t
printf 'gzip_test=PASS\n' >> "$report"

printf '\n[GSE264586 patient metadata]\n' >> "$report"
awk -F '\t' 'NR>1 {p++; samples+=$3; if($12=="Responder") r++; if($12=="Nonresponder") nr++} END {printf "patients=%d\nsamples=%d\nresponders=%d\nnonresponders=%d\n",p,samples,r,nr}' "$meta/GSE264586_patient_response_v1.tsv" >> "$report"

printf '\n[GSE294109 individual files]\n' >> "$report"
printf 'file_count=%s\n' "$(find "$g294/files" -maxdepth 1 -type f | wc -l)" >> "$report"
printf 'barcodes_count=%s\n' "$(find "$g294/files" -maxdepth 1 -type f -name '*.barcodes.tsv.gz' | wc -l)" >> "$report"
printf 'features_count=%s\n' "$(find "$g294/files" -maxdepth 1 -type f -name '*.features.tsv.gz' | wc -l)" >> "$report"
printf 'matrix_count=%s\n' "$(find "$g294/files" -maxdepth 1 -type f -name '*.matrix.mtx.gz' | wc -l)" >> "$report"

find "$g294/files" -maxdepth 1 -type f -name '*.gz' -print0 | xargs -0 -n 1 -P 8 gzip -t
printf 'gzip_test=PASS\n' >> "$report"

printf '\n[GSE294109 patient metadata]\n' >> "$report"
awk -F '\t' 'NR>1 {p++; libs+=$3; if($8=="Responder") r++; if($8=="Nonresponder") nr++} END {printf "patients=%d\nlibraries=%d\nICB_responders=%d\nICB_nonresponders=%d\n",p,libs,r,nr}' "$meta/GSE294109_patient_metadata_v1.tsv" >> "$report"

printf '\n[Matrix Market headers]\n' >> "$report"
first264=$(find "$g264/extracted" -type f -name '*.matrix.mtx.gz' | sort | head -n 1)
first294=$(find "$g294/files" -maxdepth 1 -type f -name '*.matrix.mtx.gz' | sort | head -n 1)
printf 'GSE264586_example=%s\n' "$first264" >> "$report"
gzip -cd "$first264" | sed -n '1,3p' >> "$report"
printf 'GSE294109_example=%s\n' "$first294" >> "$report"
gzip -cd "$first294" | sed -n '1,3p' >> "$report"

find "$base/downloads" -type f -printf '%P\t%s\n' | sort > "$pre/download_manifest_sizes_v1.tsv"
du -sh "$g264" "$g294" >> "$report"

