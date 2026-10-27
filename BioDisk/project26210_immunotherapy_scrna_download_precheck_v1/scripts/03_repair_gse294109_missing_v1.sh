#!/usr/bin/env bash
set -euo pipefail

base=/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1
filelist="$base/metadata/GSE294109_filelist.txt"
files="$base/downloads/GSE294109/files"
pre="$base/precheck"

mkdir -p "$files" "$pre"

download_one() {
    local url="$1"
    local filename="${url##*/}"
    curl -fL --retry 20 --retry-all-errors --retry-delay 3 \
        -C - -o "$files/$filename" "$url"
}
export -f download_one
export files

for pass in $(seq 1 10); do
    missing="$pre/GSE294109_missing_or_incomplete_pass_${pass}.tsv"
    : > "$missing"

    while IFS=$'\t' read -r filename expected; do
        actual=0
        if [[ -f "$files/$filename" ]]; then
            actual=$(stat -c %s "$files/$filename")
        fi
        if [[ "$actual" -ne "$expected" ]]; then
            gsm=${filename%%_*}
            gsm_group="${gsm:0:${#gsm}-3}nnn"
            url="https://ftp.ncbi.nlm.nih.gov/geo/samples/$gsm_group/$gsm/suppl/$filename"
            printf '%s\t%s\t%s\t%s\n' "$filename" "$expected" "$actual" "$url" >> "$missing"
        fi
    done < <(awk -F '\t' '$1 == "File" {print $2 "\t" $4}' "$filelist")

    count=$(wc -l < "$missing")
    printf 'pass=%s incomplete=%s\n' "$pass" "$count"
    if [[ "$count" -eq 0 ]]; then
        printf 'All GSE294109 files match the official byte sizes.\n'
        exit 0
    fi

    cut -f 4 "$missing" | xargs -n 1 -P 8 bash -c 'download_one "$1"' _ || true
    sleep 5
done

printf 'GSE294109 still has incomplete files after 10 repair passes.\n' >&2
exit 1
