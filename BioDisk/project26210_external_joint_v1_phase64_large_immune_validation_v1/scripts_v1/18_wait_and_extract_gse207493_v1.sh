#!/usr/bin/env bash
set -euo pipefail

root="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1"
archive="$root/downloads_v1/GSE207493/GSE207493_scRNA-seq.tar.gz"
download_pid_file="$root/logs_v1/gse207493_scrna_download_v1.pid"
target="$root/inputs_v1/GSE207493/raw_extracted_v1"
inventory="$root/inputs_v1/GSE207493/archive_inventory_v1.tsv"
extractor="$root/scripts_v1/17_safe_extract_gse207493_v1.py"
python_bin="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"

mkdir -p "$target" "$(dirname "$inventory")" "$root/logs_v1"

if [[ -s "$target/EXTRACTION_COMPLETE.txt" ]]; then
  echo "[$(date -Is)] extraction already complete"
  exit 0
fi

echo "[$(date -Is)] waiting for validated archive"
for _ in $(seq 1 720); do
  if [[ -s "$archive" ]]; then
    echo "[$(date -Is)] archive ready; starting safe extraction"
    "$python_bin" "$extractor" "$archive" "$target" "$inventory"
    echo "[$(date -Is)] extraction complete"
    exit 0
  fi

  if [[ -s "$download_pid_file" ]]; then
    download_pid=$(cat "$download_pid_file")
    if ! kill -0 "$download_pid" 2>/dev/null; then
      echo "[$(date -Is)] download process ended before archive became ready" >&2
      exit 2
    fi
  fi
  sleep 30
done

echo "[$(date -Is)] timed out waiting for archive" >&2
exit 3
