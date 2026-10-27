#!/usr/bin/env bash
set -euo pipefail

root="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1"
raw_root="$root/inputs_v1/GSE207493/raw_extracted_v1"
extract_marker="$raw_root/EXTRACTION_COMPLETE.txt"
signature="$root/inputs_v1/frozen_signature_genes_v1.tsv"
script="$root/scripts_v1/19_run_gse207493_cell_source_validation_v1.py"
output="$root/results_v1/GSE207493"
python_bin="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
minimum_available_kib=$((80 * 1024 * 1024))

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export VECLIB_MAXIMUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

mkdir -p "$output" "$root/logs_v1"
echo $$ > "$root/logs_v1/gse207493_validation_v1.pid"

if [[ -s "$output/gse207493_data_gate_v1.tsv" ]]; then
  echo "[$(date -Is)] validation outputs already exist"
  exit 0
fi

echo "[$(date -Is)] waiting for safe extraction"
for _ in $(seq 1 720); do
  [[ -s "$extract_marker" ]] && break
  sleep 30
done
[[ -s "$extract_marker" ]] || { echo "[$(date -Is)] extraction wait timed out" >&2; exit 3; }

echo "[$(date -Is)] extraction ready; waiting for at least 80 GiB available RAM"
for _ in $(seq 1 720); do
  available_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  if [[ "$available_kib" -ge "$minimum_available_kib" ]]; then
    break
  fi
  sleep 30
done
available_kib=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
[[ "$available_kib" -ge "$minimum_available_kib" ]] || {
  echo "[$(date -Is)] insufficient available RAM after wait: ${available_kib} KiB" >&2
  exit 4
}

echo "[$(date -Is)] starting 8-thread cell-source validation; MemAvailable=${available_kib} KiB"
"$python_bin" "$script" "$raw_root" "$signature" "$output"
echo "[$(date -Is)] GSE207493 validation complete"
