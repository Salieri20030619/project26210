#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 07_download_gse242299_official_h5ad_v3_segmented.sh REBUILD_ROOT PROXY_URL CONNECTIONS}
proxy_url=${2:?Usage: 07_download_gse242299_official_h5ad_v3_segmented.sh REBUILD_ROOT PROXY_URL CONNECTIONS}
connections=${3:?Usage: 07_download_gse242299_official_h5ad_v3_segmented.sh REBUILD_ROOT PROXY_URL CONNECTIONS}
url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242299/suppl/GSE242299_all_cells_50236_33538.h5ad.gz
expected_size=5350598000
output_root="$rebuild_root/01_input_inventory/public_geo_alternatives_v1"
prefix_file="$output_root/GSE242299_all_cells_50236_33538.h5ad.gz.part"
final_file="$output_root/GSE242299_all_cells_50236_33538.h5ad.gz"
work_dir="$output_root/GSE242299_segmented_download_v3"
assembly_file="$work_dir/GSE242299_all_cells_50236_33538.h5ad.gz.assembly.part"
log_file="$rebuild_root/logs/07_download_gse242299_official_h5ad_v3_segmented.log"
plan_file="$work_dir/range_plan_and_validation_v3.tsv"

exec > >(tee "$log_file") 2>&1
echo "started_at=$(date --iso-8601=seconds)"
echo "proxy_url=$proxy_url"
echo "connections=$connections"
[[ "$connections" =~ ^[1-9][0-9]*$ ]] || { echo "invalid connection count" >&2; exit 1; }
[[ -f "$prefix_file" ]] || { echo "missing prefix: $prefix_file" >&2; exit 1; }
[[ ! -e "$final_file" ]] || { echo "refusing existing final file: $final_file" >&2; exit 1; }
[[ ! -e "$work_dir" ]] || { echo "refusing existing work directory: $work_dir" >&2; exit 1; }
mkdir "$work_dir"

prefix_size=$(stat -c %s "$prefix_file")
prefix_sha256=$(sha256sum "$prefix_file" | cut -d ' ' -f 1)
[[ "$prefix_size" -gt 0 && "$prefix_size" -lt "$expected_size" ]] || { echo "invalid prefix size: $prefix_size" >&2; exit 1; }
remaining=$((expected_size - prefix_size))
chunk_size=$(((remaining + connections - 1) / connections))
echo -e "range_id\tstart_byte\tend_byte\texpected_bytes\tobserved_bytes\thttp_code\tsha256\tstatus" > "$plan_file"
echo "prefix_size=$prefix_size"
echo "prefix_sha256=$prefix_sha256"
echo "remaining_bytes=$remaining"
echo "chunk_size=$chunk_size"

pids=()
range_ids=()
for ((range_id=0; range_id<connections; range_id++)); do
  start=$((prefix_size + range_id * chunk_size))
  (( start < expected_size )) || break
  end=$((start + chunk_size - 1))
  (( end < expected_size )) || end=$((expected_size - 1))
  expected=$((end - start + 1))
  file="$work_dir/range_${range_id}.part"
  meta="$work_dir/range_${range_id}.curl.tsv"
  stderr_file="$work_dir/range_${range_id}.curl.stderr.log"
  max_size=$((expected + 1024))
  curl --silent --show-error --proxy "$proxy_url" --fail --location --retry 5 --retry-delay 10 \
    --max-filesize "$max_size" --range "$start-$end" --output "$file" \
    --write-out '%{http_code}\t%{size_download}\n' "$url" > "$meta" 2> "$stderr_file" &
  pids+=("$!")
  range_ids+=("$range_id")
done

download_failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    download_failed=1
  fi
done
[[ "$download_failed" -eq 0 ]] || { echo "one or more range downloads failed" >&2; exit 1; }

for range_id in "${range_ids[@]}"; do
  start=$((prefix_size + range_id * chunk_size))
  end=$((start + chunk_size - 1))
  (( end < expected_size )) || end=$((expected_size - 1))
  expected=$((end - start + 1))
  file="$work_dir/range_${range_id}.part"
  meta="$work_dir/range_${range_id}.curl.tsv"
  observed=$(stat -c %s "$file")
  http_code=$(cut -f 1 "$meta")
  sha256=$(sha256sum "$file" | cut -d ' ' -f 1)
  status=PASS
  [[ "$observed" -eq "$expected" && "$http_code" = 206 ]] || status=FAIL
  echo -e "$range_id\t$start\t$end\t$expected\t$observed\t$http_code\t$sha256\t$status" >> "$plan_file"
done
if grep -q FAIL "$plan_file"; then
  echo "range validation failed" >&2
  exit 1
fi

cp "$prefix_file" "$assembly_file"
for range_id in "${range_ids[@]}"; do
  cat "$work_dir/range_${range_id}.part" >> "$assembly_file"
done
observed_size=$(stat -c %s "$assembly_file")
[[ "$observed_size" -eq "$expected_size" ]] || { echo "assembled size mismatch: $observed_size" >&2; exit 1; }
gzip -t "$assembly_file"
full_sha256=$(sha256sum "$assembly_file" | cut -d ' ' -f 1)
mv "$assembly_file" "$final_file"
{
  echo -e "file\tsize_bytes\tsha256\turl\tprefix_size\tprefix_sha256\tconnections\tstatus"
  echo -e "$(basename "$final_file")\t$observed_size\t$full_sha256\t$url\t$prefix_size\t$prefix_sha256\t${#range_ids[@]}\tPASS"
} > "$output_root/GSE242299_official_h5ad_download_manifest_v3.tsv"
echo "final_size=$observed_size"
echo "full_sha256=$full_sha256"
echo "finished_at=$(date --iso-8601=seconds)"
