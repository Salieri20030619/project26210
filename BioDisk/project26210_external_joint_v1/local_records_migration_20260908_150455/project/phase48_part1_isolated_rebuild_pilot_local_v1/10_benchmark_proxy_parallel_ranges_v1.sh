#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 10_benchmark_proxy_parallel_ranges_v1.sh REBUILD_ROOT PROXY_URL}
proxy_url=${2:?Usage: 10_benchmark_proxy_parallel_ranges_v1.sh REBUILD_ROOT PROXY_URL}
url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242299/suppl/GSE242299_all_cells_50236_33538.h5ad.gz
output_dir="$rebuild_root/tmp/proxy_parallel_range_benchmark_v1"
result_file="$output_dir/benchmark_result_v1.tsv"
chunk_size=8388608
base_start=1073741824

[[ ! -e "$output_dir" ]] || { echo "refusing existing output: $output_dir" >&2; exit 1; }
mkdir "$output_dir"
echo -e "range_id\tstart_byte\tend_byte\texpected_bytes\tobserved_bytes\thttp_code\tsha256\tstatus" > "$result_file"

SECONDS=0
pids=()
for range_id in 0 1 2 3; do
  start=$((base_start + range_id * chunk_size))
  end=$((start + chunk_size - 1))
  file="$output_dir/range_${range_id}.part"
  meta="$output_dir/range_${range_id}.curl.tsv"
  curl --silent --show-error --proxy "$proxy_url" --fail --location \
    --max-filesize 9437184 --range "$start-$end" --output "$file" \
    --write-out '%{http_code}\t%{size_download}\n' "$url" > "$meta" &
  pids+=("$!")
done
for pid in "${pids[@]}"; do
  wait "$pid"
done
elapsed=$SECONDS

for range_id in 0 1 2 3; do
  start=$((base_start + range_id * chunk_size))
  end=$((start + chunk_size - 1))
  file="$output_dir/range_${range_id}.part"
  meta="$output_dir/range_${range_id}.curl.tsv"
  observed=$(stat -c %s "$file")
  http_code=$(cut -f 1 "$meta")
  sha256=$(sha256sum "$file" | cut -d ' ' -f 1)
  status=PASS
  [[ "$observed" -eq "$chunk_size" && "$http_code" = 206 ]] || status=FAIL
  echo -e "$range_id\t$start\t$end\t$chunk_size\t$observed\t$http_code\t$sha256\t$status" >> "$result_file"
done
if grep -q FAIL "$result_file"; then
  echo "benchmark_status=FAIL"
  exit 1
fi
total_bytes=$((4 * chunk_size))
aggregate_bps=$((total_bytes / (elapsed > 0 ? elapsed : 1)))
{
  echo "benchmark_status=PASS"
  echo "parallel_connections=4"
  echo "total_bytes=$total_bytes"
  echo "elapsed_seconds=$elapsed"
  echo "aggregate_bytes_per_second=$aggregate_bps"
} > "$output_dir/BENCHMARK_PASS_v1"
cat "$output_dir/BENCHMARK_PASS_v1"
