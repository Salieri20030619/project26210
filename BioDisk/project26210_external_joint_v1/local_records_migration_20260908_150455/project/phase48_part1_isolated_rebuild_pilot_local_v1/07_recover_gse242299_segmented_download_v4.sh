#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 07_recover_gse242299_segmented_download_v4.sh REBUILD_ROOT PROXY_URL SUBCONNECTIONS}
proxy_url=${2:?Usage: 07_recover_gse242299_segmented_download_v4.sh REBUILD_ROOT PROXY_URL SUBCONNECTIONS}
subconnections=${3:?Usage: 07_recover_gse242299_segmented_download_v4.sh REBUILD_ROOT PROXY_URL SUBCONNECTIONS}
url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242299/suppl/GSE242299_all_cells_50236_33538.h5ad.gz
expected_size=5350598000
prefix_size=194220032
original_chunk_size=644547246
output_root="$rebuild_root/01_input_inventory/public_geo_alternatives_v1"
prefix_file="$output_root/GSE242299_all_cells_50236_33538.h5ad.gz.part"
v3_dir="$output_root/GSE242299_segmented_download_v3"
final_file="$output_root/GSE242299_all_cells_50236_33538.h5ad.gz"
work_dir="$output_root/GSE242299_segmented_recovery_v4"
assembly_file="$work_dir/GSE242299_all_cells_50236_33538.h5ad.gz.assembly.part"
plan_file="$work_dir/recovery_subrange_plan_v4.tsv"
source_file="$work_dir/original_range_source_validation_v4.tsv"
log_file="$rebuild_root/logs/07_recover_gse242299_segmented_download_v4.log"

exec > >(tee "$log_file") 2>&1
echo "started_at=$(date --iso-8601=seconds)"
[[ "$subconnections" =~ ^[1-9][0-9]*$ ]] || { echo "invalid subconnection count" >&2; exit 1; }
[[ $(stat -c %s "$prefix_file") -eq "$prefix_size" ]] || { echo "prefix size changed" >&2; exit 1; }
[[ -d "$v3_dir" ]] || { echo "missing preserved v3 directory" >&2; exit 1; }
[[ ! -e "$final_file" ]] || { echo "refusing existing final file" >&2; exit 1; }
[[ ! -e "$work_dir" ]] || { echo "refusing existing v4 work directory" >&2; exit 1; }
mkdir "$work_dir"
echo -e "original_range\trecovery_part\tstart_byte\tend_byte\texpected_bytes\tobserved_bytes\thttp_code\tsha256\tstatus" > "$plan_file"

pids=()
for original_range in 0 1 2 3 4 5 6 7; do
  original_start=$((prefix_size + original_range * original_chunk_size))
  original_end=$((original_start + original_chunk_size - 1))
  (( original_end < expected_size )) || original_end=$((expected_size - 1))
  original_expected=$((original_end - original_start + 1))
  v3_file="$v3_dir/range_${original_range}.part"
  [[ -f "$v3_file" ]] || { echo "missing v3 range: $v3_file" >&2; exit 1; }
  observed=$(stat -c %s "$v3_file")
  (( observed <= original_expected )) || { echo "oversized v3 range $original_range" >&2; exit 1; }
  if (( observed == original_expected )); then
    continue
  fi
  remaining=$((original_expected - observed))
  subchunk=$(((remaining + subconnections - 1) / subconnections))
  for ((part=0; part<subconnections; part++)); do
    start=$((original_start + observed + part * subchunk))
    (( start <= original_end )) || break
    end=$((start + subchunk - 1))
    (( end <= original_end )) || end=$original_end
    expected=$((end - start + 1))
    file="$work_dir/range_${original_range}_recovery_${part}.part"
    meta="$work_dir/range_${original_range}_recovery_${part}.curl.tsv"
    stderr_file="$work_dir/range_${original_range}_recovery_${part}.curl.stderr.log"
    max_size=$((expected + 1024))
    echo -e "$original_range\t$part\t$start\t$end\t$expected\tPENDING\tPENDING\tPENDING\tPENDING" >> "$plan_file"
    curl --silent --show-error --proxy "$proxy_url" --fail --location --retry 5 --retry-delay 10 \
      --max-filesize "$max_size" --range "$start-$end" --output "$file" \
      --write-out '%{http_code}\t%{size_download}\n' "$url" > "$meta" 2> "$stderr_file" &
    pids+=("$!")
  done
done

download_failed=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    download_failed=1
  fi
done
[[ "$download_failed" -eq 0 ]] || { echo "one or more v4 recovery downloads failed" >&2; exit 1; }

validated_plan="$work_dir/recovery_subrange_validation_v4.tsv"
echo -e "original_range\trecovery_part\tstart_byte\tend_byte\texpected_bytes\tobserved_bytes\thttp_code\tsha256\tstatus" > "$validated_plan"
while IFS=$'\t' read -r original_range part start end expected _rest; do
  [[ "$original_range" = "original_range" ]] && continue
  file="$work_dir/range_${original_range}_recovery_${part}.part"
  meta="$work_dir/range_${original_range}_recovery_${part}.curl.tsv"
  observed=$(stat -c %s "$file")
  http_code=$(cut -f 1 "$meta")
  sha256=$(sha256sum "$file" | cut -d ' ' -f 1)
  status=PASS
  [[ "$observed" -eq "$expected" && "$http_code" = 206 ]] || status=FAIL
  echo -e "$original_range\t$part\t$start\t$end\t$expected\t$observed\t$http_code\t$sha256\t$status" >> "$validated_plan"
done < "$plan_file"
if grep -q FAIL "$validated_plan"; then
  echo "v4 recovery subrange validation failed" >&2
  exit 1
fi

echo -e "original_range\texpected_bytes\tv3_prefix_bytes\trecovered_bytes\tcompleted_bytes\tsha256\tstatus" > "$source_file"
for original_range in 0 1 2 3 4 5 6 7; do
  original_start=$((prefix_size + original_range * original_chunk_size))
  original_end=$((original_start + original_chunk_size - 1))
  (( original_end < expected_size )) || original_end=$((expected_size - 1))
  expected=$((original_end - original_start + 1))
  v3_file="$v3_dir/range_${original_range}.part"
  v3_bytes=$(stat -c %s "$v3_file")
  recovered=0
  source="$v3_file"
  if (( v3_bytes < expected )); then
    source="$work_dir/range_${original_range}.completed_v4.part"
    cp "$v3_file" "$source"
    for recovery_file in "$work_dir/range_${original_range}_recovery_"*.part; do
      cat "$recovery_file" >> "$source"
      recovered=$((recovered + $(stat -c %s "$recovery_file")))
    done
  fi
  completed=$(stat -c %s "$source")
  [[ "$completed" -eq "$expected" ]] || { echo "completed original range size mismatch: $original_range" >&2; exit 1; }
  sha256=$(sha256sum "$source" | cut -d ' ' -f 1)
  echo -e "$original_range\t$expected\t$v3_bytes\t$recovered\t$completed\t$sha256\tPASS" >> "$source_file"
done

cp "$prefix_file" "$assembly_file"
for original_range in 0 1 2 3 4 5 6 7; do
  completed="$work_dir/range_${original_range}.completed_v4.part"
  if [[ -f "$completed" ]]; then
    cat "$completed" >> "$assembly_file"
  else
    cat "$v3_dir/range_${original_range}.part" >> "$assembly_file"
  fi
done
observed_size=$(stat -c %s "$assembly_file")
[[ "$observed_size" -eq "$expected_size" ]] || { echo "assembled size mismatch: $observed_size" >&2; exit 1; }
gzip -t "$assembly_file"
full_sha256=$(sha256sum "$assembly_file" | cut -d ' ' -f 1)
mv "$assembly_file" "$final_file"
{
  echo -e "file\tsize_bytes\tsha256\turl\tprefix_size\tv3_reused_bytes\tv4_recovered_bytes\tstatus"
  v3_reused=$(awk -F '\t' 'NR>1 {s+=$3} END {print s+0}' "$source_file")
  v4_recovered=$(awk -F '\t' 'NR>1 {s+=$4} END {print s+0}' "$source_file")
  echo -e "$(basename "$final_file")\t$observed_size\t$full_sha256\t$url\t$prefix_size\t$v3_reused\t$v4_recovered\tPASS"
} > "$output_root/GSE242299_official_h5ad_download_manifest_v4.tsv"
echo "final_size=$observed_size"
echo "full_sha256=$full_sha256"
echo "finished_at=$(date --iso-8601=seconds)"
