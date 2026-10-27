#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 07_download_gse242299_official_h5ad_v2.sh REBUILD_ROOT PROXY_URL}
proxy_url=${2:?Usage: 07_download_gse242299_official_h5ad_v2.sh REBUILD_ROOT PROXY_URL}
url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242299/suppl/GSE242299_all_cells_50236_33538.h5ad.gz
expected_size=5350598000
output_dir="$rebuild_root/01_input_inventory/public_geo_alternatives_v1"
final_file="$output_dir/GSE242299_all_cells_50236_33538.h5ad.gz"
part_file="$final_file.part"
log_file="$rebuild_root/logs/07_download_gse242299_official_h5ad_v2.log"

mkdir -p "$output_dir"
exec > >(tee -a "$log_file") 2>&1
echo "attempt_started=$(date --iso-8601=seconds)"
echo "url=$url"
echo "expected_size=$expected_size"
echo "proxy_url=$proxy_url"
curl --version | head -n 1

[[ ! -e "$final_file" ]] || { echo "refusing existing final file: $final_file" >&2; exit 1; }
curl --proxy "$proxy_url" --fail --location --retry 5 --retry-delay 10 --continue-at - --output "$part_file" "$url"

observed_size=$(stat -c %s "$part_file")
echo "observed_size=$observed_size"
[[ "$observed_size" -eq "$expected_size" ]] || { echo "size mismatch" >&2; exit 1; }
gzip -t "$part_file"
sha256=$(sha256sum "$part_file" | cut -d ' ' -f 1)
mv "$part_file" "$final_file"
{
  echo -e "file\tsize_bytes\tsha256\turl\tstatus"
  echo -e "$(basename "$final_file")\t$observed_size\t$sha256\t$url\tPASS"
} > "$output_dir/GSE242299_official_h5ad_download_manifest_v2.tsv"
echo "sha256=$sha256"
echo "final_file=$final_file"
echo "attempt_finished=$(date --iso-8601=seconds)"
