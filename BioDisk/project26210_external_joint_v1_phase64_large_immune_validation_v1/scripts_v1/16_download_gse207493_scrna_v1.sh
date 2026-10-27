#!/usr/bin/env bash
set -euo pipefail

root="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1"
download_dir="$root/downloads_v1/GSE207493"
log_dir="$root/logs_v1"
url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207493/suppl/GSE207493_scRNA-seq.tar.gz"
proxy_url="http://127.0.0.1:7897"
partial="$download_dir/GSE207493_scRNA-seq.tar.gz.part"
final="$download_dir/GSE207493_scRNA-seq.tar.gz"

mkdir -p "$download_dir" "$log_dir"
echo "$$" > "$log_dir/gse207493_scrna_download_v1.pid"
date -Is > "$log_dir/gse207493_scrna_download_v1.started_at"

if [[ -e "$final" ]]; then
  echo "Final archive already exists; validating it"
else
  curl --http1.1 --proxy "$proxy_url" --fail --location \
    --retry 30 --retry-all-errors --continue-at - --user-agent "Mozilla/5.0" \
    --output "$partial" "$url"
  tar -tzf "$partial" > "$log_dir/gse207493_scrna_tar_inventory_v1.txt"
  mv "$partial" "$final"
fi

tar -tzf "$final" > "$log_dir/gse207493_scrna_tar_inventory_v1.txt"
sha256sum "$final" > "$final.sha256"
stat -c "bytes=%s" "$final"
echo "status=PASS"
date -Is > "$log_dir/gse207493_scrna_download_v1.completed_at"
