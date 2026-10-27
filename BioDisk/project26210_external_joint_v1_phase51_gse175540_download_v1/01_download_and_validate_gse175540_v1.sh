#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: bash 01_download_and_validate_gse175540_v1.sh <final_download_dir> <http_proxy> <validator.py>" >&2
  exit 2
fi

final_dir="$1"
proxy_url="$2"
validator="$3"
part_dir="${final_dir}.part"
base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE175nnn/GSE175540"
expected_archive_bytes=1122877440

if [[ -e "$final_dir" ]]; then
  echo "Refusing to overwrite existing final directory: $final_dir" >&2
  exit 3
fi
if [[ ! -f "$validator" ]]; then
  echo "Validator not found: $validator" >&2
  exit 4
fi

mkdir -p "$part_dir/downloads"

download_resumable() {
  local url="$1"
  local destination="$2"
  local partial="${destination}.part"
  curl --fail --location --retry 20 --retry-delay 5 --retry-all-errors \
    --proxy "$proxy_url" --continue-at - --output "$partial" "$url"
  mv "$partial" "$destination"
}

download_resumable "$base_url/suppl/filelist.txt" "$part_dir/downloads/filelist.txt"
download_resumable "$base_url/miniml/GSE175540_family.xml.tgz" "$part_dir/downloads/GSE175540_family.xml.tgz"
download_resumable "$base_url/suppl/GSE175540_RAW.tar" "$part_dir/downloads/GSE175540_RAW.tar"

observed_archive_bytes="$(stat -c '%s' "$part_dir/downloads/GSE175540_RAW.tar")"
if [[ "$observed_archive_bytes" -ne "$expected_archive_bytes" ]]; then
  echo "Archive length mismatch: observed=$observed_archive_bytes expected=$expected_archive_bytes" >&2
  exit 5
fi

python3 "$validator" "$part_dir"
mv "$part_dir" "$final_dir"
echo "status=PASS"
echo "final_dir=$final_dir"
