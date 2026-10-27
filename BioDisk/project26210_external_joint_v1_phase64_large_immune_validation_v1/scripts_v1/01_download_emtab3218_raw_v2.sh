#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: bash 01_download_emtab3218_raw_v2.sh <dataset_dir> <proxy_url> <api_manifest_json> <validator.py>" >&2
  exit 2
fi

dataset_dir="$1"
proxy_url="$2"
api_manifest="$3"
validator="$4"
base_url="https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-3218"
archive_names=("E-MTAB-3218.raw.1.zip" "E-MTAB-3218.raw.2.zip" "E-MTAB-3218.raw.3.zip")
expected_bytes=(148148560 152462656 107548363)
extract_final="$dataset_dir/raw_cel_v2"
extract_part="${extract_final}.part"
qc_dir="$dataset_dir/raw_qc_v2"

mkdir -p "$dataset_dir" "$qc_dir"
if [[ ! -f "$api_manifest" || ! -f "$validator" ]]; then
  echo "Required manifest or validator is missing" >&2
  exit 3
fi
if [[ -e "$extract_final" || -e "$extract_part" ]]; then
  echo "Refusing to overwrite existing v2 extraction" >&2
  exit 4
fi

for index in 0 1 2; do
  name="${archive_names[$index]}"
  expected="${expected_bytes[$index]}"
  archive="$dataset_dir/$name"
  archive_part="${archive}.part"
  if [[ -f "$archive" ]]; then
    observed="$(stat -c '%s' "$archive")"
    if [[ "$observed" -ne "$expected" ]]; then
      echo "Existing archive has wrong size: $name observed=$observed expected=$expected" >&2
      exit 5
    fi
  else
    curl --fail --location --retry 20 --retry-delay 5 --retry-all-errors \
      --proxy "$proxy_url" --continue-at - --output "$archive_part" "$base_url/$name"
    observed="$(stat -c '%s' "$archive_part")"
    if [[ "$observed" -ne "$expected" ]]; then
      echo "Archive size mismatch: $name observed=$observed expected=$expected" >&2
      exit 6
    fi
    mv "$archive_part" "$archive"
  fi
  unzip -t "$archive" >/dev/null
done

mkdir -p "$extract_part"
for name in "${archive_names[@]}"; do
  unzip -q -n "$dataset_dir/$name" -d "$extract_part"
done
python3 "$validator" "$api_manifest" "$extract_part" "$qc_dir"
mv "$extract_part" "$extract_final"
sha256sum "${archive_names[@]/#/$dataset_dir/}" > "$dataset_dir/raw_archives_sha256_v2.tsv"
printf 'status\tPASS\narchives\t3\ncel_files\t114\nproxy\t%s\n' "$proxy_url" > "$dataset_dir/raw_download_status_v2.tsv"
echo "PASS_E_MTAB_3218_RAW_DOWNLOAD_V2"
