#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 ISOLATED_ROOT SOURCE_TSV" >&2
  exit 64
fi

isolated_root=$1
source_tsv=$2

if [[ "$(basename -- "$isolated_root")" != "part1_dataset_scope_rebuild_pilot_v1" ]]; then
  echo "refusing unexpected isolated root: $isolated_root" >&2
  exit 65
fi
test -f "$source_tsv"

inventory_root="$isolated_root/01_input_inventory"
archive_dir="$inventory_root/public_geo_archives_v1"
tmp_dir="$isolated_root/tmp"
manifest_final="$archive_dir/public_geo_archive_manifest_v1.tsv"
manifest_part="$tmp_dir/public_geo_archive_manifest_v1.tsv.part"
versions_final="$archive_dir/download_tool_versions_v1.txt"
versions_part="$tmp_dir/download_tool_versions_v1.txt.part"

test -d "$inventory_root"
test -d "$tmp_dir"
mkdir -p "$archive_dir"
if [[ -e "$manifest_final" || -e "$versions_final" ]]; then
  echo "refusing to overwrite existing download manifest or version record" >&2
  exit 66
fi

printf 'study\turl\tarchive_path\texpected_bytes\tactual_bytes\tsha256\ttar_members\tdownload_status\tvalidation_status\n' > "$manifest_part"

while IFS=$'\t' read -r study url expected_bytes public_scope selection_policy; do
  [[ "$study" == study ]] && continue
  [[ -n "$study" && -n "$url" && -n "$expected_bytes" ]]
  filename=$(basename -- "$url")
  final_path="$archive_dir/$filename"
  part_path="$tmp_dir/$filename.part"
  member_final="$archive_dir/${study}_tar_members_v1.txt"
  member_part="$tmp_dir/${study}_tar_members_v1.txt.part"
  if [[ -e "$member_final" ]]; then
    echo "refusing to overwrite existing member list: $member_final" >&2
    exit 67
  fi

  download_status=downloaded
  if [[ -e "$final_path" ]]; then
    download_status=preexisting
  else
    curl --fail --location --retry 5 --retry-delay 3 --retry-all-errors \
      --connect-timeout 30 --continue-at - --output "$part_path" "$url"
    actual_part_bytes=$(stat -c '%s' -- "$part_path")
    if [[ "$actual_part_bytes" -ne "$expected_bytes" ]]; then
      echo "downloaded size mismatch for $study: $actual_part_bytes != $expected_bytes" >&2
      exit 68
    fi
    mv -- "$part_path" "$final_path"
  fi

  actual_bytes=$(stat -c '%s' -- "$final_path")
  if [[ "$actual_bytes" -ne "$expected_bytes" ]]; then
    echo "final size mismatch for $study: $actual_bytes != $expected_bytes" >&2
    exit 69
  fi
  checksum=$(sha256sum -- "$final_path" | awk '{print $1}')
  tar -tf "$final_path" > "$member_part"
  tar_members=$(wc -l < "$member_part" | awk '{print $1}')
  if [[ "$tar_members" -lt 1 ]]; then
    echo "empty tar archive for $study" >&2
    exit 70
  fi
  mv -- "$member_part" "$member_final"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "$study" "$url" "$final_path" "$expected_bytes" "$actual_bytes" "$checksum" \
    "$tar_members" "$download_status" >> "$manifest_part"
done < "$source_tsv"

row_count=$(tail -n +2 "$manifest_part" | wc -l | awk '{print $1}')
if [[ "$row_count" -ne 3 ]]; then
  echo "archive manifest row count mismatch: $row_count != 3" >&2
  exit 71
fi

{
  printf 'timestamp='; date -Is
  printf 'hostname='; hostname
  printf 'bash='; bash --version | head -n 1
  printf 'curl='; curl --version | head -n 1
  printf 'tar='; tar --version | head -n 1
  printf 'sha256sum='; sha256sum --version | head -n 1
} > "$versions_part"

mv -- "$manifest_part" "$manifest_final"
mv -- "$versions_part" "$versions_final"
cat "$manifest_final"
