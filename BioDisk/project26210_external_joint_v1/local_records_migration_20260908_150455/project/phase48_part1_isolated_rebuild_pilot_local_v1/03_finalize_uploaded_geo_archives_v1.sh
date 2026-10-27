#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "usage: $0 ISOLATED_ROOT VALIDATION_TSV" >&2
  exit 64
fi

isolated_root=$1
validation_tsv=$2
if [[ "$(basename -- "$isolated_root")" != "part1_dataset_scope_rebuild_pilot_v1" ]]; then
  echo "refusing unexpected isolated root: $isolated_root" >&2
  exit 65
fi
test -f "$validation_tsv"

archive_dir="$isolated_root/01_input_inventory/public_geo_archives_v1"
tmp_dir="$isolated_root/tmp"
manifest_final="$archive_dir/public_geo_archive_manifest_v2.tsv"
manifest_part="$tmp_dir/public_geo_archive_manifest_v2.tsv.part"
versions_final="$archive_dir/finalize_tool_versions_v1.txt"
versions_part="$tmp_dir/finalize_tool_versions_v1.txt.part"
mkdir -p "$archive_dir"
if [[ -e "$manifest_final" || -e "$versions_final" ]]; then
  echo "refusing to overwrite existing finalization output" >&2
  exit 66
fi

printf 'study\tarchive_path\texpected_bytes\tactual_bytes\texpected_sha256\tactual_sha256\texpected_tar_members\tactual_tar_members\tstatus\n' > "$manifest_part"

while IFS=$'\t' read -r study filename expected_bytes expected_sha256 expected_tar_members transfer_method; do
  [[ "$study" == study ]] && continue
  part_path="$tmp_dir/$filename.local_upload_v1.part"
  final_path="$archive_dir/$filename"
  member_final="$archive_dir/${study}_tar_members_v2.txt"
  member_part="$tmp_dir/${study}_tar_members_v2.txt.part"
  if [[ -e "$final_path" ]]; then
    test ! -e "$part_path"
    source_path="$final_path"
    finalize_status=preexisting_validated
  else
    test -f "$part_path"
    source_path="$part_path"
    finalize_status=finalized_from_local_upload
  fi

  actual_bytes=$(stat -c '%s' -- "$source_path")
  actual_sha256=$(sha256sum -- "$source_path" | awk '{print $1}')
  tar -tf "$source_path" > "$member_part"
  actual_tar_members=$(wc -l < "$member_part" | awk '{print $1}')
  if [[ "$actual_bytes" -ne "$expected_bytes" || "$actual_sha256" != "$expected_sha256" || "$actual_tar_members" -ne "$expected_tar_members" ]]; then
    echo "validation mismatch for $study" >&2
    exit 67
  fi
  if [[ -e "$member_final" ]]; then
    echo "refusing to overwrite member list: $member_final" >&2
    exit 68
  fi
  mv -- "$member_part" "$member_final"
  if [[ "$source_path" == "$part_path" ]]; then
    mv -- "$part_path" "$final_path"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$study" "$final_path" "$expected_bytes" "$actual_bytes" "$expected_sha256" "$actual_sha256" \
    "$expected_tar_members" "$actual_tar_members" "$finalize_status" >> "$manifest_part"
done < "$validation_tsv"

row_count=$(tail -n +2 "$manifest_part" | wc -l | awk '{print $1}')
[[ "$row_count" -eq 3 ]]
{
  printf 'timestamp='; date -Is
  printf 'hostname='; hostname
  printf 'bash='; bash --version | head -n 1
  printf 'tar='; tar --version | head -n 1
  printf 'sha256sum='; sha256sum --version | head -n 1
} > "$versions_part"
mv -- "$manifest_part" "$manifest_final"
mv -- "$versions_part" "$versions_final"
cat "$manifest_final"
