#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 ISOLATED_ROOT GSE178481_RDS_DIR GSE207493_RDS_DIR" >&2
  exit 64
fi

isolated_root=$1
gse178481_dir=$2
gse207493_dir=$3

if [[ "$(basename -- "$isolated_root")" != "part1_dataset_scope_rebuild_pilot_v1" ]]; then
  echo "refusing unexpected isolated root: $isolated_root" >&2
  exit 65
fi

inventory_dir="$isolated_root/01_input_inventory"
tmp_dir="$isolated_root/tmp"
test -d "$inventory_dir"
test -d "$tmp_dir"
test -d "$gse178481_dir"
test -d "$gse207493_dir"

inventory_part="$tmp_dir/input_file_inventory_v2.tsv.part"
summary_part="$tmp_dir/preflight_summary_v2.tsv.part"
versions_part="$tmp_dir/command_versions_v2.txt.part"
inventory_final="$inventory_dir/input_file_inventory_v2.tsv"
summary_final="$inventory_dir/preflight_summary_v2.tsv"
versions_final="$inventory_dir/command_versions_v2.txt"

for final_path in "$inventory_final" "$summary_final" "$versions_final"; do
  if [[ -e "$final_path" ]]; then
    echo "refusing to overwrite existing output: $final_path" >&2
    exit 66
  fi
done

printf 'dataset\tinput_path\tbytes\tsha256\n' > "$inventory_part"

inventory_dataset() {
  local dataset=$1
  local source_dir=$2
  local filename_pattern=$3
  while IFS= read -r -d '' input_path; do
    local bytes
    local checksum
    bytes=$(stat -c '%s' -- "$input_path")
    checksum=$(sha256sum -- "$input_path" | awk '{print $1}')
    printf '%s\t%s\t%s\t%s\n' "$dataset" "$input_path" "$bytes" "$checksum" >> "$inventory_part"
  done < <(find "$source_dir" -maxdepth 1 -type f -name "$filename_pattern" -print0 | sort -z)
}

inventory_dataset GSE178481 "$gse178481_dir" '*.rds'
inventory_dataset GSE207493 "$gse207493_dir" 'GSE207493__*.rds'

gse178481_n=$(awk -F '\t' '$1 == "GSE178481" {n++} END {print n+0}' "$inventory_part")
gse207493_n=$(awk -F '\t' '$1 == "GSE207493" {n++} END {print n+0}' "$inventory_part")
total_n=$((gse178481_n + gse207493_n))
unique_path_n=$(tail -n +2 "$inventory_part" | cut -f2 | sort -u | wc -l | awk '{print $1}')
unique_hash_n=$(tail -n +2 "$inventory_part" | cut -f4 | sort -u | wc -l | awk '{print $1}')

status=PASS
if [[ "$gse178481_n" -ne 26 || "$gse207493_n" -ne 19 || "$total_n" -ne 45 || "$unique_path_n" -ne 45 ]]; then
  status=FAIL
fi

printf 'check\tvalue\texpected\tstatus\n' > "$summary_part"
printf 'gse178481_rds_count\t%s\t26\t%s\n' "$gse178481_n" "$([[ "$gse178481_n" -eq 26 ]] && echo PASS || echo FAIL)" >> "$summary_part"
printf 'gse207493_rds_count\t%s\t19\t%s\n' "$gse207493_n" "$([[ "$gse207493_n" -eq 19 ]] && echo PASS || echo FAIL)" >> "$summary_part"
printf 'total_rds_count\t%s\t45\t%s\n' "$total_n" "$([[ "$total_n" -eq 45 ]] && echo PASS || echo FAIL)" >> "$summary_part"
printf 'unique_input_paths\t%s\t45\t%s\n' "$unique_path_n" "$([[ "$unique_path_n" -eq 45 ]] && echo PASS || echo FAIL)" >> "$summary_part"
printf 'unique_sha256\t%s\t45\tINFO\n' "$unique_hash_n" >> "$summary_part"
printf 'overall_status\t%s\tPASS\t%s\n' "$status" "$status" >> "$summary_part"

{
  printf 'timestamp='; date -Is
  printf 'hostname='; hostname
  printf 'bash='; bash --version | head -n 1
  printf 'sha256sum='; sha256sum --version | head -n 1
  printf 'stat='; stat --version | head -n 1
  printf 'find='; find --version | head -n 1
} > "$versions_part"

mv -- "$inventory_part" "$inventory_final"
mv -- "$summary_part" "$summary_final"
mv -- "$versions_part" "$versions_final"

cat "$summary_final"
[[ "$status" == PASS ]]
