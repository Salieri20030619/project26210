#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 13_revalidate_gse242299_defective_members_v1.sh REBUILD_ROOT PROXY_URL}
proxy_url=${2:?Usage: 13_revalidate_gse242299_defective_members_v1.sh REBUILD_ROOT PROXY_URL}
extract_root="$rebuild_root/02_complete_cohort_build/source_extract_v1.part/GSE242299"
final_dir="$rebuild_root/01_input_inventory/gse242299_defect_revalidation_v1"
part_dir="$final_dir.part"
result_file="$part_dir/direct_vs_tar_member_validation_v1.tsv"
log_file="$rebuild_root/logs/13_revalidate_gse242299_defective_members_v1.log"

exec > >(tee "$log_file") 2>&1
echo "started_at=$(date --iso-8601=seconds)"
[[ ! -e "$final_dir" ]] || { echo "refusing existing final directory: $final_dir" >&2; exit 1; }
[[ ! -e "$part_dir" ]] || { echo "refusing existing part directory: $part_dir" >&2; exit 1; }
mkdir "$part_dir"
echo -e "geo_accession\tlibrary_id\tmember\tdirect_size\ttar_size\tdirect_sha256\ttar_sha256\tbyte_identical\tgzip_status" > "$result_file"

for specification in \
  "GSM7757833 T2_1 barcodes.tsv.gz" \
  "GSM7757833 T2_1 features.tsv.gz" \
  "GSM7757833 T2_1 matrix.mtx.gz" \
  "GSM7757855 N3_1 barcodes.tsv.gz" \
  "GSM7757855 N3_1 features.tsv.gz" \
  "GSM7757855 N3_1 matrix.mtx.gz"; do
  read -r gsm library suffix <<< "$specification"
  member="${gsm}_${library}_${suffix}"
  url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7757nnn/${gsm}/suppl/${member}"
  direct="$part_dir/$member"
  tar_member="$extract_root/$member"
  curl --silent --show-error --proxy "$proxy_url" --fail --location --output "$direct" "$url"
  [[ -f "$tar_member" ]] || { echo "missing extracted tar member: $tar_member" >&2; exit 1; }
  direct_size=$(stat -c %s "$direct")
  tar_size=$(stat -c %s "$tar_member")
  direct_sha=$(sha256sum "$direct" | cut -d ' ' -f 1)
  tar_sha=$(sha256sum "$tar_member" | cut -d ' ' -f 1)
  identical=false
  cmp -s "$direct" "$tar_member" && identical=true
  gzip_status=PASS
  gzip -t "$direct" 2>/dev/null || gzip_status=FAIL
  echo -e "$gsm\t$library\t$member\t$direct_size\t$tar_size\t$direct_sha\t$tar_sha\t$identical\t$gzip_status" >> "$result_file"
done

if grep -q $'\tfalse\t' "$result_file"; then
  echo "direct and tar member mismatch" >&2
  exit 1
fi
n3_zero=$(awk -F '\t' '$2=="N3_1" && ($3 ~ /barcodes/ || $3 ~ /features/) && $4==0 {n++} END {print n+0}' "$result_file")
[[ "$n3_zero" -eq 2 ]] || { echo "N3_1 zero-byte expectation not confirmed" >&2; exit 1; }
{
  echo "status=PASS_CONFIRMED_OFFICIAL_MEMBER_DEFECTS"
  echo "direct_members=6"
  echo "byte_identical_to_series_tar=6"
  echo "N3_1_zero_byte_members=2"
  echo "T2_1_semantic_matrix_failure=confirmed_by_raw_count_structure_audit_v3"
} > "$part_dir/GSE242299_OFFICIAL_MEMBER_DEFECTS_CONFIRMED_v1"
mv "$part_dir" "$final_dir"
sha256sum "$final_dir/direct_vs_tar_member_validation_v1.tsv" "$final_dir/GSE242299_OFFICIAL_MEMBER_DEFECTS_CONFIRMED_v1"
echo "finished_at=$(date --iso-8601=seconds)"
