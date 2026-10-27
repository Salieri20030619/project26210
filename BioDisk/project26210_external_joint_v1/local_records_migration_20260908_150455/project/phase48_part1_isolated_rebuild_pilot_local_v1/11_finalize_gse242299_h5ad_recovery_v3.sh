#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 11_finalize_gse242299_h5ad_recovery_v3.sh REBUILD_ROOT PYTHON_BIN}
python_bin=${2:?Usage: 11_finalize_gse242299_h5ad_recovery_v3.sh REBUILD_ROOT PYTHON_BIN}
source_h5ad="$rebuild_root/02_complete_cohort_build/gse242299_h5ad_recovery_v2.part/GSE242299_all_cells_50236_33538.h5ad"
compressed_h5ad="$rebuild_root/01_input_inventory/public_geo_alternatives_v1/GSE242299_all_cells_50236_33538.h5ad.gz"
download_manifest="$rebuild_root/01_input_inventory/public_geo_alternatives_v1/GSE242299_official_h5ad_download_manifest_v4.tsv"
scope_file="$rebuild_root/00_scope_and_provenance/full_cohort_library_scope_v2.tsv"
audit_script="$rebuild_root/scripts/11_inventory_gse242299_h5ad_v3.py"
final_dir="$rebuild_root/02_complete_cohort_build/gse242299_h5ad_recovery_v3"
part_dir="$final_dir.part"
part_h5ad="$part_dir/GSE242299_all_cells_50236_33538.h5ad"
log_file="$rebuild_root/logs/11_finalize_gse242299_h5ad_recovery_v3.log"

exec > >(tee "$log_file") 2>&1
echo "started_at=$(date --iso-8601=seconds)"
[[ -f "$source_h5ad" ]] || { echo "missing preserved v2 decompressed H5AD" >&2; exit 1; }
[[ -f "$compressed_h5ad" ]] || { echo "missing validated compressed H5AD" >&2; exit 1; }
[[ -f "$download_manifest" ]] || { echo "missing v4 download manifest" >&2; exit 1; }
[[ ! -e "$final_dir" ]] || { echo "refusing existing final directory" >&2; exit 1; }
[[ ! -e "$part_dir" ]] || { echo "refusing existing part directory" >&2; exit 1; }
mkdir "$part_dir"

cp --reflink=auto "$source_h5ad" "$part_h5ad"
[[ $(stat -c %s "$part_h5ad") -eq 7345987312 ]] || { echo "decompressed H5AD size mismatch" >&2; exit 1; }
"$python_bin" "$audit_script" --h5ad "$part_h5ad" --scope "$scope_file" --output "$part_dir"

decompressed_sha256=$(sha256sum "$part_h5ad" | cut -d ' ' -f 1)
compressed_sha256=$(awk -F '\t' 'NR==2 {print $3}' "$download_manifest")
{
  echo -e "decompressed_file\tdecompressed_size_bytes\tdecompressed_sha256\tcompressed_file\tcompressed_size_bytes\tcompressed_sha256\tstatus"
  echo -e "$(basename "$part_h5ad")\t$(stat -c %s "$part_h5ad")\t$decompressed_sha256\t$(basename "$compressed_h5ad")\t$(stat -c %s "$compressed_h5ad")\t$compressed_sha256\tPASS"
} > "$part_dir/H5AD_SOURCE_MANIFEST_v3.tsv"

mv "$part_dir" "$final_dir"
sha256sum \
  "$final_dir/h5ad_matrix_candidate_audit_v3.tsv" \
  "$final_dir/h5ad_library_cell_counts_v3.tsv" \
  "$final_dir/h5ad_recovery_metadata_v3.json" \
  "$final_dir/GSE242299_H5AD_RECOVERY_PASS_v3" \
  "$final_dir/H5AD_SOURCE_MANIFEST_v3.tsv"
echo "finished_at=$(date --iso-8601=seconds)"
