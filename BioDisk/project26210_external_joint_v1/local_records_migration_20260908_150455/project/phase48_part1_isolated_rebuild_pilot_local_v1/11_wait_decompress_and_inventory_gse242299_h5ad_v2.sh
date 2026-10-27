#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 11_wait_decompress_and_inventory_gse242299_h5ad_v2.sh REBUILD_ROOT PYTHON_BIN}
python_bin=${2:?Usage: 11_wait_decompress_and_inventory_gse242299_h5ad_v2.sh REBUILD_ROOT PYTHON_BIN}
input_gz="$rebuild_root/01_input_inventory/public_geo_alternatives_v1/GSE242299_all_cells_50236_33538.h5ad.gz"
scope_file="$rebuild_root/00_scope_and_provenance/full_cohort_library_scope_v2.tsv"
audit_script="$rebuild_root/scripts/11_inventory_gse242299_h5ad_v2.py"
final_dir="$rebuild_root/02_complete_cohort_build/gse242299_h5ad_recovery_v2"
part_dir="$final_dir.part"
part_h5ad="$part_dir/GSE242299_all_cells_50236_33538.h5ad"
log_file="$rebuild_root/logs/11_wait_decompress_and_inventory_gse242299_h5ad_v2.log"

exec > >(tee "$log_file") 2>&1
echo "watcher_started=$(date --iso-8601=seconds)"
for ((attempt=1; attempt<=720; attempt++)); do
  if [[ -f "$input_gz" ]]; then
    break
  fi
  if (( attempt % 10 == 0 )); then
    echo "waiting_for_validated_download attempt=$attempt time=$(date --iso-8601=seconds)"
  fi
  sleep 60
done
[[ -f "$input_gz" ]] || { echo "validated download did not appear within 12 hours" >&2; exit 1; }
[[ ! -e "$final_dir" ]] || { echo "refusing existing final directory: $final_dir" >&2; exit 1; }
[[ ! -e "$part_dir" ]] || { echo "refusing existing part directory: $part_dir" >&2; exit 1; }
mkdir "$part_dir"

echo "decompression_started=$(date --iso-8601=seconds)"
gzip -dc "$input_gz" > "$part_h5ad"
echo "decompressed_size=$(stat -c %s "$part_h5ad")"
"$python_bin" "$audit_script" --h5ad "$part_h5ad" --scope "$scope_file" --output "$part_dir"
mv "$part_dir" "$final_dir"
sha256sum \
  "$final_dir/h5ad_matrix_candidate_audit_v2.tsv" \
  "$final_dir/h5ad_library_cell_counts_v2.tsv" \
  "$final_dir/h5ad_recovery_metadata_v2.json" \
  "$final_dir/GSE242299_H5AD_RECOVERY_PASS_v2"
echo "finished_at=$(date --iso-8601=seconds)"
