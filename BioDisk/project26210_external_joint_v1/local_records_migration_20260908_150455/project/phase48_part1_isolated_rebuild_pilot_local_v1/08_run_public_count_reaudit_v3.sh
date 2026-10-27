#!/usr/bin/env bash
set -euo pipefail

rebuild_root=${1:?Usage: 08_run_public_count_reaudit_v3.sh REBUILD_ROOT PYTHON_BIN}
python_bin=${2:?Usage: 08_run_public_count_reaudit_v3.sh REBUILD_ROOT PYTHON_BIN}
extract_root="$rebuild_root/02_complete_cohort_build/source_extract_v1.part"
scope_file="$rebuild_root/00_scope_and_provenance/full_cohort_library_scope_v2.tsv"
script_file="$rebuild_root/scripts/08_audit_public_geo_counts_v3.py"
output_dir="$rebuild_root/02_complete_cohort_build/raw_count_source_audit_v3"
log_file="$rebuild_root/logs/08_run_public_count_reaudit_v3.log"

exec > >(tee "$log_file") 2>&1
echo "started_at=$(date --iso-8601=seconds)"
[[ -d "$extract_root" ]] || { echo "missing preserved extraction: $extract_root" >&2; exit 1; }
[[ -f "$scope_file" ]] || { echo "missing scope: $scope_file" >&2; exit 1; }
[[ ! -e "$output_dir" ]] || { echo "refusing existing output: $output_dir" >&2; exit 1; }

"$python_bin" "$script_file" \
  --extract-root "$extract_root" \
  --scope "$scope_file" \
  --output "$output_dir"

sha256sum \
  "$output_dir/raw_count_structure_audit_v3.tsv" \
  "$output_dir/raw_count_structure_summary_v3.tsv" \
  "$output_dir/PUBLIC_GEO_RAW_COUNT_STRUCTURE_BLOCKED_v3"
echo "finished_at=$(date --iso-8601=seconds)"
