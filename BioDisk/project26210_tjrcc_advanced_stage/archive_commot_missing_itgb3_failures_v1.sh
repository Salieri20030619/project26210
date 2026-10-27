#!/usr/bin/env bash
set -euo pipefail

run_root="/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
stage_root="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage"
audit_dir="${run_root}/preflight/commot_missing_itgb3_failure_v1"
samples=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)

mkdir -p "${audit_dir}"
cp -p "${stage_root}/commot_missing_itgb3_failure_audit_v1.md" "${audit_dir}/"
for sample in "${samples[@]}"; do
  status="${run_root}/commot/samples/${sample}/status.json"
  log="${run_root}/logs/commot/${sample}.log"
  test -s "${status}"
  test -s "${log}"
  rg -q '"status": "failed"' "${status}"
  rg -q "Candidate genes missing from raw mapping H5AD: \['ITGB3'\]" "${status}"
  cp -p "${status}" "${audit_dir}/${sample}_failed_status_v1.json"
  cp -p "${log}" "${audit_dir}/${sample}_failed_run_v1.log"
done
cp -p "${run_root}/scripts/run_tjrcc_commot_candidate_v1.py" "${audit_dir}/run_tjrcc_commot_candidate_failed_v1.py"
cp -p "${stage_root}/run_tjrcc_commot_candidate_v1.py" "${audit_dir}/run_tjrcc_commot_candidate_corrected_v1.py"
sha256sum "${audit_dir}"/* >"${audit_dir}/sha256_manifest_v1.tsv"
echo "COMMOT_MISSING_ITGB3_FAILURES_ARCHIVED samples=${#samples[@]} audit_dir=${audit_dir}"
