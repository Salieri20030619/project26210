#!/usr/bin/env bash
set -euo pipefail

run_root="/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
stage_root="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage"
samples=(S15_T Y12_T Y27_T)

for sample in "${samples[@]}"; do
  status="${run_root}/misty/samples/${sample}/status.json"
  log="${run_root}/logs/misty/${sample}.log"
  audit_dir="${run_root}/preflight/misty_${sample,,}_deldir_failure_v1"
  test -s "${status}"
  test -s "${log}"
  rg -q '"status": "failed"' "${status}"
  rg -q 'Bailing out of trifnd' "${status}"
  rg -q 'Both vertex orderings are clockwise' "${log}"
  mkdir -p "${audit_dir}"
  cp -p "${status}" "${audit_dir}/${sample}_failed_status_v1.json"
  cp -p "${log}" "${audit_dir}/${sample}_failed_run_v1.log"
  cp -p "${stage_root}/misty_deldir_failure_audit_template_v1.md" "${audit_dir}/misty_deldir_failure_audit_v1.md"
  cp -p "${run_root}/misty/inputs/${sample}_misty_positions_v1.tsv" "${audit_dir}/"
  cp -p "${run_root}/scripts/run_tjrcc_misty_sample_v1.R" "${audit_dir}/run_tjrcc_misty_sample_retry_v1.R"
  sha256sum "${audit_dir}"/* >"${audit_dir}/sha256_manifest_v1.tsv"
  echo "ARCHIVED sample=${sample} audit_dir=${audit_dir}"
done
echo "MISTY_REMAINING_DELDIR_FAILURES_ARCHIVED samples=${#samples[@]}"
