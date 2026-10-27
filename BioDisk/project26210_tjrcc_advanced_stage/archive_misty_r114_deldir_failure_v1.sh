#!/usr/bin/env bash
set -euo pipefail

run_root="/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
stage_root="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage"
failed_status="${run_root}/misty/samples/R114_T/status.json"
failed_log="${run_root}/logs/misty/R114_T.log"
audit_dir="${run_root}/preflight/misty_r114_deldir_failure_v1"

test -s "${failed_status}"
test -s "${failed_log}"
rg -q '"status": "failed"' "${failed_status}"
rg -q 'Bailing out of trifnd' "${failed_status}"
rg -q 'Both vertex orderings are clockwise' "${failed_log}"
mkdir -p "${audit_dir}"
cp -p "${failed_status}" "${audit_dir}/R114_T_failed_status_v1.json"
cp -p "${failed_log}" "${audit_dir}/R114_T_failed_run_v1.log"
cp -p "${stage_root}/misty_r114_deldir_failure_audit_v1.md" "${audit_dir}/"
cp -p "${run_root}/misty/inputs/R114_T_misty_positions_v1.tsv" "${audit_dir}/"
cp -p "${run_root}/scripts/run_tjrcc_misty_sample_v1.R" "${audit_dir}/run_tjrcc_misty_sample_retry_v1.R"
sha256sum "${audit_dir}"/* >"${audit_dir}/sha256_manifest_v1.tsv"
echo "MISTY_R114_DELDIR_FAILURE_ARCHIVED audit_dir=${audit_dir}"
