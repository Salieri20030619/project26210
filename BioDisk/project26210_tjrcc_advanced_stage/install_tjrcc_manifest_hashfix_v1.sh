#!/usr/bin/env bash
set -euo pipefail

project_manifest="/mnt/BioDisk/ludingyu/project26210/RESULTS_MANIFEST.tsv"
staged_manifest="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage/records_update3/RESULTS_MANIFEST.tsv"
backup_manifest="/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2/project_records_before_update_20260803_1049/RESULTS_MANIFEST.tsv"
expected_sha="be7ca5228d45998345faa831beaae1f3df1833a3e3d6c89bb4893ba1bcbd526c"

current_sha="$(sha256sum "${project_manifest}" | awk '{print $1}')"
if [[ "${current_sha}" != "${expected_sha}" ]]; then
  echo "ABORT: RESULTS_MANIFEST.tsv changed since staging: ${current_sha}" >&2
  exit 20
fi
mkdir -p "$(dirname "${backup_manifest}")"
cp -p "${project_manifest}" "${backup_manifest}"
install -m 0644 "${staged_manifest}" "${project_manifest}"
installed_sha="$(sha256sum "${project_manifest}" | awk '{print $1}')"
staged_sha="$(sha256sum "${staged_manifest}" | awk '{print $1}')"
[[ "${installed_sha}" == "${staged_sha}" ]]
printf '%s\t%s\n' "${installed_sha}" "${project_manifest}"
echo "backup=${backup_manifest}"
