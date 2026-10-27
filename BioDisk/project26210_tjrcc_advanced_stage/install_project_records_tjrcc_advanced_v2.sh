#!/usr/bin/env bash
set -euo pipefail

project_root="/mnt/BioDisk/ludingyu/project26210"
stage_root="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage/records_update2"
backup_root="${project_root}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2/project_records_before_update_20260803_1045"

declare -A expected_source=(
  [PROJECT_STATUS.md]="46446cc1e0ad387b0db51bf4def4424d63e8942b50909902fb001106896b2430"
  [CODEX_HANDOFF.md]="8c39803ec2430c51694b7c3b2e51ded601c4099e4f92b34566d17559ebaeba9d"
  [RESULTS_MANIFEST.tsv]="b426e80fd742a5e690a77c4dd75abbb0c05acdc5a8a022622130d9f1b4201aa9"
)

for name in PROJECT_STATUS.md CODEX_HANDOFF.md RESULTS_MANIFEST.tsv; do
  current_sha="$(sha256sum "${project_root}/${name}" | awk '{print $1}')"
  if [[ "${current_sha}" != "${expected_source[${name}]}" ]]; then
    echo "ABORT: ${name} changed since staging: ${current_sha}" >&2
    exit 20
  fi
  test -s "${stage_root}/${name}"
done

mkdir -p "${backup_root}"
for name in PROJECT_STATUS.md CODEX_HANDOFF.md RESULTS_MANIFEST.tsv; do
  cp -p "${project_root}/${name}" "${backup_root}/${name}"
done
for name in PROJECT_STATUS.md CODEX_HANDOFF.md RESULTS_MANIFEST.tsv; do
  install -m 0644 "${stage_root}/${name}" "${project_root}/${name}"
  installed_sha="$(sha256sum "${project_root}/${name}" | awk '{print $1}')"
  staged_sha="$(sha256sum "${stage_root}/${name}" | awk '{print $1}')"
  [[ "${installed_sha}" == "${staged_sha}" ]]
  printf '%s\t%s\n' "${installed_sha}" "${project_root}/${name}"
done
echo "backup=${backup_root}"
