#!/usr/bin/env bash
set -euo pipefail

STAGE="/mnt/BioDisk/wangmingyuan/project26210_part1_scvi_recovery_stage"
PROJECT="/mnt/BioDisk/ludingyu/project26210"

verify_original() {
  local path="$1"
  local expected="$2"
  local actual
  actual="$(sha256sum "${path}" | awk '{print $1}')"
  if [[ "${actual}" != "${expected}" ]]; then
    echo "source changed since staging: ${path}" >&2
    echo "expected=${expected} actual=${actual}" >&2
    exit 1
  fi
}

install_one() {
  local source="$1"
  local target="$2"
  local temp="${target}.part1_recovery_20260803.tmp"
  cp "${source}" "${temp}"
  cmp -s "${source}" "${temp}"
  mv "${temp}" "${target}"
  cmp -s "${source}" "${target}"
}

verify_original "${PROJECT}/CODEX_HANDOFF.md" "4ea8334d9e2f91469b6b53aa890ae0d19fc7cdf3775026527bc98ef325681a7e"
verify_original "${PROJECT}/PROJECT_STATUS.md" "0bbc7a69b8c04fe63c8846f16522cf157ddc50623150b423e07ff978e2b9a491"
verify_original "${PROJECT}/RESULTS_MANIFEST.tsv" "49ac161fceca70278d52e280180df82832307a574913bb3c95204b58d96f322f"
verify_original "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md" "1ef40c9477e970b24ac473541a604801f6284cac45c8c6cbcd3b7530623fc3cf"

install_one "${STAGE}/CODEX_HANDOFF.md" "${PROJECT}/CODEX_HANDOFF.md"
install_one "${STAGE}/PROJECT_STATUS.md" "${PROJECT}/PROJECT_STATUS.md"
install_one "${STAGE}/RESULTS_MANIFEST.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
install_one "${STAGE}/PART1_REBUILD_STATUS.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"

echo "PART1_RECOVERY_RECORDS_INSTALL_V1_PASS"
sha256sum \
  "${PROJECT}/CODEX_HANDOFF.md" \
  "${PROJECT}/PROJECT_STATUS.md" \
  "${PROJECT}/RESULTS_MANIFEST.tsv" \
  "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
