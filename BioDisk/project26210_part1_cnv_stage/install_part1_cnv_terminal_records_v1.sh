#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage
PART1_STATUS="${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
MANUAL_DST="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv/malignancy_cnv_summary_verification_v1/part1_malignancy_cnv_manual_figure_review_v1.tsv"

"${PROJECT}/check_memory_gate.sh"

check_base_hash() {
  local path="$1"
  local expected="$2"
  local actual
  actual="$(sha256sum "${path}" | awk '{print $1}')"
  if [[ "${actual}" != "${expected}" ]]
  then
    echo "RECORD_BASE_CHANGED path=${path} expected=${expected} actual=${actual}" >&2
    exit 1
  fi
}

# Abort rather than overwrite if another Part 4/other task updated shared records
# after the snapshots used to compose this Part 1-only update.
check_base_hash "${PROJECT}/PROJECT_STATUS.md" "7858a1c3545ef419f2c4324c36c8a889acca3679a2dc9c2e92f5093f0c6d7581"
check_base_hash "${PROJECT}/CODEX_HANDOFF.md" "89903c683d723e27b9d6d952ac1ba7ca166fbc550753db5faa8af7df28b7627d"
check_base_hash "${PROJECT}/RESULTS_MANIFEST.tsv" "022fdc78bc25a0c333517ff3e6bb5f8d9372fb3d7608043d6e3ad6ca45451243"
check_base_hash "${PART1_STATUS}" "68110e7a1d942198e40a2ce3a16b2ca378fac42a45024f9c88f413f12e4d1eac"

awk -F '\t' 'NF != 18 {print "invalid manifest fields at line " NR ": " NF > "/dev/stderr"; bad=1} END {exit bad}' "${STAGE}/RESULTS_MANIFEST_20260804.tsv"
awk -F '\t' 'NR>1 {if(seen[$1]++) {print "duplicate result_id: " $1 > "/dev/stderr"; bad=1}} END {exit bad}' "${STAGE}/RESULTS_MANIFEST_20260804.tsv"
rg -q 'Part 1 malignancy CNV v2 全队列、终态恢复及保守恶性分层已完成' "${STAGE}/PROJECT_STATUS_20260804.md"
rg -q 'CNV v2终态恢复、保守恶性分层与验收' "${STAGE}/PART1_REBUILD_STATUS_20260804.md"
rg -q 'Part 1 CNV v2终态恢复与保守恶性分层' "${STAGE}/CODEX_HANDOFF_20260804.md"
rg -q '^p1_malignancy_cnv_cell_calls_v1_20260804\t' "${STAGE}/RESULTS_MANIFEST_20260804.tsv"
test -s "${STAGE}/part1_malignancy_cnv_manual_figure_review_v1.tsv"

cp "${STAGE}/PROJECT_STATUS_20260804.md" "${PROJECT}/PROJECT_STATUS.md"
cp "${STAGE}/CODEX_HANDOFF_20260804.md" "${PROJECT}/CODEX_HANDOFF.md"
cp "${STAGE}/RESULTS_MANIFEST_20260804.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cp "${STAGE}/PART1_REBUILD_STATUS_20260804.md" "${PART1_STATUS}"
cp "${STAGE}/part1_malignancy_cnv_manual_figure_review_v1.tsv" "${MANUAL_DST}"

cmp -s "${STAGE}/PROJECT_STATUS_20260804.md" "${PROJECT}/PROJECT_STATUS.md"
cmp -s "${STAGE}/CODEX_HANDOFF_20260804.md" "${PROJECT}/CODEX_HANDOFF.md"
cmp -s "${STAGE}/RESULTS_MANIFEST_20260804.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cmp -s "${STAGE}/PART1_REBUILD_STATUS_20260804.md" "${PART1_STATUS}"
cmp -s "${STAGE}/part1_malignancy_cnv_manual_figure_review_v1.tsv" "${MANUAL_DST}"

echo PART1_CNV_TERMINAL_RECORDS_INSTALL_PASS
