#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage

"${PROJECT}/check_memory_gate.sh"
awk -F '\t' 'NF != 18 {print "bad_manifest_row", NR, NF; exit 1}' "${STAGE}/RESULTS_MANIFEST.tsv"
rg -q "最近完成：Part 1 新版cNMF全rank、TYMP证据与rank选择" "${STAGE}/CODEX_HANDOFF.md"
rg -q "selected_with_study_composition_warning" "${STAGE}/PROJECT_STATUS.md"
rg -q "新版cNMF输入、全rank与TYMP rank选择终态" "${STAGE}/PART1_REBUILD_STATUS.md"
rg -q "p1_cnmf_rank_selection_verification_v1_20260804" "${STAGE}/RESULTS_MANIFEST.tsv"

cp "${STAGE}/CODEX_HANDOFF.md" "${PROJECT}/CODEX_HANDOFF.md"
cp "${STAGE}/PROJECT_STATUS.md" "${PROJECT}/PROJECT_STATUS.md"
cp "${STAGE}/RESULTS_MANIFEST.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cp "${STAGE}/PART1_REBUILD_STATUS.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"

cmp -s "${STAGE}/CODEX_HANDOFF.md" "${PROJECT}/CODEX_HANDOFF.md"
cmp -s "${STAGE}/PROJECT_STATUS.md" "${PROJECT}/PROJECT_STATUS.md"
cmp -s "${STAGE}/RESULTS_MANIFEST.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cmp -s "${STAGE}/PART1_REBUILD_STATUS.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
echo PART1_CNMF_TERMINAL_RECORDS_INSTALL_PASS
