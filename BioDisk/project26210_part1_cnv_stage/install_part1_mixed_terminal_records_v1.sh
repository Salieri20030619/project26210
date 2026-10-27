#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage
VERIFY="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv/mixed_cnv_sensitivity_verification_v1"

"${PROJECT}/check_memory_gate.sh"
test -f "${VERIFY}/RUN_COMPLETE"
test -f "${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_verification_v1/RUN_COMPLETE"
cp "${STAGE}/CODEX_HANDOFF_part1_mixed_terminal_v1.md" "${PROJECT}/CODEX_HANDOFF.md"
cp "${STAGE}/PROJECT_STATUS_part1_mixed_terminal_v1.md" "${PROJECT}/PROJECT_STATUS.md"
cp "${STAGE}/RESULTS_MANIFEST_part1_mixed_terminal_v1.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cp "${STAGE}/PART1_REBUILD_STATUS_part1_mixed_terminal_v1.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
cp "${STAGE}/part1_mixed_cnv_manual_figure_review_v1.tsv" "${VERIFY}/part1_mixed_cnv_manual_figure_review_v1.tsv"
cmp -s "${STAGE}/CODEX_HANDOFF_part1_mixed_terminal_v1.md" "${PROJECT}/CODEX_HANDOFF.md"
cmp -s "${STAGE}/PROJECT_STATUS_part1_mixed_terminal_v1.md" "${PROJECT}/PROJECT_STATUS.md"
cmp -s "${STAGE}/RESULTS_MANIFEST_part1_mixed_terminal_v1.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cmp -s "${STAGE}/PART1_REBUILD_STATUS_part1_mixed_terminal_v1.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
cmp -s "${STAGE}/part1_mixed_cnv_manual_figure_review_v1.tsv" "${VERIFY}/part1_mixed_cnv_manual_figure_review_v1.tsv"
echo PART1_MIXED_TERMINAL_RECORDS_INSTALL_PASS
