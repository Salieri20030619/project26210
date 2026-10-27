#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage

"${PROJECT}/check_memory_gate.sh"
test -f "${PROJECT}/results/part1_rebuild_v1/04_subclones_tymp/input_audit_verification_v1/RUN_COMPLETE"
cp "${STAGE}/CODEX_HANDOFF_part1_subclone_audit_v1.md" "${PROJECT}/CODEX_HANDOFF.md"
cp "${STAGE}/PROJECT_STATUS_part1_subclone_audit_v1.md" "${PROJECT}/PROJECT_STATUS.md"
cp "${STAGE}/RESULTS_MANIFEST_part1_subclone_audit_v1.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cp "${STAGE}/PART1_REBUILD_STATUS_part1_subclone_audit_v1.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
cmp -s "${STAGE}/CODEX_HANDOFF_part1_subclone_audit_v1.md" "${PROJECT}/CODEX_HANDOFF.md"
cmp -s "${STAGE}/PROJECT_STATUS_part1_subclone_audit_v1.md" "${PROJECT}/PROJECT_STATUS.md"
cmp -s "${STAGE}/RESULTS_MANIFEST_part1_subclone_audit_v1.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cmp -s "${STAGE}/PART1_REBUILD_STATUS_part1_subclone_audit_v1.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
echo PART1_SUBCLONE_AUDIT_RECORDS_INSTALL_PASS
