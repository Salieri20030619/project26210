#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage

"${PROJECT}/check_memory_gate.sh"
awk -F '\t' 'NF != 18 {print "invalid manifest fields at line " NR ": " NF > "/dev/stderr"; bad=1} END {exit bad}' "${STAGE}/RESULTS_MANIFEST_v2_running.tsv"
rg -q 'CNV v1 部分失败与 v2 统一全量重跑' "${STAGE}/PART1_REBUILD_STATUS_v2_running.md"
rg -q 'v2 门控试跑通过并已进入统一全量重跑' "${STAGE}/PROJECT_STATUS_v2_running.md"
rg -q 'CNV v1 部分失败与 v2 统一全量重跑启动' "${STAGE}/CODEX_HANDOFF_v2_running.md"
rg -q '^p1_cnv_v2_launcher_20260803\t' "${STAGE}/RESULTS_MANIFEST_v2_running.tsv"

cp "${STAGE}/PART1_REBUILD_STATUS_v2_running.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
cp "${STAGE}/PROJECT_STATUS_v2_running.md" "${PROJECT}/PROJECT_STATUS.md"
cp "${STAGE}/CODEX_HANDOFF_v2_running.md" "${PROJECT}/CODEX_HANDOFF.md"
cp "${STAGE}/RESULTS_MANIFEST_v2_running.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"

cmp -s "${STAGE}/PART1_REBUILD_STATUS_v2_running.md" "${PROJECT}/results/part1_rebuild_v1/PART1_REBUILD_STATUS.md"
cmp -s "${STAGE}/PROJECT_STATUS_v2_running.md" "${PROJECT}/PROJECT_STATUS.md"
cmp -s "${STAGE}/CODEX_HANDOFF_v2_running.md" "${PROJECT}/CODEX_HANDOFF.md"
cmp -s "${STAGE}/RESULTS_MANIFEST_v2_running.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
echo PART1_CNV_V2_RUNNING_RECORDS_INSTALL_PASS
