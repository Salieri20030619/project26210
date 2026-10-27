#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_completion_stage_v1
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
DEST="${ROOT}/09_metabolism_trajectory/finalize_part1_existing_manual_figure_review_v1.py"

"${PROJECT}/check_memory_gate.sh" part1_existing_manual_figure_review_v1
test -x "${PYTHON}"
cp "${STAGE}/finalize_part1_existing_manual_figure_review_v1.py" "${DEST}"
cmp -s "${STAGE}/finalize_part1_existing_manual_figure_review_v1.py" "${DEST}"
"${PYTHON}" -m py_compile "${DEST}"
"${PYTHON}" "${DEST}" | tee "${ROOT}/09_metabolism_trajectory/trajectory_v1/logs/finalize_manual_figure_review_v1.log"
grep -F 'PART1_EXISTING_MANUAL_FIGURE_REVIEW_V1_PASS' \
  "${ROOT}/09_metabolism_trajectory/trajectory_v1/logs/finalize_manual_figure_review_v1.log"
