#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
DEST="${PROJECT}/results/part1_rebuild_v1/09_metabolism_trajectory"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_completion_stage_v1
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

"${PROJECT}/check_memory_gate.sh" part1_slingshot_figure_v2_install
test -x "${PYTHON}"
"${PYTHON}" --version
cp "${STAGE}/repair_part1_slingshot_figure_v2.py" "${DEST}/repair_part1_slingshot_figure_v2.py"
cmp -s "${STAGE}/repair_part1_slingshot_figure_v2.py" "${DEST}/repair_part1_slingshot_figure_v2.py"
"${PYTHON}" -m py_compile "${DEST}/repair_part1_slingshot_figure_v2.py"
"${PYTHON}" "${DEST}/repair_part1_slingshot_figure_v2.py" \
  > "${DEST}/trajectory_v1/logs/repair_part1_slingshot_figure_v2.log" 2>&1
grep -F 'PART1_SLINGSHOT_FIGURE_V2_PASS' \
  "${DEST}/trajectory_v1/logs/repair_part1_slingshot_figure_v2.log"
