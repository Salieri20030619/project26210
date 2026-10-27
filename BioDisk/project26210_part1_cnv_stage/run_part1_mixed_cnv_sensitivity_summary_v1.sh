#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
OUT="${ROOT}/mixed_cnv_sensitivity_v1/summary_v1"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

"${PROJECT}/check_memory_gate.sh"
test "$(${PYTHON} -c 'import sys; print(sys.executable)')" = "${PYTHON}"
test -f "${ROOT}/mixed_cnv_sensitivity_v1/RUN_COMPLETE"
if test -d "${OUT}" && test -n "$(find "${OUT}" -mindepth 1 -print -quit)"
then
  echo "Refusing to overwrite non-empty summary directory: ${OUT}" >&2
  exit 2
fi
"${PYTHON}" "${ROOT}/summarize_part1_mixed_cnv_sensitivity_v1.py"
