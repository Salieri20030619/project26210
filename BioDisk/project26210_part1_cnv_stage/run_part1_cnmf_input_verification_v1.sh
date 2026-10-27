#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

"${PROJECT}/check_memory_gate.sh"
test "$(${PYTHON} -c 'import sys; print(sys.executable)')" = "${PYTHON}"
test -f "${ROOT}/05_cnmf/input_v1/RUN_COMPLETE"
nice -n 10 "${PYTHON}" "${ROOT}/05_cnmf/verify_part1_cnmf_input_v1.py"
