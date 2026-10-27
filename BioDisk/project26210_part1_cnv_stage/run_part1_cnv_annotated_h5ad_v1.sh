#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

"${PROJECT}/check_memory_gate.sh"
test "$(${PYTHON} -c 'import sys; print(sys.executable)')" = "${PYTHON}"
test -f "${ROOT}/malignancy_cnv_summary_v1/RUN_COMPLETE"
test -f "${ROOT}/mixed_cnv_sensitivity_v1/summary_v1/RUN_COMPLETE"
test -f "${ROOT}/mixed_cnv_sensitivity_v1/summary_v2_with_control/RUN_COMPLETE"
test -f "${ROOT}/mixed_cnv_sensitivity_verification_v1/RUN_COMPLETE"
nice -n 10 "${PYTHON}" "${ROOT}/build_part1_cnv_annotated_h5ad_v1.py"
