#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

"${PROJECT}/check_memory_gate.sh"
test "$(${PYTHON} -c 'import sys; print(sys.executable)')" = "${PYTHON}"
test -f "${ROOT}/03_malignancy_cnv/cnv_annotated_h5ad_verification_v1/RUN_COMPLETE"
nice -n 10 "${PYTHON}" "${ROOT}/04_subclones_tymp/audit_part1_subclones_cnmf_inputs_v1.py"
