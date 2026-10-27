#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
LOG="${PROJECT}/results/part1_rebuild_v1/logs/cnv_v2_pilot_then_full_v1.log"

{
  "${PROJECT}/check_memory_gate.sh"
  bash "${ROOT}/run_part1_cnv_v2_pilot.sh"
  "${PROJECT}/check_memory_gate.sh"
  "${PYTHON}" "${ROOT}/verify_part1_cnv_v2_pilot.py"
  bash "${ROOT}/run_part1_malignancy_cnv_full_v2.sh"
} 2>&1 | tee -a "${LOG}"
