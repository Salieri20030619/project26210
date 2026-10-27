#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
STAGE="/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage"
PILOT="${1:-${STAGE}/cnv_pilot_v2}"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
LOG="${PILOT}/verification_v2.log"

{
  "${PROJECT}/check_memory_gate.sh"
  "${PYTHON}" "${STAGE}/verify_part1_cnv_pilot_v2_verification_v2.py" --pilot-dir "${PILOT}"
} 2>&1 | tee "${LOG}"
