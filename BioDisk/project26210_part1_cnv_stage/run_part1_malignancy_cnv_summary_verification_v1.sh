#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
LOG="${PROJECT}/results/part1_rebuild_v1/logs/malignancy_cnv_summary_verification_v1.log"

{
  echo -e "started_at\t$(date --iso-8601=seconds)"
  "${PROJECT}/check_memory_gate.sh"
  "${PYTHON}" "${ROOT}/verify_part1_malignancy_cnv_summary_v1.py"
  echo "PART1_MALIGNANCY_CNV_SUMMARY_VERIFICATION_PIPELINE_ALL_PASS"
} 2>&1 | tee -a "${LOG}"
