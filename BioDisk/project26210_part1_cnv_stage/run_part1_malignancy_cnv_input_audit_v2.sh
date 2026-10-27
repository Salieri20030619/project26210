#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
STAGE="/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage"
OUT="${1:-${STAGE}/malignancy_cnv_input_audit_v2}"
LOG="${2:-${STAGE}/malignancy_cnv_input_audit_v2.log}"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
RSCRIPT="/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript"

mkdir -p "${OUT}" "$(dirname "${LOG}")"
{
  "${PROJECT}/check_memory_gate.sh"
  "${PYTHON}" "${STAGE}/audit_part1_malignancy_cnv_inputs_v2.py" \
    --output-dir "${OUT}" \
    --rscript "${RSCRIPT}"
} 2>&1 | tee "${LOG}"
