#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
STAGE="/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage"
OUT="${1:-${STAGE}/cnv_pilot_v2}"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
RSCRIPT="/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript"
LOG="${OUT}/cnv_pilot_v2.log"

mkdir -p "${OUT}"
{
  "${PROJECT}/check_memory_gate.sh"
  "${PYTHON}" "${STAGE}/prepare_part1_cnv_unit_inputs_v2.py" \
    --output-dir "${OUT}/inputs" \
    --unit-filter "GSE156632__RCC1,RCC1t"

  "${PROJECT}/check_memory_gate.sh"
  "${RSCRIPT}" --vanilla "${STAGE}/run_part1_infercnv_unit_v1.R" \
    --input-dir "${OUT}/inputs/infercnv/GSE156632__RCC1" \
    --output-dir "${OUT}/infercnv/GSE156632__RCC1" \
    --threads 8

  "${PROJECT}/check_memory_gate.sh"
  "${RSCRIPT}" --vanilla "${STAGE}/run_part1_copykat_unit_v1.R" \
    --input-dir "${OUT}/inputs/copykat/RCC1t" \
    --output-dir "${OUT}/copykat/RCC1t" \
    --n-cores 8

  test -f "${OUT}/inputs/RUN_COMPLETE"
  test -f "${OUT}/infercnv/GSE156632__RCC1/RUN_COMPLETE"
  test -f "${OUT}/copykat/RCC1t/RUN_COMPLETE"
  echo "PART1_CNV_PILOT_V2_ALL_PASS"
} 2>&1 | tee "${LOG}"
