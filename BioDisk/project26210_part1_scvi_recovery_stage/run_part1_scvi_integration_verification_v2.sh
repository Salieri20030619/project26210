#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
SCRIPT="${PROJECT}/results/part1_rebuild_v1/02_integration_annotation/verify_part1_scvi_integration_v2.py"
LOG="${PROJECT}/results/part1_rebuild_v1/logs/scvi_integration_verification_v2.log"

if [[ ! -x "${PYTHON}" ]]; then
  echo "missing executable: ${PYTHON}" >&2
  exit 1
fi
if [[ ! -f "${SCRIPT}" ]]; then
  echo "missing verifier: ${SCRIPT}" >&2
  exit 1
fi

"${PROJECT}/check_memory_gate.sh" | tee "${LOG}"
"${PYTHON}" --version 2>&1 | tee -a "${LOG}"
"${PYTHON}" "${SCRIPT}" 2>&1 | tee -a "${LOG}"
