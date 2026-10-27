#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
STAGE_ROOT="/mnt/BioDisk/wangmingyuan/project26210_tymp_state_stage"
PYTHON_BIN="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
SCRIPT_DIR="${PROJECT_ROOT}/results/part4_spatial"
OUTPUT_DIR="${SCRIPT_DIR}/tymp_spp1_state_v1"
RUN_LOG="$(mktemp /tmp/project26210_tymp_spp1_state_run.XXXXXX.log)"

if [[ -e "${OUTPUT_DIR}" ]]; then
  echo "Refusing to overwrite existing output: ${OUTPUT_DIR}" >&2
  exit 2
fi

install -m 0644 "${STAGE_ROOT}/run_tymp_spp1_state_analysis_v1.py" "${SCRIPT_DIR}/run_tymp_spp1_state_analysis_v1.py"
install -m 0644 "${STAGE_ROOT}/verify_tymp_spp1_state_analysis_v1.py" "${SCRIPT_DIR}/verify_tymp_spp1_state_analysis_v1.py"

"${PROJECT_ROOT}/check_memory_gate.sh"
set +e
"${PYTHON_BIN}" "${SCRIPT_DIR}/run_tymp_spp1_state_analysis_v1.py" 2>&1 | tee "${RUN_LOG}"
analysis_status=${PIPESTATUS[0]}
set -e
if [[ ${analysis_status} -ne 0 ]]; then
  echo "Analysis failed with status ${analysis_status}; temporary log retained at ${RUN_LOG}" >&2
  exit "${analysis_status}"
fi
install -m 0644 "${RUN_LOG}" "${OUTPUT_DIR}/run.log"

"${PROJECT_ROOT}/check_memory_gate.sh"
"${PYTHON_BIN}" "${SCRIPT_DIR}/verify_tymp_spp1_state_analysis_v1.py" "${OUTPUT_DIR}" 2>&1 | tee "${OUTPUT_DIR}/verification_run.log"

echo "TYMP_SPP1_STATE_FORMAL_RUN_ALL_PASS output=${OUTPUT_DIR}"
