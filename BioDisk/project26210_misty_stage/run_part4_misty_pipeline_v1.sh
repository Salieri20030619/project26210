#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
RESULT_ROOT="${PROJECT_ROOT}/results/part4_spatial"
ENV_ACTIVATE="/mnt/BioDisk/wangmingyuan/project26210_envs/activate_part4_advanced_env_v1.sh"
INPUT_DIR="${RESULT_ROOT}/advanced_patch_misty_inputs_v1"
SCRIPT_DIR="${RESULT_ROOT}"

source "${ENV_ACTIVATE}"
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${SCRIPT_DIR}/prepare_part4_misty_inputs_v1.py"

"${PROJECT_ROOT}/check_memory_gate.sh"
Rscript "${SCRIPT_DIR}/run_part4_misty_sample_v1.R" \
  GLMF1 "${INPUT_DIR}" "${RESULT_ROOT}/advanced_patch_misty_glmf1_v1" 20260821

"${PROJECT_ROOT}/check_memory_gate.sh"
Rscript "${SCRIPT_DIR}/run_part4_misty_sample_v1.R" \
  GLMF2 "${INPUT_DIR}" "${RESULT_ROOT}/advanced_patch_misty_glmf2_v1" 20260822

echo "PART4_MISTY_PIPELINE_ALL_PASS"
