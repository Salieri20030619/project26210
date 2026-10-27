#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
RESULT_ROOT="${PROJECT_ROOT}/results/part4_spatial"
ENV_ACTIVATE="/mnt/BioDisk/wangmingyuan/project26210_envs/activate_part4_advanced_env_v1.sh"

source "${ENV_ACTIVATE}"
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${RESULT_ROOT}/run_part4_commot_candidate_v2.py" GLMF1
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${RESULT_ROOT}/run_part4_commot_candidate_v2.py" GLMF2
echo "PART4_COMMOT_PIPELINE_V2_ALL_PASS"
