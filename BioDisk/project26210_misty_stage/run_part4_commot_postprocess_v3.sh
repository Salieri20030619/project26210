#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
RESULT_ROOT="${PROJECT_ROOT}/results/part4_spatial"
ENV_ACTIVATE="/mnt/BioDisk/wangmingyuan/project26210_envs/activate_part4_advanced_env_v1.sh"

source "${ENV_ACTIVATE}"
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${RESULT_ROOT}/derive_part4_commot_offdiagonal_v3.py" GLMF1
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${RESULT_ROOT}/derive_part4_commot_offdiagonal_v3.py" GLMF2
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${RESULT_ROOT}/render_part4_commot_total_score_v3.py"
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${RESULT_ROOT}/summarize_part4_commot_cross_sample_v3.py"
echo "PART4_COMMOT_POSTPROCESS_V3_ALL_PASS"
