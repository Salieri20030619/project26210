#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
RUN_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
LOG_DIR="${RUN_ROOT}/logs"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"

"${PROJECT_ROOT}/check_memory_gate.sh"
"${SCRIPT_DIR}/run_tjrcc_misty_cohort_v1.sh" |& tee "${LOG_DIR}/misty_deldir_retry_cohort_v1.log"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${SCRIPT_DIR}/run_tjrcc_commot_cohort_v1.sh" |& tee "${LOG_DIR}/commot_cohort_after_misty_retry_v1.log"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${PYTHON}" "${SCRIPT_DIR}/summarize_tjrcc_advanced_cohort_v1.py" \
  "${RUN_ROOT}/cohort_synthesis_v1" |& tee "${LOG_DIR}/cohort_synthesis_v1.log"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${PYTHON}" "${SCRIPT_DIR}/verify_tjrcc_advanced_full_v1.py" \
  "${RUN_ROOT}/verification_v1" |& tee "${LOG_DIR}/verification_v1.log"
echo "[$(date '+%F %T %z')] TJRCC_ADVANCED_DELDIR_RECOVERY_AND_QC_ALL_PASS"
