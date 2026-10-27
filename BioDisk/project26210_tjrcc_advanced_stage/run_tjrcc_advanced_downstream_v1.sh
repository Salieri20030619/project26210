#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
RUN_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
LOG_DIR="${RUN_ROOT}/logs"
samples=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)

mkdir -p "${LOG_DIR}"
echo "[$(date '+%F %T %z')] WAIT_FOR_BAYESSPACE_START"
while true; do
  complete=0
  for sample in "${samples[@]}"; do
    status="${RUN_ROOT}/bayesspace/${sample}/status.json"
    if [[ -f "${status}" ]] && rg -q '"status": "completed"' "${status}"; then
      complete=$((complete + 1))
    elif [[ -f "${status}" ]] && rg -q '"status": "failed"' "${status}"; then
      echo "BayesSpace failed for ${sample}; refusing to start downstream stages" >&2
      exit 1
    fi
  done
  echo "[$(date '+%F %T %z')] BAYESSPACE_PROGRESS ${complete}/${#samples[@]}"
  (( complete == ${#samples[@]} )) && break
  sleep 30
done

"${PROJECT_ROOT}/check_memory_gate.sh"
"${SCRIPT_DIR}/run_tjrcc_spatial_programs_cohort_v1.sh" |& tee "${LOG_DIR}/spatial_programs_cohort_v1.log"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${SCRIPT_DIR}/run_tjrcc_misty_cohort_v1.sh" |& tee "${LOG_DIR}/misty_cohort_v1.log"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${SCRIPT_DIR}/run_tjrcc_commot_cohort_v1.sh" |& tee "${LOG_DIR}/commot_cohort_v1.log"
echo "[$(date '+%F %T %z')] TJRCC_ADVANCED_SAMPLE_STAGES_ALL_PASS"
