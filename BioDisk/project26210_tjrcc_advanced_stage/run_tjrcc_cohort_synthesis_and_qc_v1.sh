#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
RUN_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
LOG_DIR="${RUN_ROOT}/logs"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
samples=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)
stages=(bayesspace spatial_programs misty commot)

mkdir -p "${LOG_DIR}"
echo "[$(date '+%F %T %z')] WAIT_FOR_ALL_SAMPLE_STAGES_START"
while true; do
  all_complete=1
  progress=()
  for stage in "${stages[@]}"; do
    complete=0
    for sample in "${samples[@]}"; do
      if [[ "${stage}" == "misty" || "${stage}" == "commot" ]]; then
        status="${RUN_ROOT}/${stage}/samples/${sample}/status.json"
      else
        status="${RUN_ROOT}/${stage}/${sample}/status.json"
      fi
      if [[ -f "${status}" ]] && rg -q '"status": "completed"|"status": "complete"' "${status}"; then
        complete=$((complete + 1))
      elif [[ -f "${status}" ]] && rg -q '"status": "failed"' "${status}"; then
        echo "Sample stage failed: stage=${stage} sample=${sample}; refusing cohort synthesis" >&2
        exit 1
      fi
    done
    progress+=("${stage}=${complete}/${#samples[@]}")
    if (( complete != ${#samples[@]} )); then
      all_complete=0
    fi
  done
  echo "[$(date '+%F %T %z')] SAMPLE_STAGE_PROGRESS ${progress[*]}"
  (( all_complete == 1 )) && break
  sleep 60
done

"${PROJECT_ROOT}/check_memory_gate.sh"
"${PYTHON}" "${SCRIPT_DIR}/summarize_tjrcc_advanced_cohort_v1.py" \
  "${RUN_ROOT}/cohort_synthesis_v1" |& tee "${LOG_DIR}/cohort_synthesis_v1.log"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${PYTHON}" "${SCRIPT_DIR}/verify_tjrcc_advanced_full_v1.py" \
  "${RUN_ROOT}/verification_v1" |& tee "${LOG_DIR}/verification_v1.log"
echo "[$(date '+%F %T %z')] TJRCC_ADVANCED_COHORT_AND_QC_ALL_PASS"
