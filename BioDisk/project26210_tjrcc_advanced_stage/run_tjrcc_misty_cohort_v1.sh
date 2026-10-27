#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
TJ_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295"
RUN_ROOT="${TJ_ROOT}/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
INPUT_DIR="${RUN_ROOT}/misty/inputs"
OUTPUT_ROOT="${RUN_ROOT}/misty/samples"
LOG_ROOT="${RUN_ROOT}/logs/misty"
MAX_WORKERS=2

source /mnt/BioDisk/wangmingyuan/project26210_envs/activate_part4_advanced_env_v1.sh
expected_rscript="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript"
[[ "$(command -v Rscript)" == "${expected_rscript}" ]] || { echo "Unexpected Rscript: $(command -v Rscript)" >&2; exit 4; }
mkdir -p "${OUTPUT_ROOT}" "${LOG_ROOT}"
"${PROJECT_ROOT}/check_memory_gate.sh"
python "${SCRIPT_DIR}/prepare_tjrcc_misty_inputs_v1.py" >"${LOG_ROOT}/prepare_inputs_v1.log" 2>&1

samples=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)

run_one() {
  local index="$1"
  local sample="${samples[$index]}"
  local output="${OUTPUT_ROOT}/${sample}"
  local log="${LOG_ROOT}/${sample}.log"
  local seed="$((20260920 + index))"
  if [[ -f "${output}/status.json" ]] && rg -q '"status": "complete"' "${output}/status.json"; then
    echo "[$(date '+%F %T %z')] SKIP completed sample=${sample}" | tee -a "${log}"
    return 0
  fi
  mkdir -p "${output}"
  {
    echo "[$(date '+%F %T %z')] START sample=${sample} seed=${seed}"
    "${PROJECT_ROOT}/check_memory_gate.sh"
    Rscript --vanilla "${SCRIPT_DIR}/run_tjrcc_misty_sample_v1.R" "${sample}" "${INPUT_DIR}" "${output}" "${seed}"
    echo "[$(date '+%F %T %z')] PASS sample=${sample}"
  } >"${log}" 2>&1
}

declare -a pids=()
fail=0
for index in "${!samples[@]}"; do
  run_one "${index}" &
  pids+=("$!")
  while (( ${#pids[@]} >= MAX_WORKERS )); do
    next=()
    for pid in "${pids[@]}"; do
      if kill -0 "${pid}" 2>/dev/null; then
        next+=("${pid}")
      elif ! wait "${pid}"; then
        fail=1
      fi
    done
    pids=("${next[@]}")
    (( ${#pids[@]} < MAX_WORKERS )) || sleep 2
  done
done
for pid in "${pids[@]}"; do
  wait "${pid}" || fail=1
done
(( fail == 0 )) || { echo "TJRCC_MISTY_COHORT_FAILED" >&2; exit 1; }
echo "TJRCC_MISTY_COHORT_ALL_PASS samples=${#samples[@]}"
