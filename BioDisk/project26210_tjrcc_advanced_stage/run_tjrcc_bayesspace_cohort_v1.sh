#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
TJ_ROOT="${PROJECT_ROOT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295"
RUN_ROOT="${TJ_ROOT}/advanced_full_v2"
SCRIPT_DIR="${RUN_ROOT}/scripts"
OUTPUT_ROOT="${RUN_ROOT}/bayesspace"
LOG_ROOT="${RUN_ROOT}/logs/bayesspace"
MAX_WORKERS=2

mkdir -p "${OUTPUT_ROOT}" "${LOG_ROOT}"

samples=(H46_T J38_T R114_T R29_T R51_T S15_T X49_T X98_T Y12_T Y27_T Y7_T Z43_T)
inputs=(
  "${TJ_ROOT}/cell2location_remaining_v1/H46_T/H46_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/J38_T/J38_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/R114_T/R114_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/R29_T/R29_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_pilot_v4/high_R51_T/R51_T_cell2location_pilot.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/S15_T/S15_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_pilot_v4/medium_X49_T/X49_T_cell2location_pilot.h5ad"
  "${TJ_ROOT}/cell2location_pilot_v4/low_X98_T/X98_T_cell2location_pilot.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/Y12_T/Y12_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/Y27_T/Y27_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/Y7_T/Y7_T_cell2location.h5ad"
  "${TJ_ROOT}/cell2location_remaining_v1/Z43_T/Z43_T_cell2location.h5ad"
)

run_one() {
  local index="$1"
  local sample="${samples[$index]}"
  local input="${inputs[$index]}"
  local output="${OUTPUT_ROOT}/${sample}"
  local log="${LOG_ROOT}/${sample}.log"
  local seed="$((20260830 + index))"

  if [[ -f "${output}/status.json" ]] && rg -q '"status": "completed"' "${output}/status.json"; then
    echo "[$(date '+%F %T %z')] SKIP completed sample=${sample}" | tee -a "${log}"
    return 0
  fi
  if [[ ! -s "${input}" ]]; then
    echo "Missing input: ${input}" >&2
    return 1
  fi
  mkdir -p "${output}"
  {
    echo "[$(date '+%F %T %z')] START sample=${sample} seed=${seed}"
    "${PROJECT_ROOT}/check_memory_gate.sh"
    "${SCRIPT_DIR}/run_tjrcc_bayesspace_sample_v1.sh" "${sample}" "${input}" "${output}" "${seed}"
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
      else
        if ! wait "${pid}"; then
          fail=1
        fi
      fi
    done
    pids=("${next[@]}")
    if (( ${#pids[@]} >= MAX_WORKERS )); then
      sleep 2
    fi
  done
done

for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    fail=1
  fi
done

if (( fail != 0 )); then
  echo "TJRCC_BAYESSPACE_COHORT_FAILED" >&2
  exit 1
fi
echo "TJRCC_BAYESSPACE_COHORT_ALL_PASS samples=${#samples[@]}"
