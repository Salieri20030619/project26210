#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 5 ]]; then
  echo "Usage: $0 STAGE_ROOT INPUT_ROOT VALIDATION_DIR RSCRIPT MEMORY_GATE" >&2
  exit 2
fi
STAGE_ROOT="$1"
INPUT_ROOT="$2"
VALIDATION_DIR="$3"
RSCRIPT="$4"
MEMORY_GATE="$5"
RESULT_ROOT="${STAGE_ROOT}/results/full_v2"
LOG_ROOT="${STAGE_ROOT}/logs/full_v2"
STATUS_ROOT="${RESULT_ROOT}/status"

if [[ ! -f "${INPUT_ROOT}/RUN_COMPLETE" || ! -f "${VALIDATION_DIR}/RUN_COMPLETE" ]]; then
  echo "Validated-v2 input gate failed" >&2
  exit 1
fi
if [[ -e "${RESULT_ROOT}" || -e "${LOG_ROOT}" ]]; then
  echo "No-overwrite gate: full_v2 result or log root already exists" >&2
  exit 1
fi
if [[ ! -x "${RSCRIPT}" || ! -x "${MEMORY_GATE}" ]]; then
  echo "Rscript or memory gate is not executable" >&2
  exit 1
fi
"${MEMORY_GATE}"
mkdir -p "${STATUS_ROOT}" "${LOG_ROOT}"

nohup "${STAGE_ROOT}/scripts/21_run_external_infercnv_cohort_v2.sh" \
  "${STAGE_ROOT}" "${INPUT_ROOT}" "${RSCRIPT}" "${MEMORY_GATE}" \
  >"${LOG_ROOT}/infercnv_cohort_v2.log" 2>&1 < /dev/null &
infer_pid=$!
printf '%s\n' "${infer_pid}" >"${STATUS_ROOT}/infercnv_cohort_v2.pid.part"
mv "${STATUS_ROOT}/infercnv_cohort_v2.pid.part" "${STATUS_ROOT}/infercnv_cohort_v2.pid"

nohup "${STAGE_ROOT}/scripts/22_run_external_copykat_cohort_v2.sh" \
  "${STAGE_ROOT}" "${INPUT_ROOT}" "${RSCRIPT}" "${MEMORY_GATE}" \
  >"${LOG_ROOT}/copykat_cohort_v2.log" 2>&1 < /dev/null &
copykat_pid=$!
printf '%s\n' "${copykat_pid}" >"${STATUS_ROOT}/copykat_cohort_v2.pid.part"
mv "${STATUS_ROOT}/copykat_cohort_v2.pid.part" "${STATUS_ROOT}/copykat_cohort_v2.pid"

printf '%s\tinfercnv_pid=%s\tcopykat_pid=%s\n' "$(date --iso-8601=seconds)" "${infer_pid}" "${copykat_pid}" \
  >"${STATUS_ROOT}/COHORTS_LAUNCHED.part"
mv "${STATUS_ROOT}/COHORTS_LAUNCHED.part" "${STATUS_ROOT}/COHORTS_LAUNCHED"
echo "EXTERNAL_CNV_COHORTS_V2_LAUNCHED infercnv_pid=${infer_pid} copykat_pid=${copykat_pid}"
