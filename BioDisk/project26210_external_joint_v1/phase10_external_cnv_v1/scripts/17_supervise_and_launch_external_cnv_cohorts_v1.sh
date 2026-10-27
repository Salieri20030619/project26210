#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 7 ]]; then
  echo "Usage: $0 STAGE_ROOT INPUT_ROOT PREFLIGHT_DIR BUILDER_PID PYTHON RSCRIPT MEMORY_GATE" >&2
  exit 2
fi
STAGE_ROOT="$1"
INPUT_ROOT="$2"
PREFLIGHT_DIR="$3"
BUILDER_PID="$4"
PYTHON="$5"
RSCRIPT="$6"
MEMORY_GATE="$7"
LOG="${STAGE_ROOT}/logs/full_cohort_supervisor_v1.log"
VALIDATION_DIR="${STAGE_ROOT}/results/full_input_validation_v1"
STATUS_ROOT="${STAGE_ROOT}/results/full_v1/status"
mkdir -p "${STATUS_ROOT}"

while kill -0 "${BUILDER_PID}" 2>/dev/null; do
  printf '%s waiting_for_builder pid=%s\n' "$(date --iso-8601=seconds)" "${BUILDER_PID}" >>"${LOG}"
  sleep 30
done

if [[ ! -f "${INPUT_ROOT}/RUN_COMPLETE" ]]; then
  printf '%s builder_exited_without_RUN_COMPLETE\n' "$(date --iso-8601=seconds)" >>"${LOG}"
  date --iso-8601=seconds >"${STATUS_ROOT}/SUPERVISOR_BLOCKED_INPUT_BUILD"
  exit 1
fi

"${PYTHON}" "${STAGE_ROOT}/scripts/14_validate_external_cnv_full_inputs_v1.py" \
  --input-root "${INPUT_ROOT}" --preflight-dir "${PREFLIGHT_DIR}" --out-dir "${VALIDATION_DIR}" \
  >>"${LOG}" 2>&1
if [[ ! -f "${VALIDATION_DIR}/RUN_COMPLETE" ]]; then
  date --iso-8601=seconds >"${STATUS_ROOT}/SUPERVISOR_BLOCKED_VALIDATION"
  exit 1
fi

nohup "${STAGE_ROOT}/scripts/15_run_external_infercnv_cohort_v1.sh" \
  "${STAGE_ROOT}" "${INPUT_ROOT}" "${RSCRIPT}" "${MEMORY_GATE}" \
  >"${STAGE_ROOT}/logs/full_infercnv_cohort_v1.log" 2>&1 < /dev/null &
infer_pid=$!
printf '%s\n' "${infer_pid}" >"${STATUS_ROOT}/infercnv_cohort_v1.pid.part"
mv "${STATUS_ROOT}/infercnv_cohort_v1.pid.part" "${STATUS_ROOT}/infercnv_cohort_v1.pid"

nohup "${STAGE_ROOT}/scripts/16_run_external_copykat_cohort_v1.sh" \
  "${STAGE_ROOT}" "${INPUT_ROOT}" "${RSCRIPT}" "${MEMORY_GATE}" \
  >"${STAGE_ROOT}/logs/full_copykat_cohort_v1.log" 2>&1 < /dev/null &
copykat_pid=$!
printf '%s\n' "${copykat_pid}" >"${STATUS_ROOT}/copykat_cohort_v1.pid.part"
mv "${STATUS_ROOT}/copykat_cohort_v1.pid.part" "${STATUS_ROOT}/copykat_cohort_v1.pid"

printf '%s\tinfercnv_pid=%s\tcopykat_pid=%s\n' "$(date --iso-8601=seconds)" "${infer_pid}" "${copykat_pid}" >>"${LOG}"
date --iso-8601=seconds >"${STATUS_ROOT}/COHORTS_LAUNCHED"
echo "EXTERNAL_CNV_COHORTS_LAUNCHED infercnv_pid=${infer_pid} copykat_pid=${copykat_pid}"
