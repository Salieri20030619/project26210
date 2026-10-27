#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 3 ]]; then
  echo "Usage: $0 STAGE_ROOT PROJECT_ROOT RSCRIPT" >&2
  exit 2
fi

STAGE_ROOT="$1"
PROJECT_ROOT="$2"
RSCRIPT="$3"
SCRIPT_ROOT="${STAGE_ROOT}/scripts"
INPUT_ROOT="${STAGE_ROOT}/pilot_inputs_v1"
RESULT_ROOT="${STAGE_ROOT}/results/pilot_v1"
LOG_ROOT="${STAGE_ROOT}/logs"

mkdir -p "${RESULT_ROOT}/infercnv" "${RESULT_ROOT}/copykat" "${LOG_ROOT}"

launch_one() {
  local name="$1"
  local script="$2"
  local input_dir="$3"
  local output_dir="$4"
  shift 4
  local log_file="${LOG_ROOT}/${name}.log"
  local pid_file="${LOG_ROOT}/${name}.pid"
  if [[ -e "${pid_file}" || -e "${output_dir}/RUN_COMPLETE" || -e "${output_dir}/RUN_NOT_ESTIMABLE" ]]; then
    echo "No-overwrite/duplicate-launch gate blocked ${name}" >&2
    exit 1
  fi
  mkdir -p "${output_dir}"
  nohup "${RSCRIPT}" --vanilla "${script}" \
    --input-dir "${input_dir}" \
    --output-dir "${output_dir}" \
    "$@" >"${log_file}" 2>&1 < /dev/null &
  local pid=$!
  printf '%s\n' "${pid}" >"${pid_file}.part"
  mv "${pid_file}.part" "${pid_file}"
  echo "LAUNCHED name=${name} pid=${pid} log=${log_file}"
}

launch_one \
  "pilot_infercnv_GSE202374_RS01985670" \
  "${SCRIPT_ROOT}/run_external_infercnv_unit_v1.R" \
  "${INPUT_ROOT}/infercnv/GSE202374__RS01985670" \
  "${RESULT_ROOT}/infercnv/GSE202374__RS01985670" \
  --threads 8

launch_one \
  "pilot_copykat_GSE202374_P10" \
  "${SCRIPT_ROOT}/run_external_copykat_unit_v1.R" \
  "${INPUT_ROOT}/copykat/GSE202374__P10" \
  "${RESULT_ROOT}/copykat/GSE202374__P10" \
  --n-cores 8

date --iso-8601=seconds >"${RESULT_ROOT}/PILOT_LAUNCHED_AT"
echo "EXTERNAL_CNV_PILOT_BACKGROUND_LAUNCH_PASS project_root=${PROJECT_ROOT}"
