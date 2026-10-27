#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 7 ]]; then
  echo "Usage: $0 STAGE_ROOT RSCRIPT CROSSWALK PREFLIGHT_DIR SOURCE_ROOT GENE_ORDER OUTPUT_ROOT" >&2
  exit 2
fi

STAGE_ROOT="$1"
RSCRIPT="$2"
CROSSWALK="$3"
PREFLIGHT_DIR="$4"
SOURCE_ROOT="$5"
GENE_ORDER="$6"
OUTPUT_ROOT="$7"
SCRIPT="${STAGE_ROOT}/scripts/12_prepare_external_cnv_full_inputs_v1.R"
LOG="${STAGE_ROOT}/logs/full_input_build_v1.log"
PID_FILE="${STAGE_ROOT}/logs/full_input_build_v1.pid"

if [[ -e "${OUTPUT_ROOT}" || -e "${PID_FILE}" ]]; then
  echo "No-overwrite/duplicate-launch gate blocked full input build" >&2
  exit 1
fi

nohup "${RSCRIPT}" --vanilla "${SCRIPT}" \
  --crosswalk "${CROSSWALK}" \
  --preflight-dir "${PREFLIGHT_DIR}" \
  --source-root "${SOURCE_ROOT}" \
  --gene-order "${GENE_ORDER}" \
  --output-root "${OUTPUT_ROOT}" \
  --unit-type all >"${LOG}" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "${pid}" >"${PID_FILE}.part"
mv "${PID_FILE}.part" "${PID_FILE}"
date --iso-8601=seconds >"${STAGE_ROOT}/results/FULL_INPUT_BUILD_LAUNCHED_AT"
echo "EXTERNAL_CNV_FULL_INPUT_BUILD_LAUNCHED pid=${pid} log=${LOG}"
