#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 5 ]]; then
  echo "usage: $0 REFERENCE_H5AD SPATIAL_SECTIONS OUTPUT_ROOT PYTHON WORKERS" >&2
  exit 64
fi

REFERENCE_H5AD="$1"
SPATIAL_SECTIONS="$2"
OUTPUT_ROOT="$3"
PYTHON_BIN="$4"
WORKERS="$5"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
mkdir -p "${OUTPUT_ROOT}/logs_v1"
exec > >(tee -a "${OUTPUT_ROOT}/logs_v1/background_cell2location_v1.log") 2>&1

on_error() {
  local code=$?
  printf 'status=FAIL\nexit_code=%s\nfailed_at=%s\n' "${code}" "$(date -Is)" > "${OUTPUT_ROOT}/PHASE55_BACKGROUND_FAIL_v1"
  exit "${code}"
}
trap on_error ERR

if [[ ! -f "${REFERENCE_H5AD}" ]]; then
  echo "Reference H5AD missing: ${REFERENCE_H5AD}" >&2
  exit 2
fi
if [[ "$(find "${SPATIAL_SECTIONS}" -mindepth 1 -maxdepth 1 -type d -name 'GSM*' | wc -l)" -ne 24 ]]; then
  echo "Spatial section directory count is not 24" >&2
  exit 3
fi

printf 'status=RUNNING\npid=%s\nstarted_at=%s\nworkers=%s\n' "$$" "$(date -Is)" "${WORKERS}" > "${OUTPUT_ROOT}/PHASE55_BACKGROUND_RUNNING_v1"
"${PYTHON_BIN}" "${SCRIPT_DIR}/02_run_cell2location_v1.py" \
  --reference "${REFERENCE_H5AD}" \
  --spatial-sections "${SPATIAL_SECTIONS}" \
  --output-root "${OUTPUT_ROOT}/run_v1" \
  --workers "${WORKERS}"
printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "${OUTPUT_ROOT}/PHASE55_BACKGROUND_PASS_v1"

