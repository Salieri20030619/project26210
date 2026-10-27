#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -ne 6 ]]; then
  echo "usage: $0 ARCHIVE SIGNATURE_TSV OUTPUT_ROOT PYTHON WORKERS EXPECTED_SHA256" >&2
  exit 64
fi

ARCHIVE="$1"
SIGNATURE_TSV="$2"
OUTPUT_ROOT="$3"
PYTHON_BIN="$4"
WORKERS="$5"
EXPECTED_SHA256="$6"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="${OUTPUT_ROOT}/run_v1"
EXTRACT_DIR="${OUTPUT_ROOT}/input_extract_v1"
EXTRACT_PART="${OUTPUT_ROOT}/input_extract_v1.part"

export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export NUMEXPR_NUM_THREADS=2

mkdir -p "${OUTPUT_ROOT}/logs_v1"
exec > >(tee -a "${OUTPUT_ROOT}/logs_v1/background_run_v1.log") 2>&1

on_error() {
  local code=$?
  printf 'status=FAIL\nexit_code=%s\nfailed_at=%s\n' "${code}" "$(date -Is)" > "${OUTPUT_ROOT}/PHASE54_BACKGROUND_FAIL_v1"
  exit "${code}"
}
trap on_error ERR

printf 'status=RUNNING\npid=%s\nstarted_at=%s\nworkers=%s\n' "$$" "$(date -Is)" "${WORKERS}" > "${OUTPUT_ROOT}/PHASE54_BACKGROUND_RUNNING_v1"

if [[ ! -f "${ARCHIVE}" ]]; then
  echo "Archive missing: ${ARCHIVE}" >&2
  exit 2
fi
if [[ ! -f "${SIGNATURE_TSV}" ]]; then
  echo "Signature table missing: ${SIGNATURE_TSV}" >&2
  exit 2
fi
OBSERVED_SHA256="$(sha256sum "${ARCHIVE}" | awk '{print $1}')"
if [[ "${OBSERVED_SHA256}" != "${EXPECTED_SHA256}" ]]; then
  echo "Archive SHA256 mismatch: ${OBSERVED_SHA256}" >&2
  exit 3
fi

if [[ ! -d "${EXTRACT_DIR}" ]]; then
  if [[ -e "${EXTRACT_PART}" ]]; then
    echo "Refusing existing extraction staging directory: ${EXTRACT_PART}" >&2
    exit 4
  fi
  mkdir -p "${EXTRACT_PART}"
  tar -xf "${ARCHIVE}" -C "${EXTRACT_PART}"
  FILTERED_COUNT="$(find "${EXTRACT_PART}" -maxdepth 1 -type f -name 'GSM*_filtered_feature_bc_matrix.h5' | wc -l)"
  RAW_COUNT="$(find "${EXTRACT_PART}" -maxdepth 1 -type f -name 'GSM*_raw_feature_bc_matrix.h5' | wc -l)"
  POSITION_COUNT="$(find "${EXTRACT_PART}" -maxdepth 1 -type f -name 'GSM*_tissue_positions_list.csv.gz' | wc -l)"
  if [[ "${FILTERED_COUNT}" -ne 24 || "${RAW_COUNT}" -ne 24 || "${POSITION_COUNT}" -ne 24 ]]; then
    echo "Extraction validation failed filtered=${FILTERED_COUNT} raw=${RAW_COUNT} positions=${POSITION_COUNT}" >&2
    exit 5
  fi
  mv "${EXTRACT_PART}" "${EXTRACT_DIR}"
fi

mkdir -p "${RUN_DIR}"
"${PYTHON_BIN}" "${SCRIPT_DIR}/01_process_gse175540_spatial_v1.py" \
  --input-dir "${EXTRACT_DIR}" \
  --output-dir "${RUN_DIR}" \
  --signatures "${SIGNATURE_TSV}" \
  --workers "${WORKERS}"

printf 'status=PASS\ncompleted_at=%s\nworkers=%s\narchive_sha256=%s\n' "$(date -Is)" "${WORKERS}" "${OBSERVED_SHA256}" > "${OUTPUT_ROOT}/PHASE54_BACKGROUND_PASS_v1"
echo "PHASE54 BACKGROUND RUN PASS"
