#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="/mnt/BioDisk/ludingyu/project26210"
if [[ "$#" -ne 4 ]]; then
  echo "Usage: $0 SAMPLE INPUT_H5AD OUTPUT_DIR SEED" >&2
  exit 2
fi
SAMPLE="$1"
INPUT_H5AD="$2"
OUTPUT_DIR="$3"
SEED="$4"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT_BIN="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript"
PYTHON_BIN="/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python"

export R_LIBS_USER="/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_R_lib_v1"
export R_PROFILE_USER="/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_Rprofile_v1"
export POT_BACKEND_DISABLE_JAX=1
export MPLCONFIGDIR="/tmp/project26210_tjrcc_bayesspace_${SAMPLE}_mpl"
export NUMBA_CACHE_DIR="/tmp/project26210_tjrcc_bayesspace_${SAMPLE}_numba"

mkdir -p "${OUTPUT_DIR}" "${MPLCONFIGDIR}" "${NUMBA_CACHE_DIR}"

"${PROJECT_ROOT}/check_memory_gate.sh"
"${RSCRIPT_BIN}" --vanilla "${SCRIPT_DIR}/run_tjrcc_bayesspace_sample_v1.R" "${SAMPLE}" "${INPUT_H5AD}" "${OUTPUT_DIR}" "${SEED}"
"${PROJECT_ROOT}/check_memory_gate.sh"
"${PYTHON_BIN}" "${SCRIPT_DIR}/postprocess_tjrcc_bayesspace_boundary_v1.py" "${SAMPLE}" "${INPUT_H5AD}" "${OUTPUT_DIR}"

echo "TJRCC_BAYESSPACE_SAMPLE_PIPELINE_ALL_PASS sample=${SAMPLE}"
