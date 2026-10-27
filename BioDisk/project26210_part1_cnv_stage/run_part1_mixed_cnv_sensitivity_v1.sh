#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
SENS_ROOT="${ROOT}/mixed_cnv_sensitivity_v1"
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
RSCRIPT=/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript
THREADS=4

"${PROJECT}/check_memory_gate.sh"
test "$(${PYTHON} -c 'import sys; print(sys.executable)')" = "${PYTHON}"
test -x "${RSCRIPT}"
"${RSCRIPT}" -e 'cat(R.version.string, "\n")'

if test -e "${SENS_ROOT}/RUN_COMPLETE"
then
  echo "Refusing to overwrite completed mixed sensitivity analysis: ${SENS_ROOT}" >&2
  exit 2
fi
if test -d "${SENS_ROOT}" && test -n "$(find "${SENS_ROOT}" -mindepth 1 -print -quit)"
then
  echo "Refusing to reuse non-empty mixed sensitivity directory: ${SENS_ROOT}" >&2
  exit 2
fi

mkdir -p "${SENS_ROOT}"
"${PYTHON}" "${ROOT}/prepare_part1_mixed_cnv_sensitivity_inputs_v1.py" \
  --output-root "${SENS_ROOT}/inputs"

run_infercnv() {
  local unit="$1"
  "${PROJECT}/check_memory_gate.sh"
  nice -n 10 "${RSCRIPT}" "${ROOT}/run_part1_infercnv_unit_v2.R" \
    --input-dir "${SENS_ROOT}/inputs/infercnv/${unit}" \
    --output-dir "${SENS_ROOT}/infercnv/${unit}" \
    --threads "${THREADS}"
}

run_copykat() {
  local unit="$1"
  "${PROJECT}/check_memory_gate.sh"
  nice -n 10 "${RSCRIPT}" "${ROOT}/run_part1_copykat_unit_v2.R" \
    --input-dir "${SENS_ROOT}/inputs/copykat/${unit}" \
    --output-dir "${SENS_ROOT}/copykat/${unit}" \
    --n-cores "${THREADS}"
}

run_infercnv GSE156632__RCC3
run_infercnv GSE156632__RCC5
run_copykat RCC3t
run_copykat RCC5t

date --iso-8601=seconds > "${SENS_ROOT}/RUN_COMPLETE"
echo PART1_MIXED_CNV_SENSITIVITY_V1_RUN_PASS
