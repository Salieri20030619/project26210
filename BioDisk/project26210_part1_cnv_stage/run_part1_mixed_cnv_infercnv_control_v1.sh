#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
INPUT_ROOT="${ROOT}/cnv_unit_inputs_v1/infercnv"
OUT_ROOT="${ROOT}/mixed_cnv_sensitivity_v1/control_no_add_4threads/infercnv"
RSCRIPT=/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript
THREADS=4

"${PROJECT}/check_memory_gate.sh"
test -x "${RSCRIPT}"
"${RSCRIPT}" -e 'cat(R.version.string, "\n")'
if test -d "${OUT_ROOT}" && test -n "$(find "${OUT_ROOT}" -mindepth 1 -print -quit)"
then
  echo "Refusing to overwrite non-empty inferCNV control directory: ${OUT_ROOT}" >&2
  exit 2
fi

run_control() {
  local unit="$1"
  "${PROJECT}/check_memory_gate.sh"
  nice -n 10 "${RSCRIPT}" "${ROOT}/run_part1_infercnv_unit_v2.R" \
    --input-dir "${INPUT_ROOT}/${unit}" \
    --output-dir "${OUT_ROOT}/${unit}" \
    --threads "${THREADS}"
}

run_control GSE156632__RCC3
run_control GSE156632__RCC5
date --iso-8601=seconds > "${ROOT}/mixed_cnv_sensitivity_v1/control_no_add_4threads/RUN_COMPLETE"
echo PART1_MIXED_CNV_INFERCNV_CONTROL_V1_PASS
