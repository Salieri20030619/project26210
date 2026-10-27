#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
RSCRIPT=/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript
PILOT_ROOT="${ROOT}/cnv_v2_pilot"
LOG="${PROJECT}/results/part1_rebuild_v1/logs/cnv_v2_pilot.log"

mkdir -p "${PILOT_ROOT}/infercnv" "${PILOT_ROOT}/copykat" "$(dirname "${LOG}")"

run_infercnv_pilots() {
  local unit input_dir output_dir
  for unit in GSE156632__RCC1 GSE178481__PR5 GSE242299__P3
  do
    input_dir="${ROOT}/cnv_unit_inputs_v1/infercnv/${unit}"
    output_dir="${PILOT_ROOT}/infercnv/${unit}"
    if test -f "${output_dir}/RUN_COMPLETE"
    then
      echo "V2_PILOT_INFERCNV_SKIP unit=${unit}"
      continue
    fi
    mkdir -p "${output_dir}"
    echo "V2_PILOT_INFERCNV_START unit=${unit} time=$(date --iso-8601=seconds)"
    "${PROJECT}/check_memory_gate.sh"
    "${RSCRIPT}" --vanilla "${ROOT}/run_part1_infercnv_unit_v2.R" \
      --input-dir "${input_dir}" --output-dir "${output_dir}" --threads 8
    echo "V2_PILOT_INFERCNV_DONE unit=${unit} time=$(date --iso-8601=seconds)"
  done
}

run_copykat_pilots() {
  local unit input_dir output_dir
  for unit in RCC1t PR5_PTumor1 T3_1
  do
    input_dir="${ROOT}/cnv_unit_inputs_v1/copykat/${unit}"
    output_dir="${PILOT_ROOT}/copykat/${unit}"
    if test -f "${output_dir}/RUN_COMPLETE"
    then
      echo "V2_PILOT_COPYKAT_SKIP unit=${unit}"
      continue
    fi
    mkdir -p "${output_dir}"
    echo "V2_PILOT_COPYKAT_START unit=${unit} time=$(date --iso-8601=seconds)"
    "${PROJECT}/check_memory_gate.sh"
    "${RSCRIPT}" --vanilla "${ROOT}/run_part1_copykat_unit_v2.R" \
      --input-dir "${input_dir}" --output-dir "${output_dir}" --n-cores 8
    echo "V2_PILOT_COPYKAT_DONE unit=${unit} time=$(date --iso-8601=seconds)"
  done
}

{
  echo "V2_PILOT_START time=$(date --iso-8601=seconds)"
  set +e
  run_infercnv_pilots &
  infer_pid=$!
  run_copykat_pilots &
  copy_pid=$!
  wait "${infer_pid}"
  infer_status=$?
  wait "${copy_pid}"
  copy_status=$?
  set -e
  if test "${infer_status}" -ne 0 || test "${copy_status}" -ne 0
  then
    echo "PART1_CNV_V2_PILOT_FAILED infercnv=${infer_status} copykat=${copy_status}"
    exit 1
  fi
  echo "PART1_CNV_V2_PILOT_RUN_PASS"
} 2>&1 | tee -a "${LOG}"
