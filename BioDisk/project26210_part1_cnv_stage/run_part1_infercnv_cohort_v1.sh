#!/usr/bin/env bash
set -uo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
INPUT_ROOT="${ROOT}/cnv_unit_inputs_v1/infercnv"
OUTPUT_ROOT="${ROOT}/infercnv_v1/patients"
RSCRIPT="/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript"
UNIT_SCRIPT="${ROOT}/run_part1_infercnv_unit_v1.R"
LOG="${PROJECT}/results/part1_rebuild_v1/logs/infercnv_cohort_v1.log"
FAIL_FILE="${ROOT}/infercnv_v1/failed_units_v1.tsv"

mkdir -p "${OUTPUT_ROOT}" "$(dirname "${LOG}")"
: > "${FAIL_FILE}"

failures=0
completed=0
skipped=0
{
  echo -e "started_at\t$(date --iso-8601=seconds)"
  for input_dir in "${INPUT_ROOT}"/*
  do
    test -d "${input_dir}" || continue
    unit_id="$(basename "${input_dir}")"
    output_dir="${OUTPUT_ROOT}/${unit_id}"
    if test -f "${output_dir}/RUN_COMPLETE"
    then
      echo "INFERCNV_SKIP_COMPLETED unit=${unit_id}"
      skipped=$((skipped + 1))
      continue
    fi
    mkdir -p "${output_dir}"
    echo "INFERCNV_START unit=${unit_id} time=$(date --iso-8601=seconds)"
    if "${PROJECT}/check_memory_gate.sh" && \
       "${RSCRIPT}" --vanilla "${UNIT_SCRIPT}" \
         --input-dir "${input_dir}" \
         --output-dir "${output_dir}" \
         --threads 8
    then
      completed=$((completed + 1))
      echo "INFERCNV_DONE unit=${unit_id} time=$(date --iso-8601=seconds)"
    else
      failures=$((failures + 1))
      echo -e "${unit_id}\tfailed\t$(date --iso-8601=seconds)" >> "${FAIL_FILE}"
      echo "INFERCNV_FAILED unit=${unit_id} time=$(date --iso-8601=seconds)"
    fi
  done
  echo "INFERCNV_COHORT_COUNTS completed=${completed} skipped=${skipped} failed=${failures}"
  if test "${failures}" -eq 0
  then
    date --iso-8601=seconds > "${ROOT}/infercnv_v1/RUN_COMPLETE"
    echo "PART1_INFERCNV_COHORT_V1_ALL_PASS"
  else
    echo "PART1_INFERCNV_COHORT_V1_FAILED failures=${failures}"
    exit 1
  fi
} 2>&1 | tee -a "${LOG}"
