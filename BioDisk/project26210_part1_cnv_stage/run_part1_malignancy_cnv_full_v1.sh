#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
PYTHON="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
LOG="${PROJECT}/results/part1_rebuild_v1/logs/malignancy_cnv_full_v1.log"
STATUS="${ROOT}/malignancy_cnv_full_v1_status.tsv"

mkdir -p "${ROOT}" "$(dirname "${LOG}")"
{
  echo -e "stage\tstatus\ttime" > "${STATUS}"
  echo -e "pipeline\trunning\t$(date --iso-8601=seconds)" >> "${STATUS}"
  "${PROJECT}/check_memory_gate.sh"
  "${PYTHON}" "${ROOT}/prepare_part1_cnv_unit_inputs_v2.py" \
    --preflight-dir "${ROOT}/malignancy_cnv_input_audit_v2" \
    --output-dir "${ROOT}/cnv_unit_inputs_v1"
  echo -e "input_preparation\tcompleted\t$(date --iso-8601=seconds)" >> "${STATUS}"

  set +e
  bash "${ROOT}/run_part1_infercnv_cohort_v1.sh" &
  infer_pid=$!
  bash "${ROOT}/run_part1_copykat_cohort_v1.sh" &
  copy_pid=$!
  wait "${infer_pid}"
  infer_status=$?
  wait "${copy_pid}"
  copy_status=$?
  set -e

  echo -e "infercnv\t${infer_status}\t$(date --iso-8601=seconds)" >> "${STATUS}"
  echo -e "copykat\t${copy_status}\t$(date --iso-8601=seconds)" >> "${STATUS}"
  if test "${infer_status}" -ne 0 || test "${copy_status}" -ne 0
  then
    echo -e "pipeline\tfailed\t$(date --iso-8601=seconds)" >> "${STATUS}"
    echo "PART1_MALIGNANCY_CNV_FULL_V1_FAILED infercnv=${infer_status} copykat=${copy_status}"
    exit 1
  fi

  date --iso-8601=seconds > "${ROOT}/CNV_COHORT_RUN_COMPLETE"
  echo -e "pipeline\tcompleted\t$(date --iso-8601=seconds)" >> "${STATUS}"
  echo "PART1_MALIGNANCY_CNV_FULL_V1_RUN_ALL_PASS"
} 2>&1 | tee -a "${LOG}"
