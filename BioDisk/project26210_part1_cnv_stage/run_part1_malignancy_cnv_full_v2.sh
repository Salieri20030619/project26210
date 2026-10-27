#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv"
LOG="${PROJECT}/results/part1_rebuild_v1/logs/malignancy_cnv_full_v2.log"
STATUS="${ROOT}/malignancy_cnv_full_v2_status.tsv"

mkdir -p "${ROOT}" "$(dirname "${LOG}")"
{
  echo -e "stage\tstatus\ttime" > "${STATUS}"
  echo -e "pipeline\trunning\t$(date --iso-8601=seconds)" >> "${STATUS}"
  test -f "${ROOT}/cnv_v2_pilot/verification/RUN_COMPLETE"
  test -f "${ROOT}/cnv_unit_inputs_v1/RUN_COMPLETE"
  "${PROJECT}/check_memory_gate.sh"

  set +e
  bash "${ROOT}/run_part1_infercnv_cohort_v2.sh" &
  infer_pid=$!
  bash "${ROOT}/run_part1_copykat_cohort_v2.sh" &
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
    echo "PART1_MALIGNANCY_CNV_FULL_V2_FAILED infercnv=${infer_status} copykat=${copy_status}"
    exit 1
  fi

  date --iso-8601=seconds > "${ROOT}/CNV_COHORT_V2_RUN_COMPLETE"
  echo -e "pipeline\tcompleted\t$(date --iso-8601=seconds)" >> "${STATUS}"
  echo "PART1_MALIGNANCY_CNV_FULL_V2_RUN_ALL_PASS"
} 2>&1 | tee -a "${LOG}"
