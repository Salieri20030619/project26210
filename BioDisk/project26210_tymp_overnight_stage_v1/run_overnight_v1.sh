#!/usr/bin/env bash
set -euo pipefail
stage="/mnt/BioDisk/wangmingyuan/project26210_tymp_overnight_stage_v1"
mkdir -p "${stage}/logs"
status="${stage}/overnight_status_v1.tsv"
printf 'time\tstage\tstatus\tdetail\n' > "${status}"
printf '%s\ttrajectory\trunning\t230k-cell reference; myeloid subset DPT\n' "$(date --iso-8601=seconds)" >> "${status}"
if /mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python "${stage}/run_monocyte_spp1_trajectory_v1.py" > "${stage}/logs/trajectory_v1.log" 2>&1; then
  printf '%s\ttrajectory\tcompleted\tRUN_COMPLETE created\n' "$(date --iso-8601=seconds)" >> "${status}"
else
  code=$?
  printf '%s\ttrajectory\tfailed\texit_code=%s; see trajectory_v1.log\n' "$(date --iso-8601=seconds)" "${code}" >> "${status}"
  exit "${code}"
fi
