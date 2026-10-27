#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_completion_stage_v1
R_SCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
DEST="${ROOT}/08_metaviper_v1"
FAILED="${DEST}/terminal_verification_preflight_attempt1_row_order_assumption_failure_v1"

"${PROJECT}/check_memory_gate.sh" part1_metaviper_terminal_verification_recovery_v1
if [[ -d "${DEST}/terminal_verification_v1" && ! -f "${DEST}/terminal_verification_v1/RUN_COMPLETE" ]]; then
  test ! -e "${FAILED}"
  mv "${DEST}/terminal_verification_v1" "${FAILED}"
fi
cp "${STAGE}/verify_part1_metaviper_terminal_v1.R" "${DEST}/verify_part1_metaviper_terminal_v1.R"
cmp -s "${STAGE}/verify_part1_metaviper_terminal_v1.R" "${DEST}/verify_part1_metaviper_terminal_v1.R"
"${R_SCRIPT}" "${DEST}/verify_part1_metaviper_terminal_v1.R" \
  > "${DEST}/terminal_verification_v1.log" 2>&1
grep -F 'PART1_METAVIPER_TERMINAL_VERIFICATION_V1_PASS' "${DEST}/terminal_verification_v1.log"
