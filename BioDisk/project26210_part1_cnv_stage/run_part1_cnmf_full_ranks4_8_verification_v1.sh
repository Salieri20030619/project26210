#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
RSCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript

"${PROJECT}/check_memory_gate.sh"
test "$(${RSCRIPT} -e 'cat(R.home())')" = "/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/lib/R"
test -f "${ROOT}/05_cnmf/full_ranks4_8_v1/RUN_COMPLETE"
nice -n 10 "${RSCRIPT}" "${ROOT}/05_cnmf/verify_part1_cnmf_full_ranks4_8_v1.R"
