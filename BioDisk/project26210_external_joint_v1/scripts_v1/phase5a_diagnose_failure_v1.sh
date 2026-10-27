#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="${ROOT}/logs/phase5a_diagnose_failure_${STAMP}.log"

conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase5a_diagnose_failure_v1.R" \
    2>&1 | tee "${LOG}"

printf '\nLOG=%s\n' "${LOG}"
