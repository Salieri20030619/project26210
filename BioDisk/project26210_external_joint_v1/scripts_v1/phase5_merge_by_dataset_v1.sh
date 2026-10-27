#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"

cd "${ROOT}"

exec conda run --no-capture-output -n R_4.5.2 \
    Rscript --vanilla \
    "${ROOT}/scripts_v1/phase5_merge_by_dataset_v1.R"
