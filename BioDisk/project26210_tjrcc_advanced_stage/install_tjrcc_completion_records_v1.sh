#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
STAGE="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage"
RUN_ROOT="${PROJECT}/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
AUDIT="${RUN_ROOT}/preflight/project_records_before_completion_v1"

test "$(sha256sum "${PROJECT}/PROJECT_STATUS.md" | awk '{print $1}')" = "4b5a4564e845ef616c38ec76b3f933f6bc754419851de5015aa6ac4b3d02be79"
test "$(sha256sum "${PROJECT}/CODEX_HANDOFF.md" | awk '{print $1}')" = "bd7bad8f0f5a858a8eafd86fca6dce3c9f73966a54c067877b8b0cb27ed2da82"
test "$(sha256sum "${PROJECT}/RESULTS_MANIFEST.tsv" | awk '{print $1}')" = "63eb75da07940ee3fc15fb898f9c900145e2d1662aeee1a1f159e70d7d97fa0e"

mkdir -p "${AUDIT}"
cp "${PROJECT}/PROJECT_STATUS.md" "${AUDIT}/PROJECT_STATUS.before.md"
cp "${PROJECT}/CODEX_HANDOFF.md" "${AUDIT}/CODEX_HANDOFF.before.md"
cp "${PROJECT}/RESULTS_MANIFEST.tsv" "${AUDIT}/RESULTS_MANIFEST.before.tsv"
cp "${STAGE}/PROJECT_STATUS.final.md" "${PROJECT}/PROJECT_STATUS.md"
cp "${STAGE}/CODEX_HANDOFF.final.md" "${PROJECT}/CODEX_HANDOFF.md"
cp "${STAGE}/RESULTS_MANIFEST.final.tsv" "${PROJECT}/RESULTS_MANIFEST.tsv"
cp "${STAGE}/manual_figure_review_v2.tsv" "${RUN_ROOT}/verification_v2/manual_figure_review_v2.tsv"
sha256sum "${AUDIT}"/* "${PROJECT}/PROJECT_STATUS.md" "${PROJECT}/CODEX_HANDOFF.md" \
  "${PROJECT}/RESULTS_MANIFEST.tsv" "${RUN_ROOT}/verification_v2/manual_figure_review_v2.tsv" \
  >"${AUDIT}/installed_sha256_v1.tsv"
echo "TJRCC_COMPLETION_RECORDS_INSTALLED audit=${AUDIT}"
