#!/usr/bin/env bash
set -euo pipefail

RUN_ROOT="/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/advanced_full_v2"
STAGED="/mnt/BioDisk/wangmingyuan/project26210_tjrcc_advanced_stage/verify_tjrcc_advanced_full_v1.py"
OLD="${RUN_ROOT}/scripts/verify_tjrcc_advanced_full_v1.py"
NEW="${RUN_ROOT}/scripts/verify_tjrcc_advanced_full_v2.py"
AUDIT="${RUN_ROOT}/preflight/verification_v1_structural_na_failure_audit"

test "$(sha256sum "${OLD}" | awk '{print $1}')" = "9dc3afc0dcecd562a40f4315d0149161bd8c90f7d99255ecc179274cbaffa19b"
mkdir -p "${AUDIT}"
cp "${OLD}" "${AUDIT}/verify_tjrcc_advanced_full_v1.py"
cp "${RUN_ROOT}/verification_v1/verification_checks_v1.tsv" "${AUDIT}/"
cp "${RUN_ROOT}/verification_v1/verification_summary_v1.json" "${AUDIT}/"
cp "${RUN_ROOT}/verification_v1/failure_traceback.txt" "${AUDIT}/"
cp "${STAGED}" "${NEW}"
chmod 0644 "${NEW}"
sha256sum "${AUDIT}"/* "${NEW}" >"${AUDIT}/sha256_v1.tsv"
echo "TJRCC_VERIFIER_V2_INSTALLED new=${NEW} audit=${AUDIT}"
