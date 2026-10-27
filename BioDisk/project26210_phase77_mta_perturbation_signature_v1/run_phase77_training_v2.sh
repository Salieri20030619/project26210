#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="/mnt/BioDisk/wangmingyuan/project26210_phase77_mta_perturbation_signature_v1"
R="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript"
PY="/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python"
DOWNLOADS="${ROOT}/downloads_v1"
TRAINING="${ROOT}/training_v1"
VALIDATION="${ROOT}/training_validation_v1"
LOGS="${ROOT}/logs_v1"

mkdir -p "${LOGS}"
exec > >(tee -a "${LOGS}/training_master_v2.log") 2>&1
on_error(){ code=$?; printf 'status=FAIL\nexit_code=%s\nfailed_at=%s\n' "${code}" "$(date -Is)" > "${ROOT}/TRAINING_FAILED"; exit "${code}"; }
trap on_error ERR

available_kib=$(awk '/MemAvailable/ {print $2}' /proc/meminfo)
minimum_kib=$((50*1024*1024))
printf 'time=%s\nMemAvailable_kib=%s\nminimum_kib=%s\n' "$(date -Is)" "${available_kib}" "${minimum_kib}" > "${LOGS}/memory_gate_training_v2.log"
[[ "${available_kib}" -ge "${minimum_kib}" ]] || exit 70

for file in GSE249082_AC-59_RawCounts.txt.gz GSE298469_RAW.tar GSE276571_MC38_raw_counts.txt.gz GSE166213_RAW.tar; do
  [[ -s "${DOWNLOADS}/${file}" ]] || { echo "Missing input ${file}" >&2; exit 72; }
done
for target in "${TRAINING}" "${VALIDATION}"; do
  [[ ! -e "${target}" && ! -e "${target}.part" ]] || { echo "Refusing existing output ${target}" >&2; exit 73; }
done

sha256sum "${DOWNLOADS}"/GSE249082_AC-59_RawCounts.txt.gz "${DOWNLOADS}"/GSE298469_RAW.tar "${DOWNLOADS}"/GSE276571_MC38_raw_counts.txt.gz "${DOWNLOADS}"/GSE166213_RAW.tar > "${ROOT}/DOWNLOAD_INPUT_SHA256_v1.tsv"
printf 'status=RUNNING\npid=%s\nstarted_at=%s\n' "$$" "$(date -Is)" > "${ROOT}/TRAINING_RUNNING"
"${R}" "${ROOT}/01_run_phase77_perturbation_training_v2.R" "${DOWNLOADS}" "${TRAINING}"
"${PY}" "${ROOT}/02_validate_phase77_training_v1.py" --training "${TRAINING}" --out "${VALIDATION}"
printf 'status=PASS\ncompleted_at=%s\n' "$(date -Is)" > "${ROOT}/TRAINING_COMPLETE"
echo "PHASE77 TRAINING PASS"
