#!/usr/bin/env bash
set -Eeuo pipefail

phase98_root="/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1"
download_dir="${phase98_root}/downloads_v1"
run_complete="FALSE"
printf '%s\n' "$$" > "${phase98_root}/download_pid_v3"
rm -f "${phase98_root}/RUN_COMPLETE_v3"

record_exit() {
  exit_code="$?"
  if [[ "${run_complete}" == "TRUE" && "${exit_code}" -eq 0 ]]; then
    printf 'status\texit_code\nPASS\t0\n' > "${phase98_root}/download_status_v3.tsv"
    touch "${phase98_root}/RUN_COMPLETE_v3"
  else
    printf 'status\texit_code\nFAIL\t%s\n' "${exit_code}" > "${phase98_root}/download_status_v3.tsv"
    printf '%s\n' "${exit_code}" > "${phase98_root}/RUN_FAILED_v3"
  fi
}
trap record_exit EXIT

download_one() {
  file_name="$1"
  source_url="$2"
  if [[ -s "${download_dir}/${file_name}" ]]; then
    return 0
  fi
  curl -L --fail --retry 30 --retry-all-errors --retry-delay 10 --continue-at - \
    --silent --show-error \
    --output "${download_dir}/${file_name}.part" "${source_url}"
  mv "${download_dir}/${file_name}.part" "${download_dir}/${file_name}"
}

mkdir -p "${download_dir}"
download_one "matched_Harmonized_RawData_RC18.csv" "https://zenodo.org/api/records/11286535/files/matched_Harmonized_RawData_RC18.csv/content"
download_one "matched_Harmonized_Met_CPTAC_val.csv" "https://zenodo.org/api/records/11286535/files/matched_Harmonized_Met_CPTAC_val.csv/content"
download_one "matched_Harmonized_RawData_RC20.csv" "https://zenodo.org/api/records/11286535/files/matched_Harmonized_RawData_RC20.csv/content"
download_one "matched_Harmonized_Met_CPTAC.csv" "https://zenodo.org/api/records/11286535/files/matched_Harmonized_Met_CPTAC.csv/content"
download_one "matched_tpm_RC18.csv" "https://zenodo.org/api/records/11286535/files/matched_tpm_RC18.csv/content"
download_one "matched_tpm_CPTAC.csv" "https://zenodo.org/api/records/11286535/files/matched_tpm_CPTAC.csv/content"
download_one "matched_tpm_RC20.csv" "https://zenodo.org/api/records/11286535/files/matched_tpm_RC20.csv/content"
download_one "matched_tpm_CPTAC_val.csv" "https://zenodo.org/api/records/11286535/files/matched_tpm_CPTAC_val.csv/content"

test ! -e /mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1/validation_v2
python3 "${phase98_root}/validate_unitedmet_download_v2.py" "${phase98_root}" \
  > "${phase98_root}/validation_run_v3.log" 2>&1
run_complete="TRUE"
