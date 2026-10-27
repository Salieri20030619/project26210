#!/usr/bin/env bash
set -u

phase98_root="/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1"
download_pid="$(cat "${phase98_root}/download_pid_v3")"
printf '%s\n' "$$" > "${phase98_root}/wait_download_pid_v3"
while kill -0 "${download_pid}" 2>/dev/null; do
  sleep 900
done
if [[ -s "${phase98_root}/download_status_v3.tsv" ]]; then
  cp "${phase98_root}/download_status_v3.tsv" "${phase98_root}/wait_download_status_v3.tsv"
else
  printf 'status\texit_code\nFAIL_STATUS_MISSING\tNA\n' > "${phase98_root}/wait_download_status_v3.tsv"
fi
