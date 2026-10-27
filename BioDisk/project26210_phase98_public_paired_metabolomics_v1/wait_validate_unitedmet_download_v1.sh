#!/usr/bin/env bash
set -u

if [[ "$#" -ne 1 ]]; then
  echo "usage: wait_validate_unitedmet_download_v1.sh DOWNLOAD_PID" >&2
  exit 64
fi

phase98_root="/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1"
download_pid="$1"
printf '%s\n' "$$" > "${phase98_root}/wait_validate_pid_v1"

while kill -0 "${download_pid}" 2>/dev/null; do
  sleep 900
done

cd "${phase98_root}" || exit 65
python3 "${phase98_root}/validate_unitedmet_download_v1.py" "${phase98_root}" \
  > "${phase98_root}/validation_run_v1.log" 2>&1
validation_exit="$?"
printf 'download_pid\tvalidation_exit\n%s\t%s\n' \
  "${download_pid}" "${validation_exit}" \
  > "${phase98_root}/wait_validate_status_v1.tsv"
exit "${validation_exit}"
