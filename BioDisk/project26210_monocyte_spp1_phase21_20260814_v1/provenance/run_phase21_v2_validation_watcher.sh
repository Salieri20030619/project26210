#!/usr/bin/env bash
set -u

root=/mnt/BioDisk/wangmingyuan/project26210_monocyte_spp1_phase21_20260814_v1
analysis_pid=2834434
watch_log="$root/logs/validation_watcher_v2.log"

echo "$(date --iso-8601=seconds) watcher_started analysis_pid=$analysis_pid" >> "$watch_log"
while ps -p "$analysis_pid" -o args= 2>/dev/null | grep -q '01_analyze_monocyte_spp1_trajectory_v2.py'; do
  sleep 60
done

status="$root/results_v2/PHASE21_MONOCYTE_SPP1_TRAJECTORY_STATUS_v1.json"
if ! test -f "$status" || ! grep -q 'PASS_ANALYSIS_COMPLETE' "$status"; then
  echo "$(date --iso-8601=seconds) analysis_not_pass; validation_not_started" >> "$watch_log"
  exit 1
fi

if test -e "$root/independent_validation_v1"; then
  echo "$(date --iso-8601=seconds) validation_target_exists; refusing_overwrite" >> "$watch_log"
  exit 1
fi

/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python \
  "$root/provenance/02_validate_monocyte_spp1_trajectory_v2.py" \
  --analysis "$root/results_v2" \
  --assignments /mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1/reference_v2/reference_v2_cell_assignments.tsv.gz \
  --output "$root/independent_validation_v1" \
  >> "$root/logs/validation_v1.log" 2>&1
rc=$?
echo "$(date --iso-8601=seconds) validation_finished exit_code=$rc" >> "$watch_log"
exit "$rc"
