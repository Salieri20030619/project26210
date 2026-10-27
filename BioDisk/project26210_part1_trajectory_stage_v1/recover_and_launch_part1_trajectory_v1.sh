#!/usr/bin/env bash
set -euo pipefail

STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_trajectory_stage_v1
PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
BASE="$ROOT/09_metabolism_trajectory"
OUT="$BASE/trajectory_v1"
PY_BASE=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

if [[ ! -d "$OUT" || ! -e "$OUT/RUN_FAILED" || -e "$OUT/input/INPUT_COMPLETE" || -e "$OUT/RUN_COMPLETE" ]]; then
  printf 'Recovery invariant failed for %s\n' "$OUT" >&2
  exit 30
fi
"$PROJECT/check_memory_gate.sh" part1_trajectory_recovery_prelaunch_v1
install -m 0755 "$STAGE/run_part1_trajectory_master_v1.sh" "$BASE/run_part1_trajectory_master_v1.sh"
attempt=1
while [[ -e "$OUT/provenance/failed_attempt_${attempt}_status_v1.tsv" ]]; do
  attempt=$((attempt + 1))
done
mv "$OUT/failure_status_v1.tsv" "$OUT/provenance/failed_attempt_${attempt}_status_v1.tsv"
mv "$OUT/RUN_FAILED" "$OUT/provenance/failed_attempt_${attempt}_RUN_FAILED"
cp "$OUT/background_status_v1.tsv" "$OUT/provenance/failed_attempt_${attempt}_background_status_v1.tsv"
cp "$OUT/logs/background_master_v1.log" "$OUT/provenance/failed_attempt_${attempt}_background_master_v1.log"
printf 'attempt\tcause\taction\n%s\tR wrapper rejected --version during provenance capture;no matrix input was read\treplaced with R -e version query and relaunched\n' "$attempt" > "$OUT/provenance/failed_attempt_${attempt}_audit_v1.tsv"

started_at=$(date --iso-8601=seconds)
nohup setsid bash "$BASE/run_part1_trajectory_master_v1.sh" > "$OUT/logs/background_master_v1.log" 2>&1 < /dev/null &
pid=$!
printf 'started_at\tpid\tstate\toutput\tlog\tenvironment\trecovered_from_attempt\n%s\t%s\trunning\t%s\t%s\t%s\t%s\n' "$started_at" "$pid" "$OUT" "$OUT/logs/background_master_v1.log" "/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1" "$attempt" > "$OUT/background_status_v1.tsv"
sleep 3
if ! kill -0 "$pid" 2>/dev/null; then
  printf 'Recovered trajectory process failed immediately; inspect %s\n' "$OUT/logs/background_master_v1.log" >&2
  exit 31
fi
"$PY_BASE" "$OUT/provenance/update_part1_records_for_trajectory_launch_v1.py" --pid "$pid" --started-at "$started_at"
ps -o pid=,ppid=,stat=,etime=,cmd= -p "$pid"
printf 'PART1_TRAJECTORY_RECOVERED_BACKGROUND_LAUNCHED pid=%s log=%s\n' "$pid" "$OUT/logs/background_master_v1.log"
