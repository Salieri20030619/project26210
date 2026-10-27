#!/usr/bin/env bash
set -euo pipefail

STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_trajectory_stage_v1
PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
BASE="$ROOT/09_metabolism_trajectory"
OUT="$BASE/trajectory_v1"
PY_BASE=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
PY_TRAJ=/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1/bin/python
RS_TRAJ=/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1/bin/Rscript

if [[ -e "$OUT" ]]; then
  printf 'Refusing existing output directory: %s\n' "$OUT" >&2
  exit 20
fi
"$PROJECT/check_memory_gate.sh" part1_trajectory_prelaunch_v1
mkdir -p "$OUT/logs" "$OUT/qc" "$OUT/provenance"
install -m 0755 "$STAGE/prepare_part1_trajectory_input_v1.py" "$BASE/prepare_part1_trajectory_input_v1.py"
install -m 0755 "$STAGE/run_part1_trajectory_r_methods_v1.R" "$BASE/run_part1_trajectory_r_methods_v1.R"
install -m 0755 "$STAGE/run_part1_trajectory_python_methods_v1.py" "$BASE/run_part1_trajectory_python_methods_v1.py"
install -m 0755 "$STAGE/verify_part1_trajectory_v1.py" "$BASE/verify_part1_trajectory_v1.py"
install -m 0755 "$STAGE/run_part1_trajectory_master_v1.sh" "$BASE/run_part1_trajectory_master_v1.sh"
install -m 0755 "$STAGE/update_part1_records_for_trajectory_launch_v1.py" "$OUT/provenance/update_part1_records_for_trajectory_launch_v1.py"
cp "$PROJECT/PROJECT_STATUS.md" "$OUT/provenance/PROJECT_STATUS_before_trajectory_launch_v1.md"
cp "$ROOT/PART1_REBUILD_STATUS.md" "$OUT/provenance/PART1_REBUILD_STATUS_before_trajectory_launch_v1.md"
cp "$PROJECT/CODEX_HANDOFF.md" "$OUT/provenance/CODEX_HANDOFF_before_trajectory_launch_v1.md"
cp "$PROJECT/RESULTS_MANIFEST.tsv" "$OUT/provenance/RESULTS_MANIFEST_before_trajectory_launch_v1.tsv"
"$PROJECT/check_memory_gate.sh" part1_trajectory_record_hashes_v1 > "$OUT/qc/memory_gate_record_hashes_v1.tsv"

started_at=$(date --iso-8601=seconds)
nohup setsid bash "$BASE/run_part1_trajectory_master_v1.sh" > "$OUT/logs/background_master_v1.log" 2>&1 < /dev/null &
pid=$!
printf 'started_at\tpid\tstate\toutput\tlog\tenvironment\n%s\t%s\trunning\t%s\t%s\t%s\n' "$started_at" "$pid" "$OUT" "$OUT/logs/background_master_v1.log" "/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1" > "$OUT/background_status_v1.tsv"
sleep 2
if ! kill -0 "$pid" 2>/dev/null; then
  printf 'Trajectory process failed immediately; inspect %s\n' "$OUT/logs/background_master_v1.log" >&2
  exit 21
fi
"$PY_BASE" "$OUT/provenance/update_part1_records_for_trajectory_launch_v1.py" --pid "$pid" --started-at "$started_at"
ps -o pid=,ppid=,stat=,etime=,cmd= -p "$pid"
printf 'PART1_TRAJECTORY_BACKGROUND_LAUNCHED pid=%s log=%s\n' "$pid" "$OUT/logs/background_master_v1.log"
