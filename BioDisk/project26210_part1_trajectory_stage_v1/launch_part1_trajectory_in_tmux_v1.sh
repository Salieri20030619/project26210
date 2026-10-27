#!/usr/bin/env bash
set -euo pipefail

STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_trajectory_stage_v1
PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
BASE="$ROOT/09_metabolism_trajectory"
OUT="$BASE/trajectory_v1"
SOCKET=/tmp/project26210_part1_trajectory_v1.sock
SESSION=part1_trajectory_v1
PY_BASE=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

if [[ ! -d "$OUT" || -e "$OUT/RUN_COMPLETE" || ! -e "$OUT/RUN_FAILED" || ! -e "$OUT/input/INPUT_COMPLETE" || -e "$OUT/R_METHODS_COMPLETE" ]]; then
  printf 'tmux launch invariant failed for %s\n' "$OUT" >&2
  exit 40
fi
if tmux -S "$SOCKET" has-session -t "$SESSION" 2>/dev/null; then
  printf 'Refusing duplicate tmux session %s\n' "$SESSION" >&2
  exit 41
fi
"$PROJECT/check_memory_gate.sh" part1_trajectory_tmux_prelaunch_v1
install -m 0755 "$STAGE/prepare_part1_trajectory_input_v1.py" "$BASE/prepare_part1_trajectory_input_v1.py"
install -m 0755 "$STAGE/run_part1_trajectory_r_methods_v1.R" "$BASE/run_part1_trajectory_r_methods_v1.R"
install -m 0755 "$STAGE/run_part1_trajectory_master_v1.sh" "$BASE/run_part1_trajectory_master_v1.sh"
install -m 0755 "$STAGE/correct_part1_trajectory_launch_records_v1.py" "$OUT/provenance/correct_part1_trajectory_launch_records_v1.py"
mv "$OUT/failure_status_v1.tsv" "$OUT/provenance/failed_attempt_2_status_v1.tsv"
mv "$OUT/RUN_FAILED" "$OUT/provenance/failed_attempt_2_RUN_FAILED"
cp "$OUT/background_status_v1.tsv" "$OUT/provenance/failed_attempt_2_background_status_v1.tsv"
cp "$OUT/logs/background_master_v1.log" "$OUT/provenance/failed_attempt_2_background_master_v1.log"
printf 'attempt\tcause\taction\n2\tfread auto-detected the blank scVI latent index header as a data row;11951 versus 11950\texplicit header=TRUE;reuse verified completed input;relaunch in dedicated tmux server\n' > "$OUT/provenance/failed_attempt_2_audit_v1.tsv"

started_at=$(date --iso-8601=seconds)
tmux -S "$SOCKET" new-session -d -s "$SESSION" "exec bash '$BASE/run_part1_trajectory_master_v1.sh' >> '$OUT/logs/background_master_v1.log' 2>&1"
pid=$(tmux -S "$SOCKET" list-panes -t "$SESSION" -F '#{pane_pid}')
printf 'started_at\tpid\tstate\tlaunch_backend\ttmux_socket\ttmux_session\toutput\tlog\tenvironment\n%s\t%s\trunning\ttmux\t%s\t%s\t%s\t%s\t%s\n' "$started_at" "$pid" "$SOCKET" "$SESSION" "$OUT" "$OUT/logs/background_master_v1.log" "/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1" > "$OUT/background_status_v1.tsv"
sleep 5
tmux -S "$SOCKET" has-session -t "$SESSION"
kill -0 "$pid"
"$PY_BASE" "$OUT/provenance/correct_part1_trajectory_launch_records_v1.py" --pid "$pid" --started-at "$started_at" --socket "$SOCKET"
tmux -S "$SOCKET" list-panes -t "$SESSION" -F 'session=#{session_name} pane_pid=#{pane_pid} pane_current_command=#{pane_current_command} pane_dead=#{pane_dead}'
printf 'PART1_TRAJECTORY_TMUX_LAUNCHED pid=%s socket=%s session=%s\n' "$pid" "$SOCKET" "$SESSION"
