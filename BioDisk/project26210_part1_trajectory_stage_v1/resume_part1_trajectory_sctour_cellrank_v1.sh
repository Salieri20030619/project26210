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

if [[ ! -e "$OUT/RUN_FAILED" || ! -e "$OUT/R_METHODS_COMPLETE" || ! -s "$OUT/source_data/CytoTRACE2_predictions_by_study_v1.tsv.gz" || -e "$OUT/PYTHON_METHODS_COMPLETE" || -e "$OUT/RUN_COMPLETE" ]]; then
  printf 'Resume invariant failed for %s\n' "$OUT" >&2
  exit 50
fi
if tmux -S "$SOCKET" has-session -t "$SESSION" 2>/dev/null; then
  printf 'Refusing duplicate tmux session %s\n' "$SESSION" >&2
  exit 51
fi
"$PROJECT/check_memory_gate.sh" part1_trajectory_sctour_cellrank_resume_prelaunch_v1
install -m 0755 "$STAGE/run_part1_trajectory_python_methods_v1.py" "$BASE/run_part1_trajectory_python_methods_v1.py"
install -m 0755 "$STAGE/run_part1_trajectory_master_v1.sh" "$BASE/run_part1_trajectory_master_v1.sh"
install -m 0755 "$STAGE/update_part1_trajectory_resume_records_v1.py" "$OUT/provenance/update_part1_trajectory_resume_records_v1.py"
mv "$OUT/failure_status_v1.tsv" "$OUT/provenance/failed_attempt_3_status_v1.tsv"
mv "$OUT/RUN_FAILED" "$OUT/provenance/failed_attempt_3_RUN_FAILED"
cp "$OUT/background_status_v1.tsv" "$OUT/provenance/failed_attempt_3_background_status_v1.tsv"
cp "$OUT/logs/background_master_v1.log" "$OUT/provenance/failed_attempt_3_background_master_v1.log"
printf 'attempt\tcause\taction\n3\ttrajectory AnnData lacked scTour-required n_genes_by_counts\tcompute the field on the 2500-HVG raw-count matrix;reuse R and CytoTRACE2 checkpoints;resume scTour and CellRank\n' > "$OUT/provenance/failed_attempt_3_audit_v1.tsv"

started_at=$(date --iso-8601=seconds)
tmux -S "$SOCKET" new-session -d -s "$SESSION" "exec bash '$BASE/run_part1_trajectory_master_v1.sh' >> '$OUT/logs/background_master_v1.log' 2>&1"
pid=$(tmux -S "$SOCKET" list-panes -t "$SESSION" -F '#{pane_pid}')
printf 'started_at\tpid\tstate\tlaunch_backend\ttmux_socket\ttmux_session\tresume_scope\toutput\tlog\tenvironment\n%s\t%s\trunning\ttmux\t%s\t%s\tscTour_and_CellRank_only\t%s\t%s\t%s\n' "$started_at" "$pid" "$SOCKET" "$SESSION" "$OUT" "$OUT/logs/background_master_v1.log" "/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1" > "$OUT/background_status_v1.tsv"
sleep 5
tmux -S "$SOCKET" has-session -t "$SESSION"
kill -0 "$pid"
"$PY_BASE" "$OUT/provenance/update_part1_trajectory_resume_records_v1.py" --pid "$pid" --started-at "$started_at" --socket "$SOCKET"
tmux -S "$SOCKET" list-panes -t "$SESSION" -F 'session=#{session_name} pane_pid=#{pane_pid} pane_current_command=#{pane_current_command} pane_dead=#{pane_dead}'
printf 'PART1_TRAJECTORY_SCTOUR_CELLRANK_RESUMED pid=%s socket=%s session=%s\n' "$pid" "$SOCKET" "$SESSION"
