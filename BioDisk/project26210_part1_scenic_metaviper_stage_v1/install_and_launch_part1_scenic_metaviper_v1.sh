#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
SCENIC="$ROOT/07_pyscenic_v1"
VIPER="$ROOT/08_metaviper_v1"
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_scenic_metaviper_stage_v1
SOCKET=/tmp/project26210_p1_scenic_metaviper_v1.sock
SESSION=p1_scenic_viper_v1

for target in "$SCENIC" "$VIPER"; do
  if [[ -e "$target" ]] && [[ -n "$(find "$target" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "Refusing to overwrite non-empty target: $target" >&2
    exit 1
  fi
done

mkdir -p "$SCENIC/input" "$SCENIC/logs" "$SCENIC/intermediate" "$SCENIC/tables" "$VIPER/logs"
install -m 0755 "$STAGE/prepare_part1_scenic_metaviper_input_v1.py" "$SCENIC/prepare_part1_scenic_metaviper_input_v1.py"
install -m 0755 "$STAGE/build_part1_pyscenic_looms_v1.py" "$SCENIC/build_part1_pyscenic_looms_v1.py"
install -m 0755 "$STAGE/run_part1_scenic_metaviper_overnight_v1.sh" "$SCENIC/run_part1_scenic_metaviper_overnight_v1.sh"
install -m 0755 "$STAGE/run_part1_metaviper_v1.R" "$VIPER/run_part1_metaviper_v1.R"

"$PROJECT/check_memory_gate.sh" launch_part1_scenic_metaviper_v1
tmux -S "$SOCKET" new-session -d -s "$SESSION" "cd '$PROJECT' && '$SCENIC/run_part1_scenic_metaviper_overnight_v1.sh' > '$SCENIC/logs/overnight_master_v1.log' 2>&1"
sleep 3
tmux -S "$SOCKET" has-session -t "$SESSION"
tmux -S "$SOCKET" list-panes -t "$SESSION" -F '#{session_name}\t#{pane_pid}\t#{pane_current_command}\t#{pane_dead}'
tail -n 30 "$SCENIC/logs/overnight_master_v1.log"
