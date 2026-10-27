#!/usr/bin/env bash
set -euo pipefail

PROJECT="/mnt/BioDisk/ludingyu/project26210"
SOCKET="project26210_part1_rebuild_v1"
SESSION="malignancy_cnv_full_v1"
RUNNER="${PROJECT}/results/part1_rebuild_v1/03_malignancy_cnv/run_part1_malignancy_cnv_full_v1.sh"

"${PROJECT}/check_memory_gate.sh"
if tmux -L "${SOCKET}" has-session -t "${SESSION}" 2>/dev/null
then
  echo "SESSION_ALREADY_EXISTS socket=${SOCKET} session=${SESSION}"
  exit 1
fi
tmux -L "${SOCKET}" new-session -d -s "${SESSION}" "bash '${RUNNER}'"
tmux -L "${SOCKET}" list-panes -t "${SESSION}" -F 'socket=project26210_part1_rebuild_v1 session=malignancy_cnv_full_v1 pane_pid=#{pane_pid} pane_dead=#{pane_dead} command=#{pane_current_command}'
