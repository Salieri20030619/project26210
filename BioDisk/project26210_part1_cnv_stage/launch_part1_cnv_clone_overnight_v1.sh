#!/usr/bin/env bash
set -euo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
WORK="${PROJECT}/results/part1_rebuild_v1/04_subclones_tymp"
STATE="${WORK}/overnight_cnv_clone_v1"
RUNNER="${WORK}/run_part1_cnv_clone_overnight_v1.sh"
LOG="${PROJECT}/results/part1_rebuild_v1/logs/part1_cnv_clone_overnight_v1.log"
SOCKET=project26210_part1_cnv_clone_overnight_v1
SESSION=cnv_clone_overnight_v1

"${PROJECT}/check_memory_gate.sh"
test -x "$(command -v tmux)"
test -f "${RUNNER}"
mkdir -p "${STATE}"
if tmux -L "${SOCKET}" has-session -t "${SESSION}" 2>/dev/null; then
    echo "Session already exists: socket=${SOCKET} session=${SESSION}"
    exit 75
fi
printf '%s\tsocket=%s\tsession=%s\trunner=%s\tlog=%s\n' "$(date --iso-8601=seconds)" "${SOCKET}" "${SESSION}" "${RUNNER}" "${LOG}" >> "${STATE}/launch_history_v1.tsv"
tmux -L "${SOCKET}" new-session -d -s "${SESSION}" "exec bash '${RUNNER}' >> '${LOG}' 2>&1"
if tmux -L "${SOCKET}" has-session -t "${SESSION}" 2>/dev/null; then
    echo "PART1_CNV_CLONE_OVERNIGHT_LAUNCH_PASS socket=${SOCKET} session=${SESSION} log=${LOG}"
elif [[ -f "${STATE}/OVERNIGHT_COMPLETE" ]]; then
    echo "PART1_CNV_CLONE_OVERNIGHT_LAUNCH_AND_FAST_COMPLETE_PASS log=${LOG}"
else
    echo "Overnight session exited before a completion marker was written" >&2
    exit 1
fi
