#!/usr/bin/env bash
set -euo pipefail
export MPLCONFIGDIR=/tmp/project26210_part1_trajectory_mplconfig
mkdir -p "$MPLCONFIGDIR"
exec /mnt/BioDisk/ludingyu/project26210/envs/part2_6_2_7_20260726_v1/python/bin/python "$@"
