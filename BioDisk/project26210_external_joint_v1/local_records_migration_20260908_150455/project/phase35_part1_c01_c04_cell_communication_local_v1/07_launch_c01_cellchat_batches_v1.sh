#!/usr/bin/env bash
set -euo pipefail

root=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase35_part1_c01_c04_cell_communication_v1
rscript=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
script="$root/provenance/04_run_c01_cellchat_v2.R"
input="$root/input_v3"
db="$root/resource_lock_v1/CellChatDB.human_v1.rds"
mkdir -p "$root/logs"

launch() {
  batch="$1"
  patients="$2"
  output="$root/c01_cellchat_${batch}_v1"
  log="$root/logs/c01_cellchat_${batch}_v1.log"
  pidfile="$root/logs/c01_cellchat_${batch}_v1.pid"
  if test -e "$output" || test -e "$pidfile"; then
    echo "Refusing existing batch target: $batch" >&2
    exit 1
  fi
  nohup "$rscript" "$script" "$input" "$db" "$output" "$patients" >"$log" 2>&1 &
  echo "$!" >"$pidfile.part"
  mv "$pidfile.part" "$pidfile"
  echo "$batch pid=$! patients=$patients"
}

launch batch01 "GSE156632__RCC2,GSE156632__RCC4,GSE156632__RCC5,GSE159115__SS_2023"
launch batch02 "GSE207493__RCC100,GSE207493__RCC104,GSE207493__RCC106,GSE207493__RCC112"
launch batch03 "GSE207493__RCC113,GSE207493__RCC86,GSE207493__RCC87"
launch batch04 "GSE207493__RCC94,GSE207493__RCC96,GSE207493__RCC99"
