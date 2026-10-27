#!/usr/bin/env bash
set -euo pipefail
STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_trajectory_env_stage_v1
PROJECT=/mnt/BioDisk/ludingyu/project26210
OUT=/mnt/BioDisk/wangmingyuan/project26210_envs/part1_trajectory_v1
[[ ! -e "$OUT/RUN_COMPLETE" ]] || { echo "Refusing to overwrite completed trajectory environment record" >&2; exit 3; }
mkdir -p "$OUT/bin"
"$PROJECT/check_memory_gate.sh" part1_trajectory_environment_install_v1 | tee "$OUT/memory_gate_v1.tsv"
install -m 0755 "$STAGE/part1_trajectory_Rscript_v1.sh" "$OUT/bin/Rscript"
install -m 0755 "$STAGE/part1_trajectory_python_v1.sh" "$OUT/bin/python"
install -m 0755 "$STAGE/verify_part1_trajectory_environment_v1.py" "$OUT/verify_part1_trajectory_environment_v1.py"
"$OUT/bin/python" "$OUT/verify_part1_trajectory_environment_v1.py"
echo PART1_TRAJECTORY_ENVIRONMENT_V1_PASS
