#!/usr/bin/env bash
set -euo pipefail
ROOT=/mnt/BioDisk/wangmingyuan/project26210_phase74_virtual_knockout_enrichment_v1
GATE=/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
"$GATE" | tee "$ROOT/logs_v1/memory_gate_validate_v2.log"
"$PYTHON" "$ROOT/03_validate_phase74_v2.py" 2>&1 | tee "$ROOT/logs_v1/validate_v2.log"
date --iso-8601=seconds > "$ROOT/RUN_COMPLETE"
echo PHASE74_VALIDATION_RECOVERY_V2_PASS
