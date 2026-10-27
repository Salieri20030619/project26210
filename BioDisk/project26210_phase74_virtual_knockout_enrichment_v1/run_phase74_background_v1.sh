#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/BioDisk/wangmingyuan/project26210_phase74_virtual_knockout_enrichment_v1
INPUT=/mnt/BioDisk/wangmingyuan/project26210_phase73_virtual_knockout_v1/dataset_runs_v1
GATE=/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh
RSCRIPT=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
PYTHON=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python

mkdir -p "$ROOT/logs_v1" "$ROOT/results_v1" "$ROOT/independent_validation_v1"

"$GATE" | tee "$ROOT/logs_v1/memory_gate_enrichment_v1.log"
"$RSCRIPT" --vanilla "$ROOT/01_run_phase74_enrichment_v1.R" "$INPUT" "$ROOT/enrichment_run_v1" 2>&1 | tee "$ROOT/logs_v1/enrichment_v1.log"

"$GATE" | tee "$ROOT/logs_v1/memory_gate_finalize_v1.log"
"$PYTHON" "$ROOT/02_finalize_phase74_v1.py" 2>&1 | tee "$ROOT/logs_v1/finalize_v1.log"

"$GATE" | tee "$ROOT/logs_v1/memory_gate_validate_v1.log"
"$PYTHON" "$ROOT/03_validate_phase74_v1.py" 2>&1 | tee "$ROOT/logs_v1/validate_v1.log"

date --iso-8601=seconds > "$ROOT/RUN_COMPLETE"
echo PHASE74_BACKGROUND_PIPELINE_PASS
