#!/usr/bin/env bash
set -euo pipefail

STAGE=/mnt/BioDisk/wangmingyuan/project26210_part1_function_metabolism_stage_v1
PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="$PROJECT/results/part1_rebuild_v1"
CNMF="$ROOT/05_cnmf"
METAB_BASE="$ROOT/09_metabolism_trajectory"
METAB="$METAB_BASE/metabolism_v1"
SCENIC="$ROOT/07_pyscenic_v1"
R452=/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript
R453=/mnt/BioDisk/ludingyu/miniconda3/envs/r453/bin/Rscript
PY_BATCH=/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python
PY_SCFEA="$PROJECT/envs/part2_6_2_7_20260726_v1/python/bin/python"
SOCKET=/tmp/project26210_part1_metabolism_v1.sock
SESSION=p1_metabolism_v1

for executable in "$R452" "$R453" "$PY_BATCH" "$PY_SCFEA"; do
  [[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 2; }
done
[[ -f "$SCENIC/input/RUN_COMPLETE" ]] || { echo "Current 11,950-cell input is incomplete" >&2; exit 2; }

mkdir -p "$CNMF" "$METAB_BASE" "$METAB" "$METAB/logs"
if [[ -e "$CNMF/cnv_clone_program_function_v1/RUN_COMPLETE" ]]; then
  echo "Refusing to overwrite completed CNV/cNMF output" >&2; exit 3
fi
if find "$CNMF/cnv_clone_program_function_v1" -type f 2>/dev/null | grep -Ev '/(memory_gate_v1.txt|pseudobulk_sample_qc_v1.tsv)$' | grep -q .; then
  echo "Refusing to overwrite an unrecognized partial CNV/cNMF output" >&2; exit 3
fi
if find "$METAB" -type f -print -quit 2>/dev/null | grep -q .; then
  echo "Refusing to overwrite nonempty metabolism output" >&2; exit 3
fi

install -m 0755 "$STAGE/run_part1_cnv_cnmf_function_v1.R" "$CNMF/run_part1_cnv_cnmf_function_v1.R"
install -m 0755 "$STAGE/verify_part1_cnv_cnmf_function_v1.py" "$CNMF/verify_part1_cnv_cnmf_function_v1.py"
for script in run_part1_metabolism_auc_v1.R prepare_part1_scfea_input_v1.py run_part1_scfea_m168_v1.py analyze_part1_scfea_flux_v1.R verify_part1_metabolism_v1.py; do
  install -m 0755 "$STAGE/$script" "$METAB_BASE/$script"
done
install -m 0755 "$STAGE/run_part1_metabolism_waiter_v1.sh" "$METAB_BASE/run_part1_metabolism_waiter_v1.sh"

"$R452" "$CNMF/run_part1_cnv_cnmf_function_v1.R" >"$CNMF/cnv_clone_program_function_v1.log" 2>&1
"$PY_BATCH" "$CNMF/verify_part1_cnv_cnmf_function_v1.py" >"$CNMF/cnv_clone_program_function_verification_v1.log" 2>&1

tmux -S "$SOCKET" has-session -t "$SESSION" 2>/dev/null && { echo "tmux session already exists" >&2; exit 4; }
tmux -S "$SOCKET" new-session -d -s "$SESSION" "bash '$METAB_BASE/run_part1_metabolism_waiter_v1.sh' 2>&1 | tee '$METAB/logs/background_master_v1.log'"
tmux -S "$SOCKET" list-sessions
echo PART1_FUNCTION_AND_METABOLISM_INSTALL_LAUNCH_V1_PASS
