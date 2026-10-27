#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: run_tjrcc_spatial_programs_sample_v1.sh SAMPLE INPUT_H5AD OUTPUT_DIR SEED BOUNDARY_TSV" >&2
  exit 2
fi

sample_id="$1"
input_h5ad="$2"
output_dir="$3"
seed="$4"
boundary_tsv="$5"

project_root="/mnt/BioDisk/ludingyu/project26210"
activate_script="/mnt/BioDisk/wangmingyuan/project26210_envs/activate_part4_advanced_env_v1.sh"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
r_script="$script_dir/run_tjrcc_sparkx_spatial_programs_v1.R"
python_script="$script_dir/run_tjrcc_hotspot_spatial_programs_v1.py"
gate_script="$project_root/check_memory_gate.sh"

if [[ -f "$output_dir/status.json" ]] && rg -q '"status": "completed"' "$output_dir/status.json"; then
  echo "SPATIAL_PROGRAMS_${sample_id}_ALREADY_COMPLETE"
  exit 0
fi
[[ -s "$input_h5ad" ]] || { echo "Missing input: $input_h5ad" >&2; exit 3; }
[[ -s "$boundary_tsv" ]] || { echo "Missing boundary: $boundary_tsv" >&2; exit 3; }

source "$activate_script"
expected_python="/mnt/BioDisk/wangmingyuan/project26210_envs/part4_advanced_py_v1/bin/python"
expected_rscript="/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript"
if [[ "$(command -v python)" != "$expected_python" ]]; then
  echo "Unexpected Python executable: $(command -v python)" >&2
  exit 4
fi
if [[ "$(command -v Rscript)" != "$expected_rscript" ]]; then
  echo "Unexpected Rscript executable: $(command -v Rscript)" >&2
  exit 4
fi

mkdir -p "$output_dir"
log_file="$output_dir/run.log"
exec > >(tee -a "$log_file") 2>&1

echo "[$(date '+%F %T %z')] SPATIAL_PROGRAMS_${sample_id}_START"
echo "OUTPUT_DIR=$output_dir"
echo "PYTHON=$(command -v python)"
echo "RSCRIPT=$(command -v Rscript)"
"$gate_script"

Rscript "$r_script" "$sample_id" "$input_h5ad" "$output_dir" "$seed"
python "$python_script" "$sample_id" "$input_h5ad" "$output_dir" "$seed" "$boundary_tsv"

"$gate_script"
echo "[$(date '+%F %T %z')] SPATIAL_PROGRAMS_${sample_id}_ALL_PASS"
