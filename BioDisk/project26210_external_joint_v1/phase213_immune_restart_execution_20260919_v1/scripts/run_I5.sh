#!/usr/bin/env bash
set -euo pipefail
O=/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase213_immune_restart_execution_20260919_v1
PY=/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
stage=I5_sample_metrics
trap 'rc=$?; printf "stage=%s exit=%s\n" "$stage" "$rc" > "$O/I5_EXIT.txt"' EXIT
"$PY" -c "import json; assert json.load(open('$O/I4_directed_communication/STATUS.json'))['status']=='COMPLETE_VALIDATED'"
/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla "$O/scripts/05a_external_metrics.R" > "$O/I5_extract.log" 2>&1
stage=I5_description
"$PY" "$O/scripts/05b_external_description.py" > "$O/I5_description.log" 2>&1
stage=COMPLETE_I5
