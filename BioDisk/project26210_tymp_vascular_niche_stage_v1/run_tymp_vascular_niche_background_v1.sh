#!/usr/bin/env bash
set -euo pipefail

stage_dir="/mnt/BioDisk/wangmingyuan/project26210_tymp_vascular_niche_stage_v1"
log_dir="${stage_dir}/logs"
status_file="${stage_dir}/background_status_v1.tsv"
python_bin="/mnt/BioDisk/wangmingyuan/miniconda3/envs/batch_bench/bin/python"
analysis_script="${stage_dir}/run_tymp_vascular_niche_analysis_v1.py"

mkdir -p "${log_dir}"
printf 'time\tstage\tstatus\tdetail\n' > "${status_file}"
printf '%s\tlauncher\trunning\tbackground pipeline started\n' "$(date --iso-8601=seconds)" >> "${status_file}"

if ! /mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh >> "${log_dir}/memory_gate_launcher_v1.log" 2>&1; then
    printf '%s\tlauncher\tfailed\tinitial memory gate failed\n' "$(date --iso-8601=seconds)" >> "${status_file}"
    exit 1
fi

printf '%s\tanalysis\trunning\t14-sample targeted mechanism analysis\n' "$(date --iso-8601=seconds)" >> "${status_file}"
if "${python_bin}" "${analysis_script}" > "${log_dir}/analysis_v1.log" 2>&1; then
    printf '%s\tanalysis\tcompleted\tRUN_COMPLETE created\n' "$(date --iso-8601=seconds)" >> "${status_file}"
    printf '%s\tlauncher\tcompleted\tanalysis outputs ready for independent QC and figure/report generation\n' "$(date --iso-8601=seconds)" >> "${status_file}"
else
    exit_code=$?
    printf '%s\tanalysis\tfailed\texit_code=%s; see analysis_v1.log\n' "$(date --iso-8601=seconds)" "${exit_code}" >> "${status_file}"
    printf '%s\tlauncher\tfailed\tanalysis failed\n' "$(date --iso-8601=seconds)" >> "${status_file}"
    exit "${exit_code}"
fi
