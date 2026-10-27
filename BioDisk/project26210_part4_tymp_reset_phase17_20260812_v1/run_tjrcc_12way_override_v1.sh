#!/usr/bin/env bash
set -euo pipefail

root=${1:?remote Phase 17 root is required}
python=${2:?cell2location Python path is required}
master_pid_file="$root/overnight_master_v1.pid"
state="$root/overnight_state_v1.tsv"
runner="$root/05_run_tjrcc12_cell2location_tymp_continuous_v1.py"
output_root="$root/tjrcc12_cell2location_tymp_continuous_v1"
preflight="$output_root/preflight_tjrcc12_tymp_continuous_v1.tsv"
override_log="$root/logs/tjrcc_12way_override_v1.log"
pid_manifest="$root/logs/tjrcc_12way_override_pids_v1.tsv"
override_pid_file="$root/tjrcc_12way_override_v1.pid"

stamp() { date -Is; }
record() { printf '%s\t%s\t%s\n' "$(stamp)" "$1" "$2" >> "$state"; }

master_pid=$(tr -cd '0-9' < "$master_pid_file")
master_cmd=$(ps -p "$master_pid" -o args= || true)
if [[ "$master_cmd" != *"run_overnight_phase17_tymp_continuous_v1.sh"* ]]; then
    echo "Refusing override: expected master PID is not the Phase 17 scheduler" >&2
    exit 2
fi
if [[ -f "$override_pid_file" ]]; then
    prior_override_pid=$(tr -cd '0-9' < "$override_pid_file")
    if kill -0 "$prior_override_pid" 2>/dev/null; then
        echo "Refusing duplicate override: PID $prior_override_pid is already active" >&2
        exit 2
    fi
fi
printf '%s\n' "$$" > "$override_pid_file"

# Stop only the scheduler. Its four existing mapping children continue running.
kill -STOP "$master_pid"
record tjrcc_parallel_override START

mapfile -t samples < <(tail -n +2 "$preflight" | cut -f1)
if [[ "${#samples[@]}" -ne 12 ]]; then
    record tjrcc_parallel_override FAIL_SAMPLE_DISCOVERY
    kill -CONT "$master_pid"
    exit 3
fi

printf 'timestamp\tsample\tpid\taction\n' > "$pid_manifest"
for sample in "${samples[@]}"; do
    if [[ -d "$output_root/$sample" ]]; then
        printf '%s\t%s\t%s\t%s\n' "$(stamp)" "$sample" "NA" "retain_existing" >> "$pid_manifest"
        continue
    fi
    "$python" "$runner" --sample "$sample" \
        > "$root/logs/tjrcc_${sample}_tymp_continuous_v1.log" 2>&1 &
    child_pid=$!
    printf '%s\t%s\t%s\t%s\n' "$(stamp)" "$sample" "$child_pid" "launched" >> "$pid_manifest"
done

while true; do
    completed=0
    failed=0
    terminal=0
    for sample in "${samples[@]}"; do
        status_file="$output_root/$sample/status.json"
        if [[ -f "$status_file" ]] && grep -q '"status": "completed"' "$status_file"; then
            ((completed+=1))
            ((terminal+=1))
        elif [[ -f "$status_file" ]] && grep -q '"status": "failed"' "$status_file"; then
            ((failed+=1))
            ((terminal+=1))
        fi
    done
    printf '%s completed=%d failed=%d terminal=%d\n' "$(stamp)" "$completed" "$failed" "$terminal" >> "$override_log"
    if (( terminal == 12 )); then
        break
    fi
    sleep 60
done

# The stopped scheduler is no longer needed once every mapping has terminated.
kill -TERM "$master_pid" 2>/dev/null || true

if (( failed > 0 )); then
    record tjrcc_parallel_override FAIL
    exit 4
fi

record tjrcc_parallel_override PASS
record tjrcc_bounded_parallel PASS_OVERRIDE_12WAY
record final_mapping_audit START
if "$python" "$root/06_audit_mapping_queue_v1.py" \
    --root "$root" --output "$root/mapping_queue_audit_v1" \
    >> "$root/logs/overnight_final_mapping_audit_v1.log" 2>&1; then
    record final_mapping_audit PASS
else
    record final_mapping_audit FAIL
    exit 5
fi

cat > "$root/OVERNIGHT_COMPLETE_v1.json.part" <<EOF
{
  "schema_version": "phase17_overnight_pipeline_v1",
  "status": "PASS",
  "completed_at": "$(stamp)",
  "design": "13-state Cell2location with one Malignant state; continuous TYMP primary exposure; continuous F1 auxiliary covariate",
  "mapping_scope": "GLMF1, GLMF2, and 12 TJ-RCC sections",
  "scheduler_override": "12-way TJ-RCC concurrency"
}
EOF
mv "$root/OVERNIGHT_COMPLETE_v1.json.part" "$root/OVERNIGHT_COMPLETE_v1.json"
record pipeline PASS
