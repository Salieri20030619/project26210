#!/usr/bin/env bash
set -euo pipefail

root=/mnt/BioDisk/wangmingyuan/project26210_part4_tymp_reset_phase17_20260812_v1
python=/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python
h5ad=/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad
old_assign=/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1/reference_v2/reference_v2_cell_assignments.tsv.gz
old_sig=/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1/cell2location_reference_v3/reference_v3_signatures_study_balanced_mean_counts.tsv.gz
bridge="$root/results_v2/tables/tymp_f1_cell_bridge_v1.tsv.gz"
ref="$root/reference_tymp_continuous_v1"
validation="$root/reference_tymp_continuous_validation_v1"
state="$root/overnight_state_v1.tsv"

stamp() { date -Is; }
record() { printf '%s\t%s\t%s\n' "$(stamp)" "$1" "$2" >> "$state"; }
fail() { record "$1" FAIL; exit 1; }
run_stage() {
    local name=$1
    shift
    record "$name" START
    if "$@" >> "$root/logs/overnight_${name}_v1.log" 2>&1; then
        record "$name" PASS
    else
        fail "$name"
    fi
}

if [[ -e "$root/OVERNIGHT_COMPLETE_v1.json" ]]; then
    echo "Completed marker already exists; refusing duplicate run" >&2
    exit 2
fi
printf 'timestamp\tstage\tstatus\n' > "$state"
record pipeline START

# The rejected two-state attempt is preserved at reference_tymp_v1.part and its log.
run_stage build_reference "$python" "$root/02b_build_continuous_tymp_reference_v1.py" \
    --h5ad "$h5ad" --old-assignments "$old_assign" --old-signature "$old_sig" \
    --bridge "$bridge" --memory-gate /mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh \
    --output "$ref"

run_stage validate_reference "$python" "$root/03b_validate_continuous_tymp_reference_v1.py" \
    --reference "$ref" --old-signature "$old_sig" --bridge "$bridge" --output "$validation"
grep -q '"status": "PASS"' "$validation/PHASE17_REFERENCE_TYMP_CONTINUOUS_INDEPENDENT_VALIDATION_v1.json" \
    || fail validation_gate
record validation_gate PASS

run_stage preflight_glmf "$python" "$root/04_run_glmf12_cell2location_tymp_continuous_v1.py" --preflight-only
run_stage preflight_tjrcc "$python" "$root/05_run_tjrcc12_cell2location_tymp_continuous_v1.py" --preflight-only

record glmf_parallel START
"$python" "$root/04_run_glmf12_cell2location_tymp_continuous_v1.py" --sample GLMF1 \
    > "$root/logs/glmf1_tymp_continuous_v1.log" 2>&1 &
p1=$!
"$python" "$root/04_run_glmf12_cell2location_tymp_continuous_v1.py" --sample GLMF2 \
    > "$root/logs/glmf2_tymp_continuous_v1.log" 2>&1 &
p2=$!
rc=0
wait "$p1" || rc=1
wait "$p2" || rc=1
[[ "$rc" -eq 0 ]] || fail glmf_parallel
record glmf_parallel PASS

mapfile -t samples < <("$python" -c \
    "import importlib.util; p='$root/05_run_tjrcc12_cell2location_tymp_continuous_v1.py'; s=importlib.util.spec_from_file_location('r',p); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); print('\\n'.join(m.SAMPLES))")
[[ "${#samples[@]}" -eq 12 ]] || fail tjrcc_sample_discovery
record tjrcc_sample_discovery PASS

record tjrcc_bounded_parallel START
running=0
failed=0
for sample in "${samples[@]}"; do
    "$python" "$root/05_run_tjrcc12_cell2location_tymp_continuous_v1.py" --sample "$sample" \
        > "$root/logs/tjrcc_${sample}_tymp_continuous_v1.log" 2>&1 &
    ((running+=1))
    if (( running >= 4 )); then
        wait -n || failed=1
        ((running-=1))
    fi
done
while (( running > 0 )); do
    wait -n || failed=1
    ((running-=1))
done
[[ "$failed" -eq 0 ]] || fail tjrcc_bounded_parallel
record tjrcc_bounded_parallel PASS

run_stage final_mapping_audit "$python" "$root/06_audit_mapping_queue_v1.py" \
    --root "$root" --output "$root/mapping_queue_audit_v1"

cat > "$root/OVERNIGHT_COMPLETE_v1.json.part" <<EOF
{
  "schema_version": "phase17_overnight_pipeline_v1",
  "status": "PASS",
  "completed_at": "$(stamp)",
  "design": "13-state Cell2location with one Malignant state; continuous TYMP primary exposure; continuous F1 auxiliary covariate",
  "mapping_scope": "GLMF1, GLMF2, and 12 TJ-RCC sections"
}
EOF
mv "$root/OVERNIGHT_COMPLETE_v1.json.part" "$root/OVERNIGHT_COMPLETE_v1.json"
record pipeline PASS
