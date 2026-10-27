#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT=/mnt/BioDisk/ludingyu/project26210
ROOT="${PROJECT}/results/part1_rebuild_v1"
WORK="${ROOT}/04_subclones_tymp"
STATE="${WORK}/overnight_cnv_clone_v1"
STATUS="${STATE}/overnight_status_v1.tsv"
LOCK="${STATE}/overnight.lock"

mkdir -p "${STATE}"
exec 9>"${LOCK}"
if ! flock -n 9; then
    echo "Another Part1 CNV-clone overnight controller is already running" >&2
    exit 75
fi
if [[ ! -f "${STATUS}" ]]; then
    printf 'timestamp\tstage\tstatus\tdetail\n' > "${STATUS}"
fi

CURRENT_STAGE=initialization
record_status() {
    local stage="$1"
    local status="$2"
    local detail="$3"
    printf '%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "${stage}" "${status}" "${detail}" >> "${STATUS}"
}
on_error() {
    local exit_code=$?
    record_status "${CURRENT_STAGE}" FAIL "exit_code=${exit_code}; stopped_without_running_later_stages"
    exit "${exit_code}"
}
trap on_error ERR

run_stage() {
    local stage="$1"
    local marker="$2"
    local command="$3"
    CURRENT_STAGE="${stage}"
    if [[ -f "${marker}" ]]; then
        record_status "${stage}" SKIP "verified_marker_exists=${marker}"
        return 0
    fi
    "${PROJECT}/check_memory_gate.sh" "overnight_${stage}"
    record_status "${stage}" RUNNING "command=${command}"
    bash "${command}"
    test -f "${marker}"
    record_status "${stage}" PASS "marker=${marker}"
}

record_status controller START "pid=$$; no_downloads; no_dependency_changes; fail_fast"
run_stage formal_cnv_clone_verification \
    "${WORK}/cnv_subclones_verification_v1/RUN_COMPLETE" \
    "${WORK}/run_part1_formal_cnv_subclones_verification_v1.sh"
run_stage cnv_clone_tymp_cnmf_association \
    "${WORK}/cnv_clone_tymp_cnmf_association_v1/RUN_COMPLETE" \
    "${WORK}/run_part1_cnv_clone_tymp_cnmf_association_v1.sh"
run_stage cnv_clone_tymp_cnmf_association_verification \
    "${WORK}/cnv_clone_tymp_cnmf_association_verification_v1/RUN_COMPLETE" \
    "${WORK}/run_part1_cnv_clone_tymp_cnmf_association_verification_v1.sh"

CURRENT_STAGE=terminal_snapshot
"${PROJECT}/check_memory_gate.sh" overnight_terminal_snapshot > "${STATE}/terminal_memory_gate_v1.txt"
test -f "${WORK}/cnv_subclones_verification_v1/RUN_COMPLETE"
test -f "${WORK}/cnv_clone_tymp_cnmf_association_v1/RUN_COMPLETE"
test -f "${WORK}/cnv_clone_tymp_cnmf_association_verification_v1/RUN_COMPLETE"
record_status terminal_snapshot PASS "all_three_required_RUN_COMPLETE_markers_present; manual_figure_review_pending"
date --iso-8601=seconds > "${STATE}/OVERNIGHT_COMPLETE"
record_status controller PASS "overnight_compute_complete; project_records_and_manual_figure_review_pending"
trap - ERR
echo PART1_CNV_CLONE_OVERNIGHT_V1_PASS
