#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path

PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
PILOT = ROOT / "cnv_v2_pilot"
OUT = PILOT / "verification"
INFER_UNITS = ["GSE156632__RCC1", "GSE178481__PR5", "GSE242299__P3"]
COPY_UNITS = ["RCC1t", "PR5_PTumor1", "T3_1"]
EXPECTED_NOT_ESTIMABLE = {("copykat", "T3_1")}
HIGH_INFER = "GSE156632__RCC1"
HIGH_COPY = "RCC1t"
HMM_FILES = [
    "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",
    "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat",
    "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat",
    "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_summary(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {row["parameter"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def normalized_tsv_rows(path: Path) -> tuple[tuple[str, ...], tuple[tuple[str, ...], ...]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = tuple(next(reader))
        rows = tuple(sorted(tuple(row) for row in reader))
    return header, rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    for method, units, retired_guard in (
        ("infercnv", INFER_UNITS, 10000),
        ("copykat", COPY_UNITS, 5000),
    ):
        for unit in units:
            unit_dir = PILOT / method / unit
            complete = unit_dir / "RUN_COMPLETE"
            not_estimable_file = unit_dir / "RUN_NOT_ESTIMABLE"
            status_file = unit_dir / "status.json"
            summary_file = unit_dir / "run_summary_v2.tsv"
            mapping_file = unit_dir / "cell_id_mapping_v2.tsv"
            gate_file = unit_dir / "memory_gate_log.tsv"
            expected_not_estimable = (method, unit) in EXPECTED_NOT_ESTIMABLE
            if expected_not_estimable:
                check(f"{method}:{unit}:not_estimable_marker", not_estimable_file.is_file(), str(not_estimable_file))
                check(f"{method}:{unit}:no_complete_marker", not complete.exists(), str(complete))
            else:
                check(f"{method}:{unit}:run_complete", complete.is_file(), str(complete))
            check(f"{method}:{unit}:status_exists", status_file.is_file(), str(status_file))
            if status_file.is_file():
                status = json.loads(status_file.read_text())
                expected_status = "not_estimable" if expected_not_estimable else "completed"
                check(f"{method}:{unit}:status_terminal", status.get("status") == expected_status, json.dumps(status, sort_keys=True))
                if expected_not_estimable:
                    check(f"{method}:{unit}:reason_class", status.get("reason_class") == "copykat_internal_dimension_failure_low_coverage", json.dumps(status, sort_keys=True))
                    check(f"{method}:{unit}:coverage_qc_status", status.get("gene_coverage_qc") == "warning_below_v1_guard", json.dumps(status, sort_keys=True))
            check(f"{method}:{unit}:summary_state", summary_file.is_file() == (not expected_not_estimable), str(summary_file))
            if summary_file.is_file():
                summary = read_summary(summary_file)
                n_genes = int(summary["n_genes_used"])
                expected_qc = "warning_below_v1_guard" if n_genes < retired_guard else "pass_v1_guard"
                check(f"{method}:{unit}:coverage_qc", summary.get("gene_coverage_qc") == expected_qc, f"n_genes={n_genes};qc={summary.get('gene_coverage_qc')};expected={expected_qc}")
                check(f"{method}:{unit}:guard_recorded", int(summary.get("retired_v1_guard_genes", -1)) == retired_guard, f"guard={summary.get('retired_v1_guard_genes')}")
            check(f"{method}:{unit}:mapping_exists", mapping_file.is_file(), str(mapping_file))
            if mapping_file.is_file():
                _, mapping_rows = normalized_tsv_rows(mapping_file)
                check(f"{method}:{unit}:mapping_nonempty", len(mapping_rows) > 0, f"rows={len(mapping_rows)}")
            check(f"{method}:{unit}:memory_gates", gate_file.is_file() and "\t0" in gate_file.read_text() and "PASS" in gate_file.read_text(), str(gate_file))

    v1_infer = ROOT / "infercnv_v1/patients" / HIGH_INFER / "infercnv_run"
    v2_infer = PILOT / "infercnv" / HIGH_INFER / "infercnv_run"
    for name in HMM_FILES:
        left = v1_infer / name
        right = v2_infer / name
        exists = left.is_file() and right.is_file()
        check(f"infercnv_high_control:{name}:exists", exists, f"v1={left};v2={right}")
        if exists:
            left_table = normalized_tsv_rows(left)
            right_table = normalized_tsv_rows(right)
            check(f"infercnv_high_control:{name}:exact_table", left_table == right_table, f"v1_rows={len(left_table[1])};v2_rows={len(right_table[1])}")

    v1_copy = ROOT / "copykat_v1/samples" / HIGH_COPY / "copykat_predictions_mapped_v1.tsv"
    v2_copy = PILOT / "copykat" / HIGH_COPY / "copykat_predictions_mapped_v2.tsv"
    copy_exists = v1_copy.is_file() and v2_copy.is_file()
    check("copykat_high_control:predictions_exist", copy_exists, f"v1={v1_copy};v2={v2_copy}")
    if copy_exists:
        check("copykat_high_control:predictions_exact", normalized_tsv_rows(v1_copy) == normalized_tsv_rows(v2_copy), f"v1_sha={sha256(v1_copy)};v2_sha={sha256(v2_copy)}")

    manifest_rows = []
    for path in sorted(p for p in PILOT.rglob("*") if p.is_file() and OUT not in p.parents):
        manifest_rows.append({"path": str(path), "size_bytes": str(path.stat().st_size), "sha256": sha256(path)})

    checks_file = OUT / "verification_checks_v1.tsv"
    with checks_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)
    with (OUT / "file_manifest_v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    failed = [row for row in checks if row["status"] != "PASS"]
    metadata = {
        "schema_version": "part1_cnv_v2_pilot_verification_v1",
        "created_at": datetime.now().astimezone().isoformat(),
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed),
        "checks_failed": len(failed),
        "files_hashed": len(manifest_rows),
        "scientific_parameters_changed_from_v1": False,
        "retired_guards": {"infercnv_min_genes": 10000, "copykat_min_genes": 5000},
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if failed:
        (OUT / "RUN_FAILED").write_text(datetime.now().astimezone().isoformat() + "\n")
        raise SystemExit(f"PART1_CNV_V2_PILOT_VERIFICATION_FAILED checks={len(checks)} failed={len(failed)}")
    prior_failed = OUT / "RUN_FAILED"
    if prior_failed.exists():
        prior_failed.replace(OUT / "RUN_FAILED_ATTEMPT1_EXPECTED_NOT_ESTIMABLE_RULE_MISSING")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(f"PART1_CNV_V2_PILOT_VERIFICATION_ALL_PASS checks={len(checks)} files={len(manifest_rows)}")


if __name__ == "__main__":
    main()
