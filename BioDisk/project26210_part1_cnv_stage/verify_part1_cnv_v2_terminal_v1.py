#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.io import mmread


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
INFER_ROOT = ROOT / "infercnv_v2/patients"
COPY_ROOT = ROOT / "copykat_v2/samples"
COPY_INPUT_ROOT = ROOT / "cnv_unit_inputs_v1/copykat"
RECOVERY = ROOT / "cnv_v2_terminal_recovery_v1"
OUT = ROOT / "cnv_v2_terminal_verification_v1"
EXPECTED_RECOVERED = {"T2_1", "T2_2"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def wrapper_gene_count(unit: str) -> int:
    input_dir = COPY_INPUT_ROOT / unit
    with gzip.open(input_dir / "cell_id_mapping_v1.tsv.gz", "rt", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(input_dir / "counts_genes_by_cells_v1.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsc()
    if matrix.shape[1] != len(mapping):
        raise RuntimeError(f"{unit}: matrix/mapping mismatch")
    keep_cells = np.asarray((matrix > 0).sum(axis=0)).ravel() >= 200
    matrix = matrix[:, keep_cells]
    return int((np.asarray((matrix > 0).sum(axis=1)).ravel() / matrix.shape[1] > 0.05).sum())


def internal_counts(unit: str) -> tuple[int, int, int]:
    path = COPY_ROOT / unit / "copykat_run" / f"{unit}_copykat_raw_results_bin_by_cell.txt"
    with path.open() as handle:
        ids = handle.readline().rstrip("\n").split("\t")[3:]
    return len(ids), sum(x.startswith("Q") for x in ids), sum(x.startswith("R") for x in ids)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})

    recovery_marker = RECOVERY / "RUN_COMPLETE"
    recovery_table = RECOVERY / "resolved_unit_status_v1.tsv"
    check("recovery:run_complete", recovery_marker.is_file(), str(recovery_marker))
    check("recovery:status_table", recovery_table.is_file(), str(recovery_table))
    rows = read_tsv(recovery_table) if recovery_table.is_file() else []
    keys = [(row.get("method", ""), row.get("unit_id", "")) for row in rows]
    check("recovery:75_unique_units", len(rows) == 75 and len(set(keys)) == 75, f"rows={len(rows)};unique={len(set(keys))}")

    infer_units = sorted(path.name for path in INFER_ROOT.iterdir() if path.is_dir())
    copy_units = sorted(path.name for path in COPY_ROOT.iterdir() if path.is_dir())
    check("infercnv:expected_unit_count", len(infer_units) == 34, f"n={len(infer_units)}")
    check("copykat:expected_unit_count", len(copy_units) == 41, f"n={len(copy_units)}")

    critical: set[Path] = {recovery_marker, recovery_table, RECOVERY / "run_metadata_v1.json", RECOVERY / "evidence_file_manifest_v1.tsv"}
    for unit in infer_units:
        unit_dir = INFER_ROOT / unit
        status_path = unit_dir / "status.json"
        marker = unit_dir / "RUN_COMPLETE"
        summary = unit_dir / "run_summary_v2.tsv"
        mapping = unit_dir / "cell_id_mapping_v2.tsv"
        hmm = unit_dir / "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"
        status = read_json(status_path) if status_path.is_file() else {}
        check(f"infercnv:{unit}:terminal", status.get("status") == "completed" and marker.is_file(), json.dumps(status, sort_keys=True))
        check(f"infercnv:{unit}:core_files", all(x.is_file() for x in (summary, mapping, hmm)), f"summary={summary.is_file()};mapping={mapping.is_file()};hmm={hmm.is_file()}")
        critical.update({status_path, marker, summary, mapping, hmm})

    completed = set()
    existing_not_estimable = set()
    failed = set()
    for unit in copy_units:
        unit_dir = COPY_ROOT / unit
        status_path = unit_dir / "status.json"
        status = read_json(status_path) if status_path.is_file() else {}
        state = status.get("status")
        if state == "completed":
            completed.add(unit)
            marker = unit_dir / "RUN_COMPLETE"
            summary = unit_dir / "run_summary_v2.tsv"
            mapping = unit_dir / "cell_id_mapping_v2.tsv"
            prediction = unit_dir / "copykat_predictions_mapped_v2.tsv"
            check(f"copykat:{unit}:completed_files", all(x.is_file() for x in (marker, summary, mapping, prediction)), f"marker={marker.is_file()};summary={summary.is_file()};mapping={mapping.is_file()};prediction={prediction.is_file()}")
            if prediction.is_file():
                pred_rows = read_tsv(prediction)
                exact = [row.get("exact_cell_id", "") for row in pred_rows]
                tools = [row.get("tool_id", "") for row in pred_rows]
                check(f"copykat:{unit}:prediction_mapping", bool(pred_rows) and all(exact) and len(exact) == len(set(exact)) and len(tools) == len(set(tools)), f"rows={len(pred_rows)};exact_unique={len(set(exact))};tool_unique={len(set(tools))}")
            critical.update({status_path, marker, summary, mapping, prediction})
        elif state == "not_estimable":
            existing_not_estimable.add(unit)
            marker = unit_dir / "RUN_NOT_ESTIMABLE"
            check(f"copykat:{unit}:not_estimable_files", marker.is_file() and (unit_dir / "ERROR.txt").is_file(), f"marker={marker.is_file()}")
            check(f"copykat:{unit}:not_estimable_reason", status.get("reason_class") == "copykat_internal_dimension_failure_low_coverage", json.dumps(status, sort_keys=True))
            critical.update({status_path, marker, unit_dir / "ERROR.txt", unit_dir / "cell_id_mapping_v2.tsv"})
        elif state == "failed":
            failed.add(unit)
            critical.update({status_path, unit_dir / "ERROR.txt", unit_dir / "cell_id_mapping_v2.tsv"})
        else:
            check(f"copykat:{unit}:known_terminal_state", False, str(state))

    check("copykat:original_terminal_totals", len(completed) == 38 and existing_not_estimable == {"T3_1"} and failed == EXPECTED_RECOVERED, f"completed={len(completed)};not_estimable={sorted(existing_not_estimable)};failed={sorted(failed)}")

    recovered_rows = {row["unit_id"]: row for row in rows if row.get("resolution") == "reclassify_v2_failure_in_recovery_layer"}
    check("recovery:exact_recovered_units", set(recovered_rows) == EXPECTED_RECOVERED, str(sorted(recovered_rows)))
    for unit in sorted(EXPECTED_RECOVERED):
        row = recovered_rows.get(unit, {})
        unit_dir = COPY_ROOT / unit
        n_genes = wrapper_gene_count(unit)
        n_internal, n_query, n_reference = internal_counts(unit)
        message = read_json(unit_dir / "status.json").get("message", "")
        cutree_error = "'k'" in message and "'h'" in message
        check(f"recovery:{unit}:low_coverage", n_genes < 5000 and row.get("n_genes_wrapper_filter") == str(n_genes), f"n_genes={n_genes};recorded={row.get('n_genes_wrapper_filter')}")
        check(f"recovery:{unit}:internal_query_insufficient", n_query < 11 and row.get("n_internal_query") == str(n_query), f"internal={n_internal};query={n_query};reference={n_reference};recorded={row.get('n_internal_query')}")
        check(f"recovery:{unit}:cutree_error", cutree_error and row.get("reason_class") == "copykat_internal_clustering_undefined_low_coverage", f"message={message};reason={row.get('reason_class')}")
        critical.update({
            unit_dir / "copykat_run" / f"{unit}_copykat_raw_results_bin_by_cell.txt",
            COPY_INPUT_ROOT / unit / "counts_genes_by_cells_v1.mtx.gz",
            COPY_INPUT_ROOT / unit / "cell_id_mapping_v1.tsv.gz",
        })

    resolved_counts: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row.get("method", ""), row.get("resolved_status", ""))
        resolved_counts[key] = resolved_counts.get(key, 0) + 1
    expected_counts = {("infercnv", "completed"): 34, ("copykat", "completed"): 38, ("copykat", "not_estimable"): 3}
    check("recovery:resolved_totals", resolved_counts == expected_counts, f"actual={resolved_counts};expected={expected_counts}")

    evidence_manifest = RECOVERY / "evidence_file_manifest_v1.tsv"
    if evidence_manifest.is_file():
        evidence_rows = read_tsv(evidence_manifest)
        evidence_ok = True
        for row in evidence_rows:
            path = Path(row["path"])
            if not path.is_file() or str(path.stat().st_size) != row["size_bytes"] or sha256(path) != row["sha256"]:
                evidence_ok = False
                break
        check("recovery:evidence_manifest_integrity", evidence_ok, f"files={len(evidence_rows)}")

    checks_file = OUT / "verification_checks_v1.tsv"
    with checks_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)

    missing_critical = sorted(str(path) for path in critical if not path.is_file())
    check("verification:critical_files_present", not missing_critical, ";".join(missing_critical))
    with checks_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "status", "detail"], delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)

    manifest_rows = []
    for path in sorted(path for path in critical if path.is_file()):
        manifest_rows.append({"path": str(path), "size_bytes": str(path.stat().st_size), "sha256": sha256(path)})
    with (OUT / "critical_file_manifest_v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    failed_checks = [row for row in checks if row["status"] != "PASS"]
    metadata = {
        "schema_version": "part1_cnv_v2_terminal_verification_v1",
        "created_at": datetime.now().astimezone().isoformat(),
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failed_checks),
        "checks_failed": len(failed_checks),
        "critical_files_hashed": len(manifest_rows),
        "resolved_totals": {"infercnv_completed": 34, "copykat_completed": 38, "copykat_not_estimable": 3},
        "biological_conclusion_generated": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    if failed_checks:
        (OUT / "RUN_FAILED").write_text(datetime.now().astimezone().isoformat() + "\n")
        raise SystemExit(f"PART1_CNV_V2_TERMINAL_VERIFICATION_FAILED checks={len(checks)} failed={len(failed_checks)}")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(f"PART1_CNV_V2_TERMINAL_VERIFICATION_ALL_PASS checks={len(checks)} files={len(manifest_rows)}")


if __name__ == "__main__":
    main()
