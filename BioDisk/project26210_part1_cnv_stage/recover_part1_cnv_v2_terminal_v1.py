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
OUT = ROOT / "cnv_v2_terminal_recovery_v1"
RECOVERABLE_CUTREE_UNITS = {"T2_1", "T2_2"}
COPYKAT_RETIRED_GENE_GUARD = 5000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def read_summary(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        return {row["parameter"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def wrapper_filter_counts(unit: str) -> tuple[int, int, int]:
    input_dir = COPY_INPUT_ROOT / unit
    with gzip.open(input_dir / "cell_id_mapping_v1.tsv.gz", "rt", newline="") as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    with gzip.open(input_dir / "counts_genes_by_cells_v1.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsc()
    if matrix.shape[1] != len(mapping):
        raise RuntimeError(f"{unit}: matrix/mapping column mismatch {matrix.shape[1]} != {len(mapping)}")
    keep_cells = np.asarray((matrix > 0).sum(axis=0)).ravel() >= 200
    kept_roles = np.asarray([row["role"] for row in mapping], dtype=object)[keep_cells]
    filtered = matrix[:, keep_cells]
    keep_genes = np.asarray((filtered > 0).sum(axis=1)).ravel() / filtered.shape[1] > 0.05
    return int(keep_genes.sum()), int((kept_roles == "query").sum()), int((kept_roles == "reference").sum())


def internal_bin_counts(unit: str) -> tuple[int, int, int]:
    path = COPY_ROOT / unit / "copykat_run" / f"{unit}_copykat_raw_results_bin_by_cell.txt"
    with path.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    if header[:3] != ["chrom", "chrompos", "abspos"]:
        raise RuntimeError(f"{unit}: unexpected CopyKAT bin header")
    cell_ids = header[3:]
    return len(cell_ids), sum(x.startswith("Q") for x in cell_ids), sum(x.startswith("R") for x in cell_ids)


def cutree_partition_undefined(message: str) -> bool:
    lowered = message.lower()
    return (
        ("'k'" in lowered and "'h'" in lowered and ("must" in lowered or "必需设定一个" in message))
        or "must specify either 'k' or 'h'" in lowered
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    infer_units = sorted(path.name for path in INFER_ROOT.iterdir() if path.is_dir())
    copy_units = sorted(path.name for path in COPY_ROOT.iterdir() if path.is_dir())
    if len(infer_units) != 34:
        raise RuntimeError(f"Expected 34 inferCNV units, found {len(infer_units)}")
    if len(copy_units) != 41:
        raise RuntimeError(f"Expected 41 CopyKAT units, found {len(copy_units)}")

    rows: list[dict[str, str]] = []
    evidence: set[Path] = set()

    for unit in infer_units:
        unit_dir = INFER_ROOT / unit
        status_path = unit_dir / "status.json"
        status = read_json(status_path)
        marker = unit_dir / "RUN_COMPLETE"
        summary_path = unit_dir / "run_summary_v2.tsv"
        if status.get("status") != "completed" or not marker.is_file() or not summary_path.is_file():
            raise RuntimeError(f"inferCNV unit not complete: {unit}")
        summary = read_summary(summary_path)
        rows.append({
            "method": "infercnv",
            "unit_id": unit,
            "original_status": "completed",
            "resolved_status": "completed",
            "resolution": "use_v2_completed_output",
            "reason_class": "completed_without_recovery",
            "n_genes_wrapper_filter": summary.get("n_genes_used", ""),
            "n_query_wrapper_filter": summary.get("n_query_cells", ""),
            "n_reference_wrapper_filter": summary.get("n_reference_cells", ""),
            "n_internal_cells": "",
            "n_internal_query": "",
            "n_internal_reference": "",
        })
        evidence.update({
            status_path,
            marker,
            summary_path,
            unit_dir / "cell_id_mapping_v2.tsv",
            unit_dir / "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings",
            unit_dir / "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat",
            unit_dir / "infercnv_run/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat",
        })

    recovered_units: set[str] = set()
    original_not_estimable: set[str] = set()
    completed_copy = 0
    for unit in copy_units:
        unit_dir = COPY_ROOT / unit
        status_path = unit_dir / "status.json"
        status = read_json(status_path)
        original_status = str(status.get("status"))
        common = {
            "method": "copykat",
            "unit_id": unit,
            "original_status": original_status,
            "n_internal_cells": "",
            "n_internal_query": "",
            "n_internal_reference": "",
        }
        if original_status == "completed":
            marker = unit_dir / "RUN_COMPLETE"
            summary_path = unit_dir / "run_summary_v2.tsv"
            predictions = unit_dir / "copykat_predictions_mapped_v2.tsv"
            if not marker.is_file() or not summary_path.is_file() or not predictions.is_file():
                raise RuntimeError(f"CopyKAT completed unit missing terminal output: {unit}")
            summary = read_summary(summary_path)
            rows.append(common | {
                "resolved_status": "completed",
                "resolution": "use_v2_completed_output",
                "reason_class": "completed_without_recovery",
                "n_genes_wrapper_filter": summary.get("n_genes_used", ""),
                "n_query_wrapper_filter": summary.get("n_query_after_filter", ""),
                "n_reference_wrapper_filter": summary.get("n_reference_after_filter", ""),
            })
            completed_copy += 1
            evidence.update({status_path, marker, summary_path, predictions, unit_dir / "cell_id_mapping_v2.tsv"})
            continue

        if original_status == "not_estimable":
            marker = unit_dir / "RUN_NOT_ESTIMABLE"
            if not marker.is_file() or status.get("reason_class") != "copykat_internal_dimension_failure_low_coverage":
                raise RuntimeError(f"Unexpected existing not-estimable state: {unit}")
            n_genes, n_query, n_reference = wrapper_filter_counts(unit)
            rows.append(common | {
                "resolved_status": "not_estimable",
                "resolution": "retain_v2_not_estimable",
                "reason_class": str(status.get("reason_class")),
                "n_genes_wrapper_filter": str(n_genes),
                "n_query_wrapper_filter": str(n_query),
                "n_reference_wrapper_filter": str(n_reference),
            })
            original_not_estimable.add(unit)
            evidence.update({status_path, marker, unit_dir / "ERROR.txt", unit_dir / "cell_id_mapping_v2.tsv"})
            continue

        if original_status == "failed" and unit in RECOVERABLE_CUTREE_UNITS:
            message = str(status.get("message", ""))
            n_genes, n_query, n_reference = wrapper_filter_counts(unit)
            n_internal, n_internal_query, n_internal_reference = internal_bin_counts(unit)
            if not cutree_partition_undefined(message):
                raise RuntimeError(f"{unit}: failure is not the audited cutree condition: {message}")
            if n_genes >= COPYKAT_RETIRED_GENE_GUARD:
                raise RuntimeError(f"{unit}: gene coverage is not low: {n_genes}")
            if n_internal_query >= 11:
                raise RuntimeError(f"{unit}: internal query count does not support recovery rule: {n_internal_query}")
            rows.append(common | {
                "resolved_status": "not_estimable",
                "resolution": "reclassify_v2_failure_in_recovery_layer",
                "reason_class": "copykat_internal_clustering_undefined_low_coverage",
                "n_genes_wrapper_filter": str(n_genes),
                "n_query_wrapper_filter": str(n_query),
                "n_reference_wrapper_filter": str(n_reference),
                "n_internal_cells": str(n_internal),
                "n_internal_query": str(n_internal_query),
                "n_internal_reference": str(n_internal_reference),
            })
            recovered_units.add(unit)
            evidence.update({
                status_path,
                unit_dir / "ERROR.txt",
                unit_dir / "cell_id_mapping_v2.tsv",
                unit_dir / "copykat_run" / f"{unit}_copykat_raw_results_gene_by_cell.txt",
                unit_dir / "copykat_run" / f"{unit}_copykat_raw_results_bin_by_cell.txt",
                COPY_INPUT_ROOT / unit / "unit_metadata_v1.json",
                COPY_INPUT_ROOT / unit / "counts_genes_by_cells_v1.mtx.gz",
                COPY_INPUT_ROOT / unit / "cell_id_mapping_v1.tsv.gz",
            })
            continue

        raise RuntimeError(f"Unexpected CopyKAT terminal state: {unit}={original_status}")

    if completed_copy != 38 or original_not_estimable != {"T3_1"} or recovered_units != RECOVERABLE_CUTREE_UNITS:
        raise RuntimeError(
            f"Unexpected CopyKAT resolution totals: completed={completed_copy}; "
            f"original_not_estimable={sorted(original_not_estimable)}; recovered={sorted(recovered_units)}"
        )

    fields = [
        "method", "unit_id", "original_status", "resolved_status", "resolution", "reason_class",
        "n_genes_wrapper_filter", "n_query_wrapper_filter", "n_reference_wrapper_filter",
        "n_internal_cells", "n_internal_query", "n_internal_reference",
    ]
    with (OUT / "resolved_unit_status_v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    missing = sorted(str(path) for path in evidence if not path.is_file())
    if missing:
        raise RuntimeError("Missing evidence files: " + "; ".join(missing))
    with (OUT / "evidence_file_manifest_v1.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        for path in sorted(evidence):
            writer.writerow({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})

    metadata = {
        "schema_version": "part1_cnv_v2_terminal_recovery_v1",
        "created_at": datetime.now().astimezone().isoformat(),
        "infercnv": {"completed": 34, "not_estimable": 0, "failed": 0},
        "copykat": {
            "completed": 38,
            "not_estimable_original_v2": 1,
            "not_estimable_recovered": 2,
            "failed_unresolved": 0,
        },
        "recovered_units": sorted(recovered_units),
        "recovery_rule": {
            "required_original_status": "failed",
            "required_error": "CopyKAT cutree partition undefined because no k=2..4 partition has all clusters >10 cells",
            "required_gene_coverage_below_retired_guard": COPYKAT_RETIRED_GENE_GUARD,
            "required_internal_query_cells_below": 11,
            "scientific_parameters_changed": False,
            "v2_outputs_modified": False,
        },
        "interpretation": "not_estimable is missing method output, not a biological negative",
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print("PART1_CNV_V2_TERMINAL_RECOVERY_ALL_PASS infercnv=34 copykat_completed=38 copykat_not_estimable=3")


if __name__ == "__main__":
    main()
