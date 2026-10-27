#!/usr/bin/env python3
"""Validate every external inferCNV/CopyKAT input before cohort launch."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "unit"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def mtx_shape(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="ascii") as handle:
        first = handle.readline().strip()
        if not first.startswith("%%MatrixMarket matrix coordinate"):
            raise RuntimeError(f"Unexpected MatrixMarket header: {path}")
        for line in handle:
            if not line.startswith("%"):
                rows, cols, nnz = (int(x) for x in line.split())
                return rows, cols, nnz
    raise RuntimeError(f"Missing MatrixMarket dimensions: {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--preflight-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    input_root = Path(args.input_root).resolve()
    preflight_dir = Path(args.preflight_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    if out_dir.exists():
        raise FileExistsError(f"No-overwrite gate: {out_dir}")
    out_dir.mkdir(parents=True)

    if not (input_root / "RUN_COMPLETE").exists():
        raise RuntimeError("Full input build has no RUN_COMPLETE")
    build_status = json.loads((input_root / "FULL_INPUT_BUILD_STATUS.json").read_text(encoding="utf-8"))
    if build_status.get("status") != "PASS" or int(build_status.get("n_units", -1)) != 111:
        raise RuntimeError("Full input build status/count gate failed")

    patient = pd.read_csv(preflight_dir / "external_infercnv_patient_preflight_v1.tsv", sep="\t")
    sample = pd.read_csv(preflight_dir / "external_copykat_sample_preflight_v1.tsv", sep="\t")
    patient = patient.loc[patient["eligible_infercnv"].astype(bool)].copy()
    sample = sample.loc[sample["eligible_copykat"].astype(bool)].copy()
    expected = [("infercnv", row) for _, row in patient.iterrows()] + [("copykat", row) for _, row in sample.iterrows()]
    if len(expected) != 111:
        raise RuntimeError(f"Expected 111 eligible units, found {len(expected)}")

    rows: list[dict[str, object]] = []
    failures: list[str] = []
    for unit_type, row in expected:
        unit_id = str(row["unit_id"])
        unit_dir = input_root / unit_type / safe_name(unit_id)
        metadata_path = unit_dir / "unit_metadata_v1.json"
        matrix_path = unit_dir / "counts_genes_by_cells_v1.mtx.gz"
        mapping_path = unit_dir / "cell_id_mapping_v1.tsv.gz"
        gene_path = unit_dir / "gene_order_v1.tsv.gz"
        issues: list[str] = []
        for path in (metadata_path, matrix_path, mapping_path, gene_path):
            if not path.is_file() or path.stat().st_size <= 0:
                issues.append(f"missing_or_empty:{path.name}")
        record: dict[str, object] = {"unit_type": unit_type, "unit_id": unit_id}
        if not issues:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            mapping = pd.read_csv(mapping_path, sep="\t", compression="gzip", dtype=str)
            genes = pd.read_csv(gene_path, sep="\t", compression="gzip")
            m_rows, m_cols, m_nnz = mtx_shape(matrix_path)
            n_query = int(mapping["role"].eq("query").sum())
            n_reference = int(mapping["role"].eq("reference").sum())
            expected_query = int(row["n_tumor_epithelial"])
            checks = {
                "metadata_status": metadata.get("status") == "completed",
                "metadata_unit_type": metadata.get("unit_type") == unit_type,
                "metadata_unit_id": metadata.get("unit_id") == unit_id,
                "query_count": n_query == expected_query == int(metadata.get("n_query_cells", -1)),
                "reference_count": n_reference >= 50 and n_reference == int(metadata.get("n_reference_cells", -1)),
                "mapping_unique_exact_id": mapping["exact_cell_id"].is_unique,
                "mapping_unique_tool_id": mapping["tool_id"].is_unique,
                "gene_count": len(genes) >= 10000 and len(genes) == int(metadata.get("n_genes", -1)),
                "gene_unique": genes["gene"].is_unique,
                "matrix_shape": m_rows == len(genes) and m_cols == len(mapping),
                "matrix_nnz": m_nnz == int(metadata.get("matrix_nnz", -1)),
                "matrix_size": matrix_path.stat().st_size == int(metadata.get("matrix_mtx_size_bytes", -1)),
                "mapping_size": mapping_path.stat().st_size == int(metadata.get("mapping_tsv_size_bytes", -1)),
                "gene_size": gene_path.stat().st_size == int(metadata.get("gene_order_tsv_size_bytes", -1)),
                "matrix_sha256": sha256(matrix_path) == metadata.get("matrix_mtx_sha256"),
                "mapping_sha256": sha256(mapping_path) == metadata.get("mapping_tsv_sha256"),
                "gene_sha256": sha256(gene_path) == metadata.get("gene_order_tsv_sha256"),
            }
            issues.extend(name for name, passed in checks.items() if not passed)
            record.update(
                n_query=n_query,
                n_reference=n_reference,
                n_genes=len(genes),
                matrix_nnz=m_nnz,
                matrix_bytes=matrix_path.stat().st_size,
            )
        record["status"] = "PASS" if not issues else "FAIL"
        record["issues"] = ";".join(issues)
        rows.append(record)
        if issues:
            failures.append(f"{unit_type}:{unit_id}:{','.join(issues)}")

    incomplete = sorted(str(path) for path in input_root.rglob("*") if path.is_file() and ".part" in path.name)
    if incomplete:
        failures.append(f"incomplete_files:{len(incomplete)}")
    result = pd.DataFrame(rows)
    result_path = out_dir / "external_cnv_input_validation_v1.tsv"
    result.to_csv(result_path, sep="\t", index=False)
    status = "PASS" if not failures and len(result) == 111 and result["status"].eq("PASS").all() else "FAIL"
    payload = {
        "status": status,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(input_root),
        "expected_units": 111,
        "validated_units": int(len(result)),
        "passed_units": int(result["status"].eq("PASS").sum()),
        "infercnv_units": int(result["unit_type"].eq("infercnv").sum()),
        "copykat_units": int(result["unit_type"].eq("copykat").sum()),
        "incomplete_files": incomplete,
        "failures": failures,
        "validation_table_sha256": sha256(result_path),
    }
    status_path = out_dir / "FULL_INPUT_VALIDATION_STATUS.json"
    status_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if status == "PASS":
        (out_dir / "RUN_COMPLETE").write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
