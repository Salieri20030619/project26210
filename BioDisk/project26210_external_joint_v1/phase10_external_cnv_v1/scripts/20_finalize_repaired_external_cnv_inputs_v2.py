#!/usr/bin/env python3
"""Finalize already-written v2 repairs after validating type-scoped unit identity."""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-v1", required=True)
    parser.add_argument("--input-v2", required=True)
    args = parser.parse_args()
    v1 = Path(args.input_v1).resolve()
    v2 = Path(args.input_v2).resolve()
    if (v2 / "FULL_INPUT_BUILD_STATUS.json").exists() or (v2 / "RUN_COMPLETE").exists():
        raise FileExistsError("No-overwrite gate: v2 terminal output exists")

    rows: list[dict[str, object]] = []
    for unit_type, expected_count in (("infercnv", 49), ("copykat", 62)):
        v2_dirs = sorted(path for path in (v2 / unit_type).iterdir() if path.is_dir())
        if len(v2_dirs) != expected_count:
            raise RuntimeError(f"{unit_type} directory count {len(v2_dirs)} != {expected_count}")
        for target_dir in v2_dirs:
            source_dir = v1 / unit_type / target_dir.name
            source_mapping = pd.read_csv(source_dir / "cell_id_mapping_v1.tsv.gz", sep="\t", compression="gzip", dtype=str)
            target_mapping = pd.read_csv(target_dir / "cell_id_mapping_v1.tsv.gz", sep="\t", compression="gzip", dtype=str)
            metadata = json.loads((target_dir / "unit_metadata_v1.json").read_text(encoding="utf-8"))
            if not source_mapping["exact_cell_id"].equals(target_mapping["exact_cell_id"]):
                raise RuntimeError(f"Cell order changed: {unit_type}/{target_dir.name}")
            if not source_mapping["role"].equals(target_mapping["role"]):
                raise RuntimeError(f"Role order changed: {unit_type}/{target_dir.name}")
            if not target_mapping["exact_cell_id"].is_unique or not target_mapping["tool_id"].is_unique:
                raise RuntimeError(f"v2 identity uniqueness failed: {unit_type}/{target_dir.name}")
            query = target_mapping["role"].eq("query")
            reference = target_mapping["role"].eq("reference")
            expected_q = [f"Q{i:07d}" for i in range(1, int(query.sum()) + 1)]
            expected_r = [f"R{i:07d}" for i in range(1, int(reference.sum()) + 1)]
            if target_mapping.loc[query, "tool_id"].tolist() != expected_q:
                raise RuntimeError(f"Query tool ID sequence failed: {unit_type}/{target_dir.name}")
            if target_mapping.loc[reference, "tool_id"].tolist() != expected_r:
                raise RuntimeError(f"Reference tool ID sequence failed: {unit_type}/{target_dir.name}")
            for filename in ("counts_genes_by_cells_v1.mtx.gz", "gene_order_v1.tsv.gz"):
                source_path = source_dir / filename
                target_path = target_dir / filename
                if source_path.stat().st_size != target_path.stat().st_size:
                    raise RuntimeError(f"Byte-size mismatch: {unit_type}/{target_dir.name}/{filename}")
            if metadata.get("input_version") != "external_cnv_full_inputs_v2":
                raise RuntimeError(f"Missing v2 provenance: {unit_type}/{target_dir.name}")
            rows.append(
                {
                    "unit_type": unit_type,
                    "unit_id": str(metadata["unit_id"]),
                    "n_query": int(query.sum()),
                    "n_reference": int(reference.sum()),
                    "n_cells": len(target_mapping),
                    "status": "PASS",
                }
            )

    result = pd.DataFrame(rows)
    if len(result) != 111 or result.duplicated(["unit_type", "unit_id"]).any():
        raise RuntimeError("Type-scoped unit identity gate failed")
    incomplete = sorted(str(path) for path in v2.rglob("*") if path.is_file() and ".part" in path.name)
    if incomplete:
        raise RuntimeError(f"Incomplete v2 files: {len(incomplete)}")
    audit_path = v2 / "repair_finalization_audit_v2.tsv"
    audit_part = audit_path.with_name(audit_path.name + ".part")
    result.to_csv(audit_part, sep="\t", index=False)
    os.replace(audit_part, audit_path)

    source_build = json.loads((v1 / "FULL_INPUT_BUILD_STATUS.json").read_text(encoding="utf-8"))
    status = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_input_root": str(v1),
        "source_build_status": source_build,
        "n_units": len(result),
        "n_infercnv_units": int(result["unit_type"].eq("infercnv").sum()),
        "n_copykat_units": int(result["unit_type"].eq("copykat").sum()),
        "all_tool_ids_unique_within_unit": True,
        "all_exact_cell_and_role_order_identical_to_v1": True,
        "matrix_and_gene_file_sizes_identical_to_v1": True,
        "matrix_content_policy": "v2 copied v1 matrices and gene order byte-for-byte; only tool_id values and dependent mapping metadata changed",
        "v1_validation_failure": "v1 repeated Q0000001/R0000001 within each role",
        "v2_finalization_note": "unit_id uniqueness is scoped to (unit_type, unit_id), allowing the same biological unit to be processed by both methods",
    }
    status_part = v2 / "FULL_INPUT_BUILD_STATUS.part.json"
    status_part.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(status_part, v2 / "FULL_INPUT_BUILD_STATUS.json")
    run_part = v2 / "RUN_COMPLETE.part"
    run_part.write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")
    os.replace(run_part, v2 / "RUN_COMPLETE")
    print("EXTERNAL_CNV_INPUT_V2_FINALIZATION_PASS units=111")


if __name__ == "__main__":
    main()
