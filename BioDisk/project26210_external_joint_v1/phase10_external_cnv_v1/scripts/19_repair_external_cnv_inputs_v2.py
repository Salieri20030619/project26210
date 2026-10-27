#!/usr/bin/env python3
"""Create non-overwriting v2 inputs with unique R-tool IDs; matrices/cell order stay exact."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_copy(source: Path, final: Path) -> None:
    part = final.with_name(final.name + ".part")
    shutil.copy2(source, part)
    os.replace(part, final)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-v1", required=True)
    parser.add_argument("--output-v2", required=True)
    args = parser.parse_args()
    input_root = Path(args.input_v1).resolve()
    output_root = Path(args.output_v2).resolve()
    if output_root.exists():
        raise FileExistsError(f"No-overwrite gate: {output_root}")
    if not (input_root / "RUN_COMPLETE").exists():
        raise RuntimeError("v1 input build is incomplete")
    output_root.mkdir(parents=True)

    records: list[dict[str, object]] = []
    for unit_type in ("infercnv", "copykat"):
        unit_dirs = sorted(path for path in (input_root / unit_type).iterdir() if path.is_dir())
        for index, source_dir in enumerate(unit_dirs, start=1):
            target_dir = output_root / unit_type / source_dir.name
            target_dir.mkdir(parents=True)
            source_matrix = source_dir / "counts_genes_by_cells_v1.mtx.gz"
            source_mapping = source_dir / "cell_id_mapping_v1.tsv.gz"
            source_gene = source_dir / "gene_order_v1.tsv.gz"
            source_metadata = source_dir / "unit_metadata_v1.json"
            for path in (source_matrix, source_mapping, source_gene, source_metadata):
                if not path.is_file() or path.stat().st_size <= 0:
                    raise RuntimeError(f"Missing source input: {path}")

            matrix_path = target_dir / source_matrix.name
            gene_path = target_dir / source_gene.name
            mapping_path = target_dir / source_mapping.name
            metadata_path = target_dir / source_metadata.name
            atomic_copy(source_matrix, matrix_path)
            atomic_copy(source_gene, gene_path)

            mapping = pd.read_csv(source_mapping, sep="\t", compression="gzip", dtype=str)
            if mapping["exact_cell_id"].duplicated().any():
                raise RuntimeError(f"Duplicate exact cell IDs: {source_dir}")
            query_mask = mapping["role"].eq("query")
            reference_mask = mapping["role"].eq("reference")
            if not (query_mask | reference_mask).all():
                raise RuntimeError(f"Unexpected roles: {source_dir}")
            mapping.loc[query_mask, "tool_id"] = [f"Q{i:07d}" for i in range(1, int(query_mask.sum()) + 1)]
            mapping.loc[reference_mask, "tool_id"] = [f"R{i:07d}" for i in range(1, int(reference_mask.sum()) + 1)]
            if not mapping["tool_id"].is_unique:
                raise RuntimeError(f"Repaired tool IDs remain duplicated: {source_dir}")
            mapping_part = mapping_path.with_name(mapping_path.name + ".part")
            mapping.to_csv(mapping_part, sep="\t", index=False, compression="gzip", lineterminator="\n")
            os.replace(mapping_part, mapping_path)

            metadata = json.loads(source_metadata.read_text(encoding="utf-8"))
            original_mapping_sha = metadata.get("mapping_tsv_sha256")
            metadata.update(
                matrix_mtx=str(matrix_path),
                matrix_mtx_size_bytes=matrix_path.stat().st_size,
                matrix_mtx_sha256=sha256(matrix_path),
                mapping_tsv=str(mapping_path),
                mapping_tsv_size_bytes=mapping_path.stat().st_size,
                mapping_tsv_sha256=sha256(mapping_path),
                gene_order_tsv=str(gene_path),
                gene_order_tsv_size_bytes=gene_path.stat().st_size,
                gene_order_tsv_sha256=sha256(gene_path),
                input_version="external_cnv_full_inputs_v2",
                repair_provenance={
                    "created_at_utc": datetime.now(timezone.utc).isoformat(),
                    "source_unit_dir": str(source_dir),
                    "source_mapping_sha256": original_mapping_sha,
                    "reason": "v1 assigned one repeated tool ID per role; v2 assigns deterministic sequential Q/R IDs without changing row/cell/matrix-column order",
                },
            )
            metadata_part = metadata_path.with_name(metadata_path.name + ".part")
            metadata_part.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            os.replace(metadata_part, metadata_path)
            records.append(
                {
                    "unit_type": unit_type,
                    "unit_id": metadata["unit_id"],
                    "n_query": int(query_mask.sum()),
                    "n_reference": int(reference_mask.sum()),
                    "n_cells": len(mapping),
                    "n_genes": int(metadata["n_genes"]),
                    "matrix_bytes": matrix_path.stat().st_size,
                    "mapping_sha256": metadata["mapping_tsv_sha256"],
                    "status": "PASS",
                }
            )
            pd.DataFrame(records).to_csv(output_root / "repair_progress_v2.tsv.part", sep="\t", index=False)
            os.replace(output_root / "repair_progress_v2.tsv.part", output_root / "repair_progress_v2.tsv")
            print(f"REPAIRED {unit_type} {index}/{len(unit_dirs)} {source_dir.name}", flush=True)

    result = pd.DataFrame(records)
    if len(result) != 111 or result["unit_id"].duplicated().any():
        raise RuntimeError("Repaired unit count/uniqueness gate failed")
    build_v1 = json.loads((input_root / "FULL_INPUT_BUILD_STATUS.json").read_text(encoding="utf-8"))
    status = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_input_root": str(input_root),
        "source_build_status": build_v1,
        "n_units": len(result),
        "n_infercnv_units": int(result["unit_type"].eq("infercnv").sum()),
        "n_copykat_units": int(result["unit_type"].eq("copykat").sum()),
        "all_tool_ids_unique_within_unit": True,
        "matrix_content_policy": "byte-identical copies of v1 matrices and gene-order files; mapping row order and exact cell IDs unchanged",
    }
    status_part = output_root / "FULL_INPUT_BUILD_STATUS.part.json"
    status_part.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(status_part, output_root / "FULL_INPUT_BUILD_STATUS.json")
    run_part = output_root / "RUN_COMPLETE.part"
    run_part.write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")
    os.replace(run_part, output_root / "RUN_COMPLETE")
    print("EXTERNAL_CNV_INPUT_REPAIR_V2_PASS units=111", flush=True)


if __name__ == "__main__":
    main()
