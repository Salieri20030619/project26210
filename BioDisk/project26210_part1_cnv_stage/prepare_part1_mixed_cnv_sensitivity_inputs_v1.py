#!/usr/bin/env python3
"""Append the two tumor mixed-lineage cells to their frozen primary CNV inputs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread, mmwrite


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1/03_malignancy_cnv"
INPUT_H5AD = PROJECT / (
    "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/"
    "part1_scvi_integrated_cell_annotated_v1.h5ad"
)
PRIMARY_INPUT_ROOT = ROOT / "cnv_unit_inputs_v1"
DEFAULT_OUTPUT_ROOT = ROOT / "mixed_cnv_sensitivity_v1/inputs"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
MIXED_LABEL = "Epithelial/Myeloid mixed"
EXPECTED = {
    ("infercnv", "GSE156632__RCC3"): "RCC3t_GAGGTGATCAGGCGAA.1",
    ("infercnv", "GSE156632__RCC5"): "RCC5t_ACTTGTTGTTCAGCGC.1",
    ("copykat", "RCC3t"): "RCC3t_GAGGTGATCAGGCGAA.1",
    ("copykat", "RCC5t"): "RCC5t_ACTTGTTGTTCAGCGC.1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=Path, default=INPUT_H5AD)
    parser.add_argument("--primary-input-root", type=Path, default=PRIMARY_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def memory_gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append(
        {"stage": stage, "exit_status": str(result.returncode), "gate_output": output}
    )
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def read_matrix(path: Path) -> sp.csc_matrix:
    with gzip.open(path, "rb") as handle:
        matrix = mmread(handle)
    return sp.csc_matrix(matrix)


def append_unit(
    adata: ad.AnnData,
    obs: pd.DataFrame,
    unit_type: str,
    unit_id: str,
    mixed_cell_id: str,
    primary_root: Path,
    output_root: Path,
) -> dict[str, object]:
    source_dir = primary_root / unit_type / unit_id
    target_dir = output_root / unit_type / unit_id
    source_matrix_path = source_dir / "counts_genes_by_cells_v1.mtx.gz"
    source_mapping_path = source_dir / "cell_id_mapping_v1.tsv.gz"
    source_gene_path = source_dir / "gene_order_v1.tsv.gz"
    source_metadata_path = source_dir / "unit_metadata_v1.json"
    required = [
        source_matrix_path,
        source_mapping_path,
        source_gene_path,
        source_metadata_path,
    ]
    if not all(path.is_file() and path.stat().st_size > 0 for path in required):
        raise RuntimeError(f"Primary inputs incomplete for {unit_type} {unit_id}")
    if target_dir.exists() and any(target_dir.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty sensitivity input: {target_dir}")
    target_dir.mkdir(parents=True, exist_ok=True)

    source_meta = json.loads(source_metadata_path.read_text(encoding="utf-8"))
    if source_meta.get("status") != "completed" or source_meta.get("unit_type") != unit_type:
        raise RuntimeError(f"Unexpected primary metadata for {unit_type} {unit_id}")
    mapping = pd.read_csv(source_mapping_path, sep="\t", dtype=str)
    gene_order = pd.read_csv(source_gene_path, sep="\t", dtype=str)
    matrix = read_matrix(source_matrix_path)
    if matrix.shape != (len(gene_order), len(mapping)):
        raise RuntimeError(f"Primary dimensions mismatch for {unit_type} {unit_id}")
    if mapping["exact_cell_id"].duplicated().any() or mixed_cell_id in set(mapping["exact_cell_id"]):
        raise RuntimeError(f"Primary mapping is invalid for {unit_type} {unit_id}")

    if mixed_cell_id not in obs.index:
        raise RuntimeError(f"Mixed cell absent from H5AD: {mixed_cell_id}")
    mixed = obs.loc[mixed_cell_id]
    if str(mixed["source_tissue"]) != "Tumor" or str(mixed["cell_type_broad_v1"]) != MIXED_LABEL:
        raise RuntimeError(f"Unexpected tissue/label for {mixed_cell_id}")
    if str(mixed["curated_patient_id"]) != str(source_meta["curated_patient_id"]):
        raise RuntimeError(f"Patient mismatch for {unit_type} {unit_id}")
    if unit_type == "copykat" and str(mixed["sample_id"]) != unit_id:
        raise RuntimeError(f"Sample mismatch for {unit_type} {unit_id}")

    var_lookup = pd.Series(np.arange(adata.n_vars, dtype=np.int64), index=adata.var_names.astype(str))
    if gene_order["gene"].duplicated().any() or not gene_order["gene"].isin(var_lookup.index).all():
        raise RuntimeError(f"Gene order cannot be mapped exactly for {unit_type} {unit_id}")
    row_index = int(obs.index.get_loc(mixed_cell_id))
    gene_indices = var_lookup.loc[gene_order["gene"]].to_numpy(dtype=np.int64)
    vector = adata.layers["counts"][row_index, gene_indices]
    if not sp.issparse(vector):
        vector = sp.csc_matrix(np.asarray(vector).reshape(1, -1))
    vector = sp.csc_matrix(vector).transpose().tocsc()
    if vector.shape != (len(gene_order), 1):
        raise RuntimeError(f"Unexpected mixed count vector shape for {unit_type} {unit_id}")
    if vector.nnz == 0 or np.any(~np.isfinite(vector.data)) or np.any(vector.data < 0):
        raise RuntimeError(f"Invalid mixed count vector for {unit_type} {unit_id}")
    if np.any(np.abs(vector.data - np.rint(vector.data)) > 1e-6):
        raise RuntimeError(f"Noninteger mixed count vector for {unit_type} {unit_id}")

    q_numbers = mapping.loc[mapping["tool_id"].str.startswith("Q"), "tool_id"].str[1:].astype(int)
    new_tool_id = f"Q{int(q_numbers.max()) + 1:07d}"
    new_row = {
        "tool_id": new_tool_id,
        "exact_cell_id": mixed_cell_id,
        "role": "query",
        "sample_id": str(mixed["sample_id"]),
        "curated_patient_id": str(mixed["curated_patient_id"]),
        "study": str(mixed["study"]),
        "source_tissue": str(mixed["source_tissue"]),
        "cell_type_broad_v1": str(mixed["cell_type_broad_v1"]),
        "_row_index": str(row_index),
    }
    if set(mapping.columns) != set(new_row):
        raise RuntimeError(f"Unexpected primary mapping columns for {unit_type} {unit_id}")
    new_mapping = pd.concat(
        [mapping, pd.DataFrame([{column: new_row[column] for column in mapping.columns}])],
        ignore_index=True,
    )
    new_matrix = sp.hstack([matrix, vector], format="csc")
    if new_matrix.shape != (matrix.shape[0], matrix.shape[1] + 1):
        raise RuntimeError(f"Append failed for {unit_type} {unit_id}")
    if (new_matrix[:, :-1] != matrix).nnz != 0 or not new_mapping.iloc[:-1].equals(mapping):
        raise RuntimeError(f"Primary input preservation failed for {unit_type} {unit_id}")

    matrix_path = target_dir / "counts_genes_by_cells_v1.mtx.gz"
    mapping_path = target_dir / "cell_id_mapping_v1.tsv.gz"
    gene_path = target_dir / "gene_order_v1.tsv.gz"
    metadata_path = target_dir / "unit_metadata_v1.json"
    with gzip.open(matrix_path, "wb", compresslevel=6) as handle:
        mmwrite(handle, new_matrix.astype(np.int64), field="integer", symmetry="general")
    new_mapping.to_csv(
        mapping_path, sep="\t", index=False, compression="gzip", lineterminator="\n"
    )
    gene_order.to_csv(
        gene_path, sep="\t", index=False, compression="gzip", lineterminator="\n"
    )
    metadata = dict(source_meta)
    metadata.update(
        {
            "completed_at": datetime.now().astimezone().isoformat(),
            "n_query_cells": int(source_meta["n_query_cells"]) + 1,
            "n_input_cells": int(source_meta["n_input_cells"]) + 1,
            "matrix_nnz": int(new_matrix.nnz),
            "matrix_mtx": str(matrix_path.resolve()),
            "matrix_mtx_size_bytes": matrix_path.stat().st_size,
            "matrix_mtx_sha256": sha256(matrix_path),
            "mapping_tsv": str(mapping_path.resolve()),
            "mapping_tsv_size_bytes": mapping_path.stat().st_size,
            "mapping_tsv_sha256": sha256(mapping_path),
            "gene_order_tsv": str(gene_path.resolve()),
            "gene_order_tsv_size_bytes": gene_path.stat().st_size,
            "gene_order_tsv_sha256": sha256(gene_path),
            "sensitivity_design": "frozen_primary_input_plus_one_tumor_mixed_query",
            "sensitivity_mixed_cell_id": mixed_cell_id,
            "sensitivity_mixed_tool_id": new_tool_id,
            "source_primary_input_dir": str(source_dir.resolve()),
            "source_primary_matrix_sha256": sha256(source_matrix_path),
            "source_primary_mapping_sha256": sha256(source_mapping_path),
            "source_primary_gene_order_sha256": sha256(source_gene_path),
            "primary_columns_preserved_exactly": True,
            "primary_mapping_rows_preserved_exactly": True,
        }
    )
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return {
        "unit_type": unit_type,
        "unit_id": unit_id,
        "mixed_cell_id": mixed_cell_id,
        "mixed_tool_id": new_tool_id,
        "n_primary_query": int(source_meta["n_query_cells"]),
        "n_sensitivity_query": int(metadata["n_query_cells"]),
        "n_reference": int(metadata["n_reference_cells"]),
        "n_genes": int(metadata["n_genes"]),
        "primary_matrix_sha256": metadata["source_primary_matrix_sha256"],
        "sensitivity_matrix_sha256": metadata["matrix_mtx_sha256"],
        "primary_columns_preserved_exactly": True,
        "primary_mapping_rows_preserved_exactly": True,
    }


def main() -> None:
    args = parse_args()
    if args.output_root.resolve() == args.primary_input_root.resolve():
        raise RuntimeError("Sensitivity output root must differ from primary input root")
    args.output_root.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    memory_gate("read_annotated_h5ad", gates)
    adata = ad.read_h5ad(args.input_h5ad, backed="r")
    obs = adata.obs.copy()
    required = {
        "sample_id", "curated_patient_id", "study", "source_tissue",
        "cell_type_broad_v1", "total_counts_raw", "n_genes_raw",
    }
    if missing := sorted(required.difference(obs.columns)):
        raise RuntimeError(f"Annotated H5AD missing columns: {missing}")
    mixed = obs.loc[obs["cell_type_broad_v1"].astype(str).eq(MIXED_LABEL)].copy()
    tumor_mixed = mixed.loc[mixed["source_tissue"].astype(str).eq("Tumor")].copy()
    normal_mixed = mixed.loc[mixed["source_tissue"].astype(str).eq("Normal")].copy()
    if len(mixed) != 918 or len(tumor_mixed) != 2 or len(normal_mixed) != 916:
        raise RuntimeError(
            f"Unexpected mixed registry: total={len(mixed)} tumor={len(tumor_mixed)} normal={len(normal_mixed)}"
        )
    if set(tumor_mixed.index.astype(str)) != set(EXPECTED.values()):
        raise RuntimeError("Tumor mixed cell IDs differ from frozen expectation")

    registry = mixed[
        [
            "sample_id", "curated_patient_id", "study", "source_tissue",
            "cell_type_broad_v1", "total_counts_raw", "n_genes_raw",
        ]
    ].copy()
    registry.insert(0, "exact_cell_id", registry.index.astype(str))
    registry["sensitivity_scope_v1"] = np.where(
        registry["source_tissue"].astype(str).eq("Tumor"),
        "tumor_query_rerun",
        "normal_tissue_not_a_tumor_query",
    )
    registry.to_csv(
        args.output_root / "mixed_cell_registry_v1.tsv", sep="\t", index=False, lineterminator="\n"
    )

    results: list[dict[str, object]] = []
    for unit_type, unit_id in [
        ("infercnv", "GSE156632__RCC3"),
        ("infercnv", "GSE156632__RCC5"),
        ("copykat", "RCC3t"),
        ("copykat", "RCC5t"),
    ]:
        memory_gate(f"prepare_{unit_type}_{unit_id}", gates)
        results.append(
            append_unit(
                adata,
                obs,
                unit_type,
                unit_id,
                EXPECTED[(unit_type, unit_id)],
                args.primary_input_root,
                args.output_root,
            )
        )
    adata.file.close()
    pd.DataFrame(results).to_csv(
        args.output_root / "input_comparison_v1.tsv", sep="\t", index=False, lineterminator="\n"
    )
    pd.DataFrame(gates).to_csv(
        args.output_root / "memory_gate_log_v1.tsv", sep="\t", index=False, lineterminator="\n"
    )
    metadata = {
        "schema_version": "part1_mixed_cnv_sensitivity_inputs_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "input_h5ad": str(args.input_h5ad.resolve()),
        "input_h5ad_sha256": sha256(args.input_h5ad),
        "primary_input_root": str(args.primary_input_root.resolve()),
        "output_root": str(args.output_root.resolve()),
        "n_mixed_total": 918,
        "n_mixed_tumor": 2,
        "n_mixed_normal": 916,
        "n_sensitivity_units": 4,
        "design": "append one tumor mixed cell to each affected frozen primary unit",
        "normal_mixed_policy": "not run as tumor queries",
    }
    (args.output_root / "run_metadata_v1.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_root / "RUN_COMPLETE").write_text(
        datetime.now().astimezone().isoformat() + "\n", encoding="utf-8"
    )
    print("PART1_MIXED_CNV_SENSITIVITY_INPUTS_V1_PASS")


if __name__ == "__main__":
    main()
