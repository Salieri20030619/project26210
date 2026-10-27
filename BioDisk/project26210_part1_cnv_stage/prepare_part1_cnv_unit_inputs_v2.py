#!/usr/bin/env python3
"""Prepare per-patient inferCNV and per-sample CopyKAT inputs with safe cell IDs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT_H5AD = PROJECT / (
    "results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/"
    "part1_scvi_integrated_cell_annotated_v1.h5ad"
)
GENE_ORDER = PROJECT / (
    "results/infercnv_rerun_seurat_v5/"
    "gene_order_EnsDb.Hsapiens.v86_seurat_genes.txt"
)
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
DEFAULT_PREFLIGHT = PROJECT / (
    "results/part1_rebuild_v1/03_malignancy_cnv/"
    "malignancy_cnv_input_audit_v2"
)
DEFAULT_OUTPUT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_part1_cnv_stage/cnv_unit_inputs_v1"
)

SEED = 26210
REFERENCE_TARGET = {"infercnv": 200, "copykat": 1000}
CLEAR_EPITHELIAL = "Epithelial cell"
LYMPHOID = ("T cell", "NK cell", "B cell")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=Path, default=INPUT_H5AD)
    parser.add_argument("--gene-order", type=Path, default=GENE_ORDER)
    parser.add_argument("--preflight-dir", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--unit-type", choices=("all", "infercnv", "copykat"), default="all")
    parser.add_argument("--unit-filter", default="", help="Comma-separated exact patient/sample IDs")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    return cleaned.strip("_") or "unit"


def stable_seed(unit_type: str, unit_id: str) -> int:
    digest = hashlib.sha256(f"{SEED}|{unit_type}|{unit_id}".encode()).digest()
    return SEED + int.from_bytes(digest[:4], "little") % 1_000_000_000


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def run_memory_gate(unit_type: str, unit_id: str) -> dict[str, str]:
    result = subprocess.run(
        [str(MEMORY_GATE)], check=False, text=True, capture_output=True
    )
    output = result.stdout.strip()
    print(output, flush=True)
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate blocked {unit_type} {unit_id}: {output} {result.stderr}")
    fields = output.split("\t")
    record = {"unit_type": unit_type, "unit_id": unit_id, "raw": output}
    for field in fields[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            record[key] = value
    record["timestamp"] = fields[0] if fields else ""
    return record


def balanced_reference(
    candidates: pd.DataFrame, target: int, seed: int
) -> pd.DataFrame:
    candidates = candidates.loc[~candidates.index.duplicated()].copy()
    if len(candidates) <= target:
        return candidates.sort_values("_row_index")
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    groups = list(candidates.groupby("sample_id", sort=True, observed=True))
    quota = max(1, target // len(groups))
    for _, group in groups:
        positions = group["_row_index"].to_numpy(dtype=np.int64)
        take = min(quota, len(positions))
        chosen = rng.choice(positions, size=take, replace=False)
        selected.extend(int(x) for x in chosen)
    selected = list(dict.fromkeys(selected))
    if len(selected) < target:
        remaining = candidates.loc[~candidates["_row_index"].isin(selected), "_row_index"].to_numpy(dtype=np.int64)
        take = min(target - len(selected), len(remaining))
        if take:
            selected.extend(int(x) for x in rng.choice(remaining, size=take, replace=False))
    selected = sorted(selected[:target])
    return candidates.loc[candidates["_row_index"].isin(selected)].sort_values("_row_index")


def choose_reference(obs: pd.DataFrame, unit_row: pd.Series, unit_type: str) -> pd.DataFrame:
    mode = str(unit_row["reference_mode"])
    patient = str(unit_row["curated_patient_id"])
    study = str(unit_row["study"])
    if mode == "matched_patient_normal_epithelial":
        mask = (
            (obs["curated_patient_id"] == patient)
            & (obs["source_tissue"] == "Normal")
            & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
        )
    elif mode == "same_study_normal_epithelial_pool":
        mask = (
            (obs["study"] == study)
            & (obs["source_tissue"] == "Normal")
            & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
        )
    elif mode == "matched_patient_tumor_lymphoid":
        mask = (
            (obs["curated_patient_id"] == patient)
            & (obs["source_tissue"] == "Tumor")
            & obs["cell_type_broad_v1"].isin(LYMPHOID)
        )
    elif mode == "balanced_cross_study_normal_epithelial_pool":
        mask = (
            (obs["source_tissue"] == "Normal")
            & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
        )
    else:
        raise RuntimeError(f"Unsupported reference mode: {mode}")
    candidates = obs.loc[mask].copy()
    if len(candidates) < 50:
        raise RuntimeError(f"Too few reference candidates for {unit_type}: {len(candidates)}")
    return balanced_reference(
        candidates,
        target=REFERENCE_TARGET[unit_type],
        seed=stable_seed(unit_type, str(unit_row["unit_id"])),
    )


def prepare_unit(
    adata: ad.AnnData,
    obs: pd.DataFrame,
    gene_order: pd.DataFrame,
    gene_indices: np.ndarray,
    unit_type: str,
    unit_row: pd.Series,
    output_root: Path,
    force: bool,
) -> dict[str, object]:
    unit_id = str(unit_row["unit_id"])
    unit_dir = output_root / unit_type / safe_name(unit_id)
    unit_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = unit_dir / "counts_genes_by_cells_v1.mtx.gz"
    mapping_path = unit_dir / "cell_id_mapping_v1.tsv.gz"
    gene_order_path = unit_dir / "gene_order_v1.tsv.gz"
    metadata_path = unit_dir / "unit_metadata_v1.json"

    if matrix_path.exists() and mapping_path.exists() and gene_order_path.exists() and metadata_path.exists() and not force:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("status") == "completed":
            return metadata

    if unit_type == "infercnv":
        query_mask = (
            (obs["curated_patient_id"] == str(unit_row["curated_patient_id"]))
            & (obs["source_tissue"] == "Tumor")
            & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
        )
    else:
        query_mask = (
            (obs["sample_id"] == unit_id)
            & (obs["source_tissue"] == "Tumor")
            & (obs["cell_type_broad_v1"] == CLEAR_EPITHELIAL)
        )
    query = obs.loc[query_mask].copy().sort_values("_row_index")
    reference = choose_reference(obs, unit_row, unit_type)
    if set(query.index).intersection(reference.index):
        raise RuntimeError(f"Query/reference overlap for {unit_type} {unit_id}")

    selected = pd.concat(
        [query.assign(role="query"), reference.assign(role="reference")],
        axis=0,
    ).sort_values("_row_index")
    role_counter = {"query": 0, "reference": 0}
    tool_ids: list[str] = []
    for role in selected["role"]:
        role_counter[role] += 1
        prefix = "Q" if role == "query" else "R"
        tool_ids.append(f"{prefix}{role_counter[role]:07d}")
    selected.insert(0, "tool_id", tool_ids)
    selected.insert(1, "exact_cell_id", selected.index.astype(str))

    if len(query) != int(unit_row["n_tumor_epithelial"]):
        raise RuntimeError(
            f"Query count mismatch for {unit_type} {unit_id}: "
            f"{len(query)} != {unit_row['n_tumor_epithelial']}"
        )

    row_indices = selected["_row_index"].to_numpy(dtype=np.int64)
    block = adata.layers["counts"][row_indices, :]
    if not sp.issparse(block):
        block = sp.csr_matrix(block)
    else:
        block = block.tocsr()
    block = block[:, gene_indices]
    matrix = block.transpose().tocsc()
    if matrix.shape != (len(gene_order), len(selected)):
        raise RuntimeError(f"Unexpected matrix shape for {unit_type} {unit_id}: {matrix.shape}")
    if matrix.nnz == 0 or np.any(~np.isfinite(matrix.data)) or np.any(matrix.data < 0):
        raise RuntimeError(f"Invalid count matrix for {unit_type} {unit_id}")
    if np.any(np.abs(matrix.data - np.rint(matrix.data)) > 1e-6):
        raise RuntimeError(f"Noninteger count values for {unit_type} {unit_id}")

    temporary = matrix_path.with_suffix(".mtx.gz.tmp")
    with gzip.open(temporary, "wb", compresslevel=6) as handle:
        mmwrite(handle, matrix.astype(np.int64), field="integer", symmetry="general")
    temporary.replace(matrix_path)

    mapping_columns = [
        "tool_id",
        "exact_cell_id",
        "role",
        "sample_id",
        "curated_patient_id",
        "study",
        "source_tissue",
        "cell_type_broad_v1",
        "_row_index",
    ]
    selected[mapping_columns].to_csv(
        mapping_path, sep="\t", index=False, compression="gzip", lineterminator="\n"
    )
    gene_order.to_csv(
        gene_order_path, sep="\t", index=False, compression="gzip", lineterminator="\n"
    )
    metadata: dict[str, object] = {
        "schema_version": "part1_cnv_unit_metadata_v1",
        "status": "completed",
        "completed_at": datetime.now().astimezone().isoformat(),
        "unit_type": unit_type,
        "unit_id": unit_id,
        "safe_unit_id": safe_name(unit_id),
        "curated_patient_id": str(unit_row["curated_patient_id"]),
        "study": str(unit_row["study"]),
        "reference_mode": str(unit_row["reference_mode"]),
        "seed": stable_seed(unit_type, unit_id),
        "n_query_cells": int(len(query)),
        "n_reference_cells": int(len(reference)),
        "n_input_cells": int(len(selected)),
        "n_genes": int(len(gene_order)),
        "matrix_nnz": int(matrix.nnz),
        "matrix_mtx": str(matrix_path),
        "matrix_mtx_size_bytes": matrix_path.stat().st_size,
        "matrix_mtx_sha256": sha256_file(matrix_path),
        "mapping_tsv": str(mapping_path),
        "mapping_tsv_size_bytes": mapping_path.stat().st_size,
        "mapping_tsv_sha256": sha256_file(mapping_path),
        "gene_order_tsv": str(gene_order_path),
        "gene_order_tsv_size_bytes": gene_order_path.stat().st_size,
        "gene_order_tsv_sha256": sha256_file(gene_order_path),
        "exact_cell_id_policy": "all external IDs preserved in mapping; only safe surrogate IDs passed to R tools",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return metadata


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    requested = {x.strip() for x in args.unit_filter.split(",") if x.strip()}

    patient_manifest = pd.read_csv(
        args.preflight_dir / "infercnv_patient_manifest_preflight_v2.tsv", sep="\t"
    )
    patient_manifest = patient_manifest[patient_manifest["eligible_infercnv"]].copy()
    patient_manifest["unit_id"] = patient_manifest["curated_patient_id"]
    sample_manifest = pd.read_csv(
        args.preflight_dir / "copykat_sample_manifest_preflight_v2.tsv", sep="\t"
    )
    sample_manifest = sample_manifest[sample_manifest["eligible_copykat"]].copy()
    sample_manifest["unit_id"] = sample_manifest["sample_id"]

    if requested:
        patient_manifest = patient_manifest[patient_manifest["unit_id"].isin(requested)]
        sample_manifest = sample_manifest[sample_manifest["unit_id"].isin(requested)]
    if args.unit_type == "infercnv":
        sample_manifest = sample_manifest.iloc[0:0]
    elif args.unit_type == "copykat":
        patient_manifest = patient_manifest.iloc[0:0]
    if requested and len(patient_manifest) + len(sample_manifest) == 0:
        raise RuntimeError(f"No eligible units matched filter: {sorted(requested)}")

    print(f"[{datetime.now().astimezone().isoformat()}] open backed H5AD", flush=True)
    adata = ad.read_h5ad(args.input_h5ad, backed="r")
    obs = adata.obs.copy()
    obs.index = adata.obs_names.astype(str)
    obs["_row_index"] = np.arange(adata.n_obs, dtype=np.int64)
    for column in (
        "sample_id",
        "curated_patient_id",
        "study",
        "source_tissue",
        "cell_type_broad_v1",
    ):
        obs[column] = obs[column].astype(str)

    order = pd.read_csv(
        args.gene_order,
        sep="\t",
        header=None,
        names=["gene", "chromosome", "start", "end"],
    )
    var_index = pd.Series(np.arange(adata.n_vars, dtype=np.int64), index=adata.var_names.astype(str))
    order = order[~order["gene"].duplicated()].copy()
    order = order[order["gene"].isin(var_index.index)].copy().reset_index(drop=True)
    gene_indices = var_index.loc[order["gene"]].to_numpy(dtype=np.int64)
    if len(order) != 25263:
        raise RuntimeError(f"Unexpected gene-order overlap: {len(order)} != 25263")

    records: list[dict[str, object]] = []
    gates: list[dict[str, str]] = []
    for unit_type, manifest in (("infercnv", patient_manifest), ("copykat", sample_manifest)):
        for _, unit_row in manifest.iterrows():
            unit_id = str(unit_row["unit_id"])
            print(f"[{datetime.now().astimezone().isoformat()}] prepare {unit_type} {unit_id}", flush=True)
            gates.append(run_memory_gate(unit_type, unit_id))
            records.append(
                prepare_unit(
                    adata,
                    obs,
                    order,
                    gene_indices,
                    unit_type,
                    unit_row,
                    args.output_dir,
                    args.force,
                )
            )
            pd.DataFrame(records).to_csv(
                args.output_dir / "cnv_unit_input_manifest_v1.tsv",
                sep="\t",
                index=False,
                lineterminator="\n",
            )
            pd.DataFrame(gates).to_csv(
                args.output_dir / "memory_gate_log_v1.tsv",
                sep="\t",
                index=False,
                lineterminator="\n",
            )

    metadata = {
        "schema_version": "part1_cnv_unit_input_preparation_v2",
        "status": "completed",
        "completed_at": datetime.now().astimezone().isoformat(),
        "input_h5ad": str(args.input_h5ad),
        "gene_order": str(args.gene_order),
        "output_dir": str(args.output_dir),
        "unit_type_requested": args.unit_type,
        "unit_filter": sorted(requested),
        "n_infercnv_units": int(sum(x["unit_type"] == "infercnv" for x in records)),
        "n_copykat_units": int(sum(x["unit_type"] == "copykat" for x in records)),
        "n_memory_gates": len(gates),
        "seed": SEED,
        "gene_count": len(order),
        "reference_targets": REFERENCE_TARGET,
    }
    (args.output_dir / "run_metadata_v1.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (args.output_dir / "RUN_COMPLETE").write_text(
        datetime.now().astimezone().isoformat() + "\n", encoding="utf-8"
    )
    print(
        "PART1_CNV_UNIT_INPUT_PREPARATION_V2_PASS "
        f"infercnv={metadata['n_infercnv_units']} copykat={metadata['n_copykat_units']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
