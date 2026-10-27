#!/usr/bin/env python3
"""Prepare frozen raw-count and metadata inputs for rebuilt Part 1 SCENIC/metaVIPER."""

from __future__ import annotations

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
from scipy.io import mmwrite


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
ROOT = PROJECT / "results/part1_rebuild_v1"
H5AD = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
CNMF_INPUT = ROOT / "05_cnmf/input_v1"
CNMF_FULL = ROOT / "05_cnmf/full_ranks4_8_v1"
CLONES = ROOT / "04_subclones_tymp/cnv_subclones_v1/formal_cell_cnv_clone_assignments_v1.tsv.gz"
OUT = ROOT / "07_pyscenic_v1/input"
GATE = PROJECT / "check_memory_gate.sh"

EXPECTED_CELLS = 11950
EXPECTED_GENES = 37795
EXPECTED_NNZ = 35713877
SEED = 26210


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(GATE), stage], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty input directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []

    for required in (
        H5AD,
        CNMF_INPUT / "RUN_COMPLETE",
        CNMF_FULL / "RUN_COMPLETE",
        CNMF_INPUT / "cnmf_cell_metadata_v1.tsv.gz",
        CNMF_FULL / "consensus_usage_rank7_v1.tsv.gz",
        CNMF_INPUT / "cnmf_hvg_genes_v1.tsv",
        CLONES,
    ):
        if not required.exists():
            raise FileNotFoundError(required)

    gate("scenic_input_open_h5ad_v1", gates)
    adata = ad.read_h5ad(H5AD, backed="r")
    selected = adata.obs["eligible_for_subclone_selection_v1"].astype(bool).to_numpy()
    idx = np.flatnonzero(selected)
    cells = adata.obs_names[idx].astype(str)
    genes = adata.var_names.astype(str).to_numpy()
    if len(idx) != EXPECTED_CELLS or len(genes) != EXPECTED_GENES:
        raise RuntimeError(f"Unexpected frozen input dimensions: {len(idx)} x {len(genes)}")

    gate("scenic_input_read_selected_counts_v1", gates)
    counts = sp.csr_matrix(adata.layers["counts"][idx, :], dtype=np.float64)
    adata.file.close()
    if counts.shape != (EXPECTED_CELLS, EXPECTED_GENES) or counts.nnz != EXPECTED_NNZ:
        raise RuntimeError(f"Unexpected counts shape/nnz: {counts.shape}, {counts.nnz}")
    if (
        np.any(~np.isfinite(counts.data))
        or np.any(counts.data < 0)
        or np.any(np.abs(counts.data - np.rint(counts.data)) > 1e-8)
    ):
        raise RuntimeError("Counts are not finite non-negative integers")

    base_meta = pd.read_csv(
        CNMF_INPUT / "cnmf_cell_metadata_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str}
    )
    usage = pd.read_csv(
        CNMF_FULL / "consensus_usage_rank7_v1.tsv.gz", sep="\t", dtype={"exact_cell_id": str}
    )
    clones = pd.read_csv(CLONES, sep="\t", dtype={"exact_cell_id": str})
    expected_order = cells.to_numpy(str)
    if not np.array_equal(base_meta["exact_cell_id"].to_numpy(str), expected_order):
        raise RuntimeError("cNMF metadata order differs from selected H5AD cells")
    if not np.array_equal(usage["exact_cell_id"].to_numpy(str), expected_order):
        raise RuntimeError("rank7 usage order differs from selected H5AD cells")
    factors = [f"F{i}" for i in range(1, 8)]
    if list(usage.columns) != ["exact_cell_id", *factors]:
        raise RuntimeError("Unexpected rank7 usage schema")
    usage_matrix = usage[factors].to_numpy(float)
    if np.any(~np.isfinite(usage_matrix)) or np.any(usage_matrix < 0):
        raise RuntimeError("Invalid rank7 usage values")
    dominant = np.asarray(factors, dtype=object)[np.argmax(usage_matrix, axis=1)]

    clone_cols = [
        "exact_cell_id",
        "formal_inferred_cnv_clone_v1",
        "clone_number",
        "clone_cells",
        "clone_size_class_v1",
        "eligible_for_primary_within_patient_clone_analysis_v1",
    ]
    meta = base_meta.merge(usage, on="exact_cell_id", how="left", validate="one_to_one")
    meta = meta.merge(clones[clone_cols], on="exact_cell_id", how="left", validate="one_to_one")
    meta["dominant_factor_rank7_v1"] = dominant
    meta["selected_TYMP_factor_v1"] = "F1"
    meta["selected_TYMP_factor_score_v1"] = usage_matrix[:, 0]
    if not np.array_equal(meta["exact_cell_id"].to_numpy(str), expected_order):
        raise RuntimeError("Merged metadata order differs from matrix")
    if meta["curated_patient_id"].nunique() != 26 or meta["study"].nunique() != 5:
        raise RuntimeError("Unexpected patient/study counts")

    gate("scenic_input_write_raw_counts_v1", gates)
    matrix_path = OUT / "raw_counts_genes_by_cells_v1.mtx.gz"
    with gzip.open(matrix_path, "wb") as handle:
        mmwrite(handle, counts.T.tocsc(), field="integer")
    pd.DataFrame({"gene": genes}).to_csv(OUT / "genes_v1.tsv", sep="\t", index=False)
    pd.DataFrame({"exact_cell_id": cells}).to_csv(OUT / "cells_v1.tsv", sep="\t", index=False)
    meta.to_csv(OUT / "cell_metadata_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.read_csv(CNMF_INPUT / "cnmf_hvg_genes_v1.tsv", sep="\t")[["gene"]].to_csv(
        OUT / "cnmf_rank7_model_features_v1.tsv", sep="\t", index=False, header=False
    )

    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    totals = np.asarray(counts.sum(axis=0)).ravel()
    gene_qc = pd.DataFrame({"gene": genes, "detected_cells": detected, "total_counts": totals})
    gene_qc.to_csv(OUT / "raw_gene_qc_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    factor_counts = (
        meta.groupby(["study", "sample_id", "dominant_factor_rank7_v1"], observed=False)
        .size().rename("n_cells").reset_index()
    )
    factor_counts.to_csv(OUT / "study_sample_factor_counts_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gates).to_csv(OUT / "memory_gate_log_v1.tsv", sep="\t", index=False)

    run_meta = {
        "schema_version": "part1_scenic_metaviper_input_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "input_h5ad": str(H5AD),
        "selection_rule": "eligible_for_subclone_selection_v1 == True",
        "matrix_orientation": "genes_by_cells",
        "n_genes": int(counts.shape[1]),
        "n_cells": int(counts.shape[0]),
        "nnz": int(counts.nnz),
        "n_patients": int(meta["curated_patient_id"].nunique()),
        "n_samples": int(meta["sample_id"].nunique()),
        "n_studies": int(meta["study"].nunique()),
        "selected_cnmf_rank": 7,
        "selected_tymp_factor": "F1",
        "seed": SEED,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(run_meta, indent=2) + "\n")
    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(json.dumps(run_meta, indent=2), flush=True)
    print("PART1_SCENIC_METAVIPER_INPUT_V1_PASS", flush=True)


if __name__ == "__main__":
    main()
