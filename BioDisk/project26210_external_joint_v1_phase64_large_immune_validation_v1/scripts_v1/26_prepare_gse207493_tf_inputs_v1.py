#!/usr/bin/env python3
"""Build frozen GSE207493 target-lineage inputs for external TF/regulon validation."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path

for variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "8")

import anndata as ad
import loompy
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite


SEED = 26210
EXPECTED_SAMPLES = [
    "RCC81", "RCC84", "RCC86", "RCC87", "RCC94", "RCC96", "RCC99",
    "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112",
    "RCC113", "RCC114", "RCC115", "RCC116", "RCC119", "RCC120",
]
MIN_CELLS = {
    "RCC81": 8, "RCC84": 7, "RCC86": 9, "RCC87": 9, "RCC94": 8,
    "RCC96": 8, "RCC99": 9, "RCC100": 10, "RCC101": 10,
    "RCC103": 10, "RCC104": 10, "RCC106": 10, "RCC112": 10,
    "RCC113": 9, "RCC114": 7, "RCC115": 9, "RCC116": 6,
    "RCC119": 10, "RCC120": 10,
}
MAX_FEATURES = {
    "RCC81": 3500, "RCC84": 4700, "RCC86": 5000, "RCC87": 2900,
    "RCC94": 3700, "RCC96": 3300, "RCC99": 2800, "RCC100": 3800,
    "RCC101": 2400, "RCC103": 3300, "RCC104": 5000,
    "RCC106": 4000, "RCC112": 4000, "RCC113": 4500,
    "RCC114": 3200, "RCC115": 2400, "RCC116": 4500,
    "RCC119": 4000, "RCC120": 3500,
}
TARGET_TYPES = {"Tumor", "Myeloid", "Dendritic", "Treg"}
EXPECTED_ALL_CELLS = 102723
EXPECTED_TARGET_CELLS = 39912


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    project = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
    raw_root = project / "inputs_v1/GSE207493/raw_extracted_v1/scRNA-seq"
    result_root = project / "results_v1/GSE207493"
    output = project / "inputs_v1/GSE207493/tf_regulon_v1"
    output.mkdir(parents=True, exist_ok=True)
    if (output / "RUN_COMPLETE").exists():
        print("GSE207493_TF_INPUTS_V1_REUSED")
        return

    cycle = pd.read_csv(
        result_root / "cell_cycle_v1/gse207493_cell_cycle_scores_v1.tsv.gz", sep="\t",
        dtype={"cell_barcode": str, "sample": str},
    ).set_index("cell_barcode", drop=False)
    if len(cycle) != EXPECTED_ALL_CELLS or cycle.index.has_duplicates:
        raise RuntimeError(f"Frozen metadata mismatch: {cycle.shape}")

    selected: list[ad.AnnData] = []
    qc_rows: list[dict[str, object]] = []
    total_qc = 0
    for sample in EXPECTED_SAMPLES:
        directory = raw_root / sample
        one = sc.read_10x_mtx(directory, var_names="gene_symbols", make_unique=True)
        raw_cells, raw_genes = one.n_obs, one.n_vars
        sc.pp.filter_genes(one, min_cells=MIN_CELLS[sample])
        one.var["mt"] = one.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(one, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            (one.obs["n_genes_by_counts"] > 500)
            & (one.obs["n_genes_by_counts"] < MAX_FEATURES[sample])
            & (one.obs["pct_counts_mt"] < 10)
        )
        one = one[keep].copy()
        one.obs_names = [f"{sample}_{barcode}" for barcode in one.obs_names]
        total_qc += one.n_obs
        missing = one.obs_names.difference(cycle.index)
        if len(missing):
            raise RuntimeError(f"{sample}: {len(missing)} QC cells absent from frozen metadata")
        sample_meta = cycle.loc[one.obs_names]
        target = sample_meta["cell_type"].isin(TARGET_TYPES).to_numpy()
        target_one = one[target].copy()
        target_meta = sample_meta.loc[target_one.obs_names]
        for column in [
            "sample", "leiden", "cell_type", "phase", "S_score", "G2M_score",
            "cycling_score", "f1_score", "TYMP_log1p",
        ]:
            target_one.obs[column] = target_meta[column].to_numpy()
        selected.append(target_one)
        qc_rows.append({
            "sample": sample, "raw_cells": raw_cells, "raw_genes": raw_genes,
            "qc_cells": one.n_obs, "target_cells": target_one.n_obs,
            "tumor_cells": int((target_meta.cell_type == "Tumor").sum()),
            "myeloid_dc_cells": int(target_meta.cell_type.isin(["Myeloid", "Dendritic"]).sum()),
            "treg_cells": int((target_meta.cell_type == "Treg").sum()),
        })

    if total_qc != EXPECTED_ALL_CELLS:
        raise RuntimeError(f"QC cell count changed: {total_qc}")
    merged = ad.concat(selected, join="outer", merge="same", fill_value=0)
    merged.var_names = merged.var_names.astype(str)
    if merged.n_obs != EXPECTED_TARGET_CELLS or merged.obs_names.has_duplicates:
        raise RuntimeError(f"Target cell count mismatch: {merged.n_obs}")
    if not np.all(np.isfinite(merged.X.data)) or np.any(merged.X.data < 0):
        raise RuntimeError("Invalid counts")

    counts = sparse.csc_matrix(merged.X.T, dtype=np.float32)
    genes = merged.var_names.to_numpy(str)
    cells = merged.obs_names.to_numpy(str)
    matrix_path = output / "raw_counts_genes_by_cells_v1.mtx.gz"
    with gzip.open(matrix_path, "wb") as handle:
        mmwrite(handle, counts, field="integer")
    pd.DataFrame({"gene": genes}).to_csv(output / "genes_v1.tsv", sep="\t", index=False)
    pd.DataFrame({"cell_barcode": cells}).to_csv(output / "cells_v1.tsv", sep="\t", index=False)
    metadata_columns = [
        "sample", "leiden", "cell_type", "phase", "S_score", "G2M_score",
        "cycling_score", "f1_score", "TYMP_log1p",
    ]
    metadata = merged.obs[metadata_columns].reset_index().rename(columns={"index": "cell_barcode"})
    metadata.to_csv(output / "cell_metadata_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    detected = np.asarray((counts > 0).sum(axis=1)).ravel()
    totals = np.asarray(counts.sum(axis=1)).ravel()
    min_cells = max(20, int(np.ceil(0.01 * merged.n_obs)))
    upper = np.char.upper(genes)
    keep_gene = (
        (detected >= min_cells) & (totals >= 3 * min_cells)
        & (~np.char.startswith(upper, "MT-")) & (~np.char.startswith(upper, "MTRNR"))
    )
    scenic_counts = counts[keep_gene, :].tocsc()
    scenic_genes = genes[keep_gene]
    loompy.create(
        str(output / "target_lineages_filtered_counts_v1.loom"),
        scenic_counts,
        row_attrs={"Gene": scenic_genes.astype(object)},
        col_attrs={
            "CellID": cells.astype(object),
            "sample": metadata["sample"].astype(str).to_numpy(object),
            "cell_type": metadata["cell_type"].astype(str).to_numpy(object),
            "F1_score": metadata["f1_score"].to_numpy(float),
            "TYMP_log1p": metadata["TYMP_log1p"].to_numpy(float),
            "S_score": metadata["S_score"].to_numpy(float),
            "G2M_score": metadata["G2M_score"].to_numpy(float),
        },
    )
    pd.DataFrame(qc_rows).to_csv(output / "sample_target_qc_v1.tsv", sep="\t", index=False)
    run_meta = {
        "schema_version": "gse207493_tf_regulon_input_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "seed": SEED, "n_all_qc_cells": total_qc, "n_target_cells": merged.n_obs,
        "n_raw_union_genes": merged.n_vars, "n_scenic_filtered_genes": int(keep_gene.sum()),
        "target_cell_types": sorted(TARGET_TYPES),
        "purpose": "fixed-network external validation; no network relearning",
    }
    (output / "run_metadata_v1.json").write_text(json.dumps(run_meta, indent=2) + "\n")
    manifest = []
    for path in sorted(output.iterdir()):
        if path.is_file() and path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(output / "file_manifest_v1.tsv", sep="\t", index=False)
    (output / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(json.dumps(run_meta, indent=2))
    print("GSE207493_TF_INPUTS_V1_PASS")


if __name__ == "__main__":
    main()
