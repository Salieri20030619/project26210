#!/usr/bin/env python3
"""Build a fresh sparse cNMF input from the strict CNV-selected H5AD."""

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
INPUT = ROOT / "03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
PRIMARY = ROOT / "03_malignancy_cnv/malignancy_cnv_summary_v1/tumor_epithelial_malignancy_calls_v1.tsv.gz"
AUDIT_REGISTRY = ROOT / "04_subclones_tymp/input_audit_v1/subclone_cnmf_cell_registry_v1.tsv.gz"
OUT = ROOT / "05_cnmf/input_v1"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
EXPECTED_SHA256 = "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"
N_HVG = 2500
SEED = 26210


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gate(stage: str, records: list[dict[str, str]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    output = " ".join((result.stdout + " " + result.stderr).split())
    records.append({"stage": stage, "exit_status": str(result.returncode), "gate_output": output})
    print(output, flush=True)
    if result.returncode != 0 or "PASS" not in output:
        raise RuntimeError(f"Memory gate blocked {stage}: {output}")


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise RuntimeError(f"Refusing to overwrite non-empty cNMF input directory: {OUT}")
    OUT.mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    gate("open_h5ad", gates)
    adata = ad.read_h5ad(INPUT, backed="r")
    selected = adata.obs["eligible_for_subclone_selection_v1"].astype(bool).to_numpy()
    idx = np.flatnonzero(selected)
    if len(idx) != 11950:
        raise RuntimeError(f"Expected 11950 strict selected cells, found {len(idx)}")
    cells = adata.obs_names[idx].astype(str)
    genes = adata.var_names.astype(str)
    gate("read_selected_counts", gates)
    counts = sp.csr_matrix(adata.layers["counts"][idx, :], dtype=np.float64)
    adata.file.close()
    if counts.shape != (11950, 37795) or counts.nnz != 35713877:
        raise RuntimeError(f"Unexpected selected count matrix: shape={counts.shape}, nnz={counts.nnz}")
    if (
        np.any(~np.isfinite(counts.data))
        or np.any(counts.data < 0)
        or np.any(np.abs(counts.data - np.rint(counts.data)) > 1e-8)
    ):
        raise RuntimeError("Selected counts are not finite non-negative integers")
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(library_size <= 0):
        raise RuntimeError("Selected registry contains zero-library cells")

    gate("normalize_and_select_hvg", gates)
    log_cp10k = counts.multiply((1e4 / library_size)[:, None]).tocsr()
    log_cp10k.data = np.log1p(log_cp10k.data)
    detected = np.asarray((log_cp10k > 0).sum(axis=0)).ravel().astype(np.int64)
    mean = np.asarray(log_cp10k.mean(axis=0)).ravel()
    squared = log_cp10k.copy()
    squared.data **= 2
    variance = np.asarray(squared.mean(axis=0)).ravel() - mean**2
    variance[variance < 0] = 0
    dispersion = variance / (mean + 1e-6)
    dispersion[~np.isfinite(dispersion)] = 0
    upper = pd.Index(genes).str.upper()
    mito_unprefixed = {
        "COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8",
        "ND1", "ND2", "ND3", "ND4", "ND5", "ND6", "ND4L",
    }
    excluded = (
        upper.str.startswith("MT-")
        | upper.str.startswith("MTRNR")
        | upper.str.startswith("RPL")
        | upper.str.startswith("RPS")
        | upper.str.startswith("HBA")
        | upper.str.startswith("HBB")
        | upper.str.startswith("HBM")
        | upper.str.startswith("HBD")
        | upper.str.startswith("HBG")
        | upper.isin(mito_unprefixed)
    )
    detection_floor = max(30, int(np.floor(0.01 * counts.shape[0])))
    eligible_gene = (detected >= detection_floor) & ~np.asarray(excluded)
    eligible_idx = np.flatnonzero(eligible_gene)
    order = eligible_idx[np.argsort(-dispersion[eligible_idx], kind="stable")]
    hvg_idx = order[:N_HVG]
    if len(hvg_idx) != N_HVG or len(np.unique(hvg_idx)) != N_HVG:
        raise RuntimeError("Failed to select 2500 unique HVGs")
    hvg_genes = genes[hvg_idx]

    scaled = log_cp10k[:, hvg_idx].tocsc().T.tocsr()
    gene_rms = np.sqrt(np.asarray(scaled.power(2).mean(axis=1)).ravel())
    gene_rms[~np.isfinite(gene_rms) | (gene_rms <= 0)] = 1
    scaled = sp.diags(1.0 / gene_rms).dot(scaled).tocsr()
    scaled.data = np.minimum(scaled.data, 10.0)
    scaled.eliminate_zeros()
    if scaled.shape != (2500, 11950) or np.any(~np.isfinite(scaled.data)) or np.any(scaled.data < 0):
        raise RuntimeError(f"Invalid scaled cNMF matrix: {scaled.shape}")

    gene_qc = pd.DataFrame(
        {
            "gene": genes,
            "detected_cells": detected,
            "mean_log_cp10k": mean,
            "variance_log_cp10k": variance,
            "dispersion": dispersion,
            "excluded_technical": np.asarray(excluded),
            "eligible": eligible_gene,
        }
    )
    gene_qc["dispersion_rank_among_eligible"] = pd.NA
    gene_qc.loc[order, "dispersion_rank_among_eligible"] = np.arange(1, len(order) + 1)
    gene_qc["selected_hvg"] = False
    gene_qc.loc[hvg_idx, "selected_hvg"] = True
    gene_qc["gene_rms_if_hvg"] = np.nan
    gene_qc.loc[hvg_idx, "gene_rms_if_hvg"] = gene_rms
    gene_qc.to_csv(OUT / "cnmf_gene_qc_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame({"gene": hvg_genes, "gene_index_1based": np.arange(1, N_HVG + 1)}).to_csv(
        OUT / "cnmf_hvg_genes_v1.tsv", sep="\t", index=False
    )
    pd.DataFrame({"exact_cell_id": cells, "cell_index_1based": np.arange(1, len(cells) + 1)}).to_csv(
        OUT / "cnmf_cells_v1.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    gate("write_matrix_market", gates)
    with gzip.open(OUT / "cnmf_hvg_scaled_genes_by_cells_v1.mtx.gz", "wb") as handle:
        mmwrite(handle, scaled, field="real", precision=12)

    registry = pd.read_csv(AUDIT_REGISTRY, sep="\t", dtype={"exact_cell_id": str})
    primary = pd.read_csv(PRIMARY, sep="\t", dtype={"exact_cell_id": str})
    meta = registry.merge(
        primary[
            [
                "exact_cell_id", "infercnv_unit", "infercnv_group", "infercnv_cnv_gene_burden",
                "infercnv_high_cnv_threshold",
            ]
        ],
        on="exact_cell_id",
        how="left",
        validate="one_to_one",
    )
    if not np.array_equal(meta["exact_cell_id"].to_numpy(str), cells.to_numpy(str)):
        raise RuntimeError("cNMF metadata order differs from matrix cell order")
    tymp_idx = int(np.flatnonzero(genes == "TYMP")[0])
    tymp_log = np.asarray(log_cp10k[:, tymp_idx].toarray()).ravel()
    meta["TYMP_log1p_CP10K_v1"] = tymp_log
    meta["infercnv_burden_over_threshold_v1"] = (
        pd.to_numeric(meta["infercnv_cnv_gene_burden"], errors="coerce")
        / pd.to_numeric(meta["infercnv_high_cnv_threshold"], errors="coerce")
    )
    if meta[["infercnv_cnv_gene_burden", "infercnv_high_cnv_threshold"]].isna().any().any():
        raise RuntimeError("Strict cNMF cells have missing inferCNV burden metadata")
    meta.to_csv(OUT / "cnmf_cell_metadata_v1.tsv.gz", sep="\t", index=False, compression="gzip")

    qc = pd.DataFrame(
        {
            "metric": [
                "n_cells", "n_genes_input", "n_count_nnz", "library_size_min",
                "library_size_median", "library_size_max", "detection_floor",
                "n_eligible_genes", "n_hvg", "scaled_rows", "scaled_columns",
                "scaled_nnz", "scaled_min_nonzero", "scaled_max_nonzero",
                "scaled_nonfinite", "scaled_negative", "TYMP_selected_hvg",
                "TYMP_positive_cells", "metadata_missing_cnv_burden",
            ],
            "value": [
                counts.shape[0], counts.shape[1], counts.nnz, library_size.min(),
                np.median(library_size), library_size.max(), detection_floor,
                int(eligible_gene.sum()), len(hvg_idx), scaled.shape[0], scaled.shape[1],
                scaled.nnz, scaled.data.min(), scaled.data.max(),
                int(np.count_nonzero(~np.isfinite(scaled.data))), int(np.count_nonzero(scaled.data < 0)),
                bool("TYMP" in set(hvg_genes)), int(np.count_nonzero(tymp_log > 0)),
                int(meta["infercnv_cnv_gene_burden"].isna().sum()),
            ],
        }
    )
    qc.to_csv(OUT / "cnmf_input_qc_v1.tsv", sep="\t", index=False)
    pd.DataFrame(gates).to_csv(OUT / "memory_gate_log_v1.tsv", sep="\t", index=False)
    run_metadata = {
        "schema_version": "part1_cnmf_input_v1",
        "completed_at": datetime.now().astimezone().isoformat(),
        "input_h5ad": str(INPUT),
        "input_h5ad_sha256": EXPECTED_SHA256,
        "selection_rule": "eligible_for_subclone_selection_v1 == True",
        "normalization": "log1p(count / cell_library_size * 10000)",
        "technical_gene_exclusion": "mitochondrial,ribosomal,hemoglobin",
        "minimum_detected_cells": detection_floor,
        "hvg_method": "variance/mean dispersion on log1p CP10K",
        "n_hvg": N_HVG,
        "scaling": "per-gene RMS without centering; cap at 10",
        "matrix_orientation": "genes_by_cells",
        "seed": SEED,
        "n_cells": 11950,
        "n_patients": int(meta["curated_patient_id"].nunique()),
        "dependency_changes": False,
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(run_metadata, indent=2) + "\n", encoding="utf-8")
    manifest = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name not in {"file_manifest_v1.tsv", "RUN_COMPLETE"}:
            manifest.append({"file": path.name, "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(manifest).to_csv(OUT / "file_manifest_v1.tsv", sep="\t", index=False)
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n", encoding="utf-8")
    print(json.dumps(run_metadata, indent=2))
    print("PART1_CNMF_INPUT_V1_PASS")


if __name__ == "__main__":
    main()
