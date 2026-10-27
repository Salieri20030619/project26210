#!/usr/bin/env python3
"""Build the frozen 11,950-cell raw-count scDRS input without loading source X."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    part = path.with_name(path.name + ".part")
    if part.exists():
        raise FileExistsError(part)
    part.write_text(text, encoding="utf-8", newline="\n")
    os.replace(part, path)


def bool_series(series: pd.Series) -> pd.Series:
    true_values = {"true", "t", "1", "yes"}
    false_values = {"false", "f", "0", "no"}
    lowered = series.astype(str).str.strip().str.lower()
    unknown = sorted(set(lowered) - true_values - false_values)
    if unknown:
        raise ValueError(f"unrecognized Boolean values: {unknown}")
    return lowered.isin(true_values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--cell-registry", required=True, type=Path)
    parser.add_argument("--cnmf-usage", required=True, type=Path)
    parser.add_argument("--gene-set", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cells", type=int, default=11950)
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()

    for source in (args.h5ad, args.cell_registry, args.cnmf_usage, args.gene_set):
        if not source.is_file():
            raise FileNotFoundError(source)
    if args.output.exists():
        raise FileExistsError(f"refusing existing output: {args.output}")
    part_dir = args.output.with_name(args.output.name + ".part")
    if part_dir.exists():
        raise FileExistsError(f"refusing stale partial output: {part_dir}")
    part_dir.mkdir(parents=True)

    registry = pd.read_csv(args.cell_registry, sep="\t", low_memory=False)
    required_registry = {
        "exact_cell_id", "study", "curated_patient_id", "sample_id",
        "eligible_for_subclone_selection_v1", "raw_library_size_v1",
        "raw_detected_genes_v1", "TYMP_raw_count_v1",
    }
    missing = required_registry - set(registry.columns)
    if missing:
        raise ValueError(f"registry columns missing: {sorted(missing)}")
    if registry["exact_cell_id"].duplicated().any():
        raise ValueError("duplicate exact_cell_id in registry")
    registry = registry.loc[bool_series(registry["eligible_for_subclone_selection_v1"])].copy()
    if len(registry) != args.expected_cells:
        raise ValueError(f"eligible registry cells {len(registry)} != {args.expected_cells}")
    ids = registry["exact_cell_id"].astype(str).tolist()

    usage = pd.read_csv(args.cnmf_usage, sep="\t")
    factor_cols = [f"F{i}" for i in range(1, 8)]
    if not {"exact_cell_id", *factor_cols}.issubset(usage.columns):
        raise ValueError("rank-7 cNMF usage columns are incomplete")
    if usage["exact_cell_id"].duplicated().any():
        raise ValueError("duplicate exact_cell_id in cNMF usage")
    usage = usage.set_index("exact_cell_id").reindex(ids)
    if usage[factor_cols].isna().any().any():
        raise ValueError("eligible registry cells missing from cNMF usage")
    vals = usage[factor_cols].to_numpy(float)
    if not np.isfinite(vals).all() or (vals < 0).any():
        raise ValueError("invalid cNMF usage values")

    src = ad.read_h5ad(args.h5ad, backed="r")
    if not src.obs_names.is_unique or not src.var_names.is_unique:
        raise ValueError("source H5AD identifiers are not unique")
    if "counts" not in src.layers:
        raise ValueError("source H5AD lacks counts layer")
    source_index = pd.Index(src.obs_names.astype(str))
    positions = source_index.get_indexer(ids)
    if (positions < 0).any():
        raise ValueError(f"{int((positions < 0).sum())} registry cells missing from H5AD")

    blocks = []
    for start in range(0, len(positions), args.chunk_size):
        block = src.layers["counts"][positions[start:start + args.chunk_size], :]
        block = sp.csr_matrix(block)
        if not np.isfinite(block.data).all():
            raise ValueError("non-finite raw counts")
        if (block.data < 0).any():
            raise ValueError("negative raw counts")
        if not np.equal(block.data, np.floor(block.data)).all():
            raise ValueError("non-integer raw counts")
        blocks.append(block)
    counts = sp.vstack(blocks, format="csr")
    del blocks

    library_size = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    detected_genes = np.diff(counts.indptr).astype(int)
    if not np.allclose(
        library_size,
        pd.to_numeric(registry["raw_library_size_v1"]).to_numpy(float),
        rtol=0,
        atol=1e-8,
    ):
        raise ValueError("raw library sizes disagree with frozen registry")
    if not np.array_equal(
        detected_genes,
        pd.to_numeric(registry["raw_detected_genes_v1"]).to_numpy(int),
    ):
        raise ValueError("detected-gene counts disagree with frozen registry")

    obs = registry.set_index("exact_cell_id", drop=False)
    for col in factor_cols:
        obs[col] = usage[col].to_numpy(float)
    obs["dominant_factor_rank7_v1"] = np.array(factor_cols)[np.argmax(vals, axis=1)]
    obs["scdrs_const_v1"] = 1.0
    obs["scdrs_n_genes_v1"] = detected_genes
    obs["scdrs_log_n_counts_v1"] = np.log1p(library_size)

    obsm = {}
    if "X_umap" in src.obsm:
        obsm["X_umap"] = np.asarray(src.obsm["X_umap"][positions, :])
    out = ad.AnnData(
        X=counts,
        obs=obs,
        var=src.var.copy(),
        obsm=obsm,
    )
    out.obs_names = ids
    out.uns["phase25_input_contract_v1"] = {
        "cell_order": "eligible rebuilt registry order",
        "matrix": "unmodified source counts layer",
        "expected_cells": int(args.expected_cells),
        "cnmf_rank": 7,
    }
    h5ad_path = part_dir / "malignant_eligible_raw_counts_phase25_v1.h5ad"
    out.write_h5ad(h5ad_path, compression="gzip")
    src.file.close()

    cov = pd.DataFrame(
        {
            "cell_id": ids,
            "const": 1.0,
            "n_genes": detected_genes,
            "log_n_counts": np.log1p(library_size),
        }
    ).set_index("cell_id")
    cov.to_csv(part_dir / "malignant_eligible_phase25_v1.cov", sep="\t", index=True)
    shutil.copy2(args.gene_set, part_dir / "RCC_GWAS_Catalog_loci_phase25_v1.gs")

    gs_lines = args.gene_set.read_text(encoding="utf-8-sig").splitlines()
    if len(gs_lines) < 2:
        raise ValueError("gene-set file has no data row")
    tokens = gs_lines[1].split("\t", 1)[1].split(",")
    genes = [token.split(":", 1)[0] for token in tokens]
    expressed = set(map(str, out.var_names))
    coverage = pd.DataFrame({"gene": genes, "present_in_matrix": [g in expressed for g in genes]})
    coverage.to_csv(part_dir / "gene_set_coverage_v1.tsv", sep="\t", index=False)
    if int(coverage["present_in_matrix"].sum()) < 10:
        raise ValueError("fewer than 10 disease genes are represented")

    versions = {
        "python": platform.python_version(),
        "anndata": version("anndata"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "scipy": version("scipy"),
    }
    atomic_text(part_dir / "software_versions_v1.json", json.dumps(versions, indent=2) + "\n")
    rows = []
    for file in sorted(p for p in part_dir.iterdir() if p.is_file()):
        rows.append({"file": file.name, "bytes": file.stat().st_size, "sha256": sha256(file)})
    pd.DataFrame(rows).to_csv(part_dir / "INPUT_FILE_MANIFEST_v1.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_INPUT_COMPLETE",
        "n_cells": int(out.n_obs),
        "n_genes": int(out.n_vars),
        "n_patients": int(obs["curated_patient_id"].nunique()),
        "n_studies": int(obs["study"].nunique()),
        "gene_set_size": len(genes),
        "gene_set_present": int(coverage["present_in_matrix"].sum()),
        "counts_nnz": int(counts.nnz),
    }
    atomic_text(part_dir / "status.json", json.dumps(status, indent=2) + "\n")
    os.replace(part_dir, args.output)


if __name__ == "__main__":
    main()

