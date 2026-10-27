#!/usr/bin/env python3
"""Run official scDRS compute-score operations without the optional Fire CLI parser."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


EXPECTED_COMMIT = "638b377bbe906a9116d3894b13386a2dc67301d8"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git(source: Path, *args: str) -> str:
    command = ["git", "-c", f"safe.directory={source}", "-C", str(source), *args]
    return subprocess.check_output(command, text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scdrs-source", required=True, type=Path)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--cov", required=True, type=Path)
    parser.add_argument("--gene-set", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-cells", type=int, default=11950)
    parser.add_argument("--n-ctrl", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-genes", type=int, default=1)
    parser.add_argument("--min-cells", type=int, default=50)
    args = parser.parse_args()

    for source in (args.scdrs_source, args.h5ad, args.cov, args.gene_set):
        if not source.exists():
            raise FileNotFoundError(source)
    if args.output.exists():
        raise FileExistsError(f"refusing existing output: {args.output}")
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(f"refusing stale partial output: {part}")
    part.mkdir(parents=True)
    score_dir = part / "scores"
    score_dir.mkdir()

    commit = git(args.scdrs_source, "rev-parse", "HEAD")
    describe = git(args.scdrs_source, "describe", "--tags", "--always")
    worktree = git(args.scdrs_source, "status", "--short")
    if commit != EXPECTED_COMMIT:
        raise ValueError(f"unexpected scDRS commit: {commit}")
    if worktree:
        raise ValueError("scDRS source worktree is not clean")

    sys.path.insert(0, str(args.scdrs_source))
    import scdrs  # noqa: PLC0415

    start = time.time()
    adata = scdrs.util.load_h5ad(
        str(args.h5ad),
        flag_filter_data=True,
        flag_raw_count=True,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
    )
    if adata.n_obs != args.expected_cells:
        raise ValueError(f"scDRS filtering retained {adata.n_obs}, expected {args.expected_cells}")
    cov = pd.read_csv(args.cov, sep="\t", index_col=0)
    cov.index = cov.index.astype(str)
    if set(cov.index) != set(adata.obs_names.astype(str)):
        raise ValueError("covariate/H5AD cell identity mismatch")
    cov = cov.reindex(adata.obs_names.astype(str))
    gene_sets = scdrs.util.load_gs(
        str(args.gene_set),
        src_species="human",
        dst_species="human",
        to_intersect=adata.var_names,
    )
    if len(gene_sets) != 1:
        raise ValueError(f"expected one gene set, found {len(gene_sets)}")
    trait = next(iter(gene_sets))
    genes, weights = gene_sets[trait]
    if len(genes) < 10:
        raise ValueError(f"gene set has only {len(genes)} genes after filtering")

    scdrs.preprocess(adata, cov=cov, adj_prop=None, n_mean_bin=20, n_var_bin=20, copy=False)
    np.random.seed(args.seed)
    result = scdrs.score_cell(
        adata,
        genes,
        gene_weight=weights,
        ctrl_match_key="mean_var",
        n_ctrl=args.n_ctrl,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
    if len(result) != args.expected_cells or set(result.index.astype(str)) != set(adata.obs_names.astype(str)):
        raise ValueError("scDRS result cell identity/count gate failed")
    if result.iloc[:, :6].isna().any().any():
        raise ValueError("primary scDRS output contains missing values")

    score_path = score_dir / f"{trait}.score.gz"
    full_path = score_dir / f"{trait}.full_score.gz"
    result.iloc[:, :6].to_csv(score_path, sep="\t", index=True, compression="gzip")
    result.to_csv(full_path, sep="\t", index=True, compression="gzip")
    fdr = multipletests(result["pval"].to_numpy(float), method="fdr_bh")[1]

    source_hash_rows = []
    for path in sorted(args.scdrs_source.rglob("*")):
        if path.is_file() and ".git" not in path.parts:
            source_hash_rows.append(
                {"relative_file": str(path.relative_to(args.scdrs_source)), "sha256": sha256(path)}
            )
    pd.DataFrame(source_hash_rows).to_csv(part / "scdrs_source_sha256_v1.tsv", sep="\t", index=False)
    run = {
        "status": "PASS_SCDRS_COMPUTE",
        "scdrs_version": str(scdrs.__version__),
        "source_commit": commit,
        "source_describe": describe,
        "source_worktree_clean": True,
        "python": platform.python_version(),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "n_cells": int(len(result)),
        "n_genes_after_filter": int(adata.n_vars),
        "trait": trait,
        "trait_genes_after_filter": int(len(genes)),
        "n_ctrl": int(args.n_ctrl),
        "seed": int(args.seed),
        "flag_filter_data": True,
        "flag_raw_count": True,
        "min_genes": int(args.min_genes),
        "min_cells": int(args.min_cells),
        "ctrl_match_opt": "mean_var",
        "weight_opt": "vs",
        "n_cells_fdr_below_0_10_descriptive": int((fdr < 0.10).sum()),
        "n_cells_fdr_below_0_20_descriptive": int((fdr < 0.20).sum()),
        "elapsed_seconds": float(time.time() - start),
    }
    (part / "run_metadata_v1.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    manifest_rows = []
    for path in sorted(p for p in part.rglob("*") if p.is_file()):
        manifest_rows.append(
            {"relative_file": str(path.relative_to(part)), "bytes": path.stat().st_size, "sha256": sha256(path)}
        )
    pd.DataFrame(manifest_rows).to_csv(part / "COMPUTE_FILE_MANIFEST_v1.tsv", sep="\t", index=False)
    (part / "status.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output)


if __name__ == "__main__":
    main()
