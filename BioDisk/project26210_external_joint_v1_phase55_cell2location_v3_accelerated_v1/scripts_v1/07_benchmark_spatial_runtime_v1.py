#!/usr/bin/env python3
"""Short, non-output benchmark for estimating Phase 55 spatial training time."""

from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyro
import scanpy as sc
import scvi
import torch
from cell2location.models import Cell2location
from scipy import sparse


def benchmark(path: Path, epochs: int, genes: int, cell_types: int) -> dict:
    started = time.perf_counter()
    spatial = sc.read_h5ad(path)
    if "counts" not in spatial.layers:
        raise RuntimeError(f"Counts layer missing from {path}")
    spatial.X = sparse.csr_matrix(spatial.layers["counts"])
    spatial.var_names = pd.Index(spatial.var_names.astype(str).str.upper())
    spatial.var_names_make_unique()
    spatial = spatial[:, : min(genes, spatial.n_vars)].copy()
    spatial.layers["counts"] = spatial.X.copy()
    sample = path.parent.name
    spatial.obs["section_batch"] = sample
    rng = np.random.default_rng(26210)
    signatures = pd.DataFrame(
        rng.uniform(0.5, 4.0, size=(spatial.n_vars, cell_types)),
        index=spatial.var_names,
        columns=[f"type_{index + 1:02d}" for index in range(cell_types)],
    )
    load_and_setup_seconds = time.perf_counter() - started
    Cell2location.setup_anndata(spatial, layer="counts", batch_key="section_batch")
    model = Cell2location(
        spatial,
        cell_state_df=signatures,
        N_cells_per_location=10,
        detection_alpha=20,
    )
    train_started = time.perf_counter()
    model.train(
        max_epochs=epochs,
        batch_size=None,
        train_size=1,
        lr=0.002,
        accelerator="cpu",
        enable_progress_bar=False,
    )
    train_seconds = time.perf_counter() - train_started
    result = {
        "sample_id": sample,
        "spots": int(spatial.n_obs),
        "genes": int(spatial.n_vars),
        "cell_types": cell_types,
        "benchmark_epochs": epochs,
        "load_and_setup_seconds": round(load_and_setup_seconds, 3),
        "train_seconds": round(train_seconds, 3),
        "seconds_per_epoch": train_seconds / epochs,
        "projected_20000_training_hours": train_seconds / epochs * 20000 / 3600,
    }
    del model, spatial, signatures
    pyro.clear_param_store()
    gc.collect()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sections", nargs="+", type=Path)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--genes", type=int, default=13805)
    parser.add_argument("--cell-types", type=int, default=12)
    args = parser.parse_args()
    scvi.settings.seed = 26210
    torch.set_num_threads(4)
    results = [benchmark(path, args.epochs, args.genes, args.cell_types) for path in args.sections]
    print(json.dumps({"status": "PASS", "results": results}, sort_keys=True))


if __name__ == "__main__":
    main()
