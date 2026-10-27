#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import loompy
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix


def create_loom(path: Path, matrix: csc_matrix, genes: np.ndarray, meta: pd.DataFrame) -> None:
    loompy.create(
        str(path),
        matrix.astype(np.float32),
        row_attrs={"Gene": genes.astype(object)},
        col_attrs={
            "CellID": meta["cell_id"].astype(str).to_numpy(dtype=object),
            "patient": meta["patient"].astype(str).to_numpy(dtype=object),
            "study": meta["study_v1"].astype(str).to_numpy(dtype=object),
            "state": meta["state_v2"].astype(str).to_numpy(dtype=object),
            "DPT": meta["dpt_pseudotime"].to_numpy(float),
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    raw = csc_matrix(mmread(args.input / "raw_counts_genes_by_cells_v1.mtx.gz"))
    genes = pd.read_csv(args.input / "genes_v1.tsv", sep="\t")["gene"].astype(str).to_numpy()
    cells = pd.read_csv(args.input / "cells_v1.tsv", sep="\t")["cell_id"].astype(str).to_numpy()
    meta = pd.read_csv(args.input / "cell_metadata_v1.tsv.gz", sep="\t", dtype={"cell_id": str})
    if raw.shape != (len(genes), len(cells)) or not np.array_equal(meta["cell_id"].to_numpy(str), cells):
        raise RuntimeError("Matrix and metadata alignment failed")
    detected = np.asarray((raw > 0).sum(axis=1)).ravel()
    totals = np.asarray(raw.sum(axis=1)).ravel()
    min_cells = max(50, int(np.ceil(0.01 * raw.shape[1])))
    min_counts = 3 * min_cells
    upper = np.char.upper(genes.astype(str))
    mito = np.char.startswith(upper, "MT-") | np.char.startswith(upper, "MTRNR")
    keep = (detected >= min_cells) & (totals >= min_counts) & (~mito)
    filtered = raw[keep, :].tocsc()
    filtered_genes = genes[keep]
    discovery = meta["selected_for_grn_discovery_v1"].astype(bool).to_numpy()
    if filtered.shape[0] < 5000 or discovery.sum() < 5000:
        raise RuntimeError("Filtered gene or discovery-cell scope unexpectedly small")
    create_loom(args.output / "all_myeloid_filtered_counts_v1.loom", filtered, filtered_genes, meta)
    create_loom(args.output / "balanced_discovery_counts_v1.loom", filtered[:, discovery], filtered_genes, meta.loc[discovery].reset_index(drop=True))
    qc = pd.DataFrame({
        "metric": ["n_input_genes", "n_filtered_genes", "n_all_cells", "n_discovery_cells", "minimum_detected_cells", "minimum_total_counts"],
        "value": [len(genes), len(filtered_genes), len(cells), int(discovery.sum()), min_cells, min_counts],
    })
    qc.to_csv(args.output / "loom_qc_v1.tsv", sep="\t", index=False)
    (args.output / "status.json").write_text(json.dumps({
        "status": "PASS_Y03_PYSCENIC_LOOM_BUILD",
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_filtered_genes": int(len(filtered_genes)),
        "n_all_cells": int(len(cells)),
        "n_discovery_cells": int(discovery.sum()),
    }, indent=2) + "\n", encoding="utf-8")
    print(qc.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
