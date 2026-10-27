#!/usr/bin/env python3
"""Build pySCENIC loom from the frozen target-lineage Matrix Market input."""

from __future__ import annotations

import gzip
from pathlib import Path

import loompy
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix


def main() -> None:
    project = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
    root = project / "inputs_v1/GSE207493/tf_regulon_v1"
    output = root / "target_lineages_filtered_counts_v1.loom"
    if output.exists() and output.stat().st_size > 0:
        print("GSE207493_SCENIC_LOOM_V1_REUSED")
        return
    if not (root / "RUN_COMPLETE").exists():
        raise RuntimeError("Frozen TF input is incomplete")
    with gzip.open(root / "raw_counts_genes_by_cells_v1.mtx.gz", "rb") as handle:
        counts = csc_matrix(mmread(handle), dtype=np.float32)
    genes = pd.read_csv(root / "genes_v1.tsv", sep="\t")["gene"].astype(str).to_numpy()
    cells = pd.read_csv(root / "cells_v1.tsv", sep="\t")["cell_barcode"].astype(str).to_numpy()
    meta = pd.read_csv(root / "cell_metadata_v1.tsv.gz", sep="\t", dtype={"cell_barcode": str})
    if counts.shape != (len(genes), len(cells)) or not np.array_equal(meta.cell_barcode, cells):
        raise RuntimeError("Frozen matrix order/dimensions changed")
    detected = np.asarray((counts > 0).sum(axis=1)).ravel()
    totals = np.asarray(counts.sum(axis=1)).ravel()
    min_cells = max(20, int(np.ceil(0.01 * len(cells))))
    upper = np.char.upper(genes)
    keep = (
        (detected >= min_cells) & (totals >= 3 * min_cells)
        & (~np.char.startswith(upper, "MT-")) & (~np.char.startswith(upper, "MTRNR"))
    )
    scenic_counts = counts[keep, :].tocsc()
    scenic_genes = genes[keep]
    loompy.create(
        str(output), scenic_counts,
        row_attrs={"Gene": scenic_genes.astype(object)},
        col_attrs={
            "CellID": cells.astype(object),
            "sample": meta["sample"].astype(str).to_numpy(object),
            "cell_type": meta["cell_type"].astype(str).to_numpy(object),
            "F1_score": meta["f1_score"].to_numpy(float),
            "TYMP_log1p": meta["TYMP_log1p"].to_numpy(float),
            "S_score": meta["S_score"].to_numpy(float),
            "G2M_score": meta["G2M_score"].to_numpy(float),
        },
    )
    pd.DataFrame({"gene": genes, "detected_cells": detected, "total_counts": totals,
                  "scenic_keep": keep}).to_csv(
        root / "gene_filter_qc_v1.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    print(f"GSE207493_SCENIC_LOOM_V1_PASS genes={len(scenic_genes)} cells={len(cells)}")


if __name__ == "__main__":
    main()
