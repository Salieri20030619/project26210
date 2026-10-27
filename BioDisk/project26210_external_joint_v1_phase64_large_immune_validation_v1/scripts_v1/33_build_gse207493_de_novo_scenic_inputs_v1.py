#!/usr/bin/env python3
"""Create lineage-specific, sample-balanced de novo pySCENIC inputs."""

from __future__ import annotations

import gzip
import json
from datetime import datetime
from pathlib import Path

import loompy
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csc_matrix


PROJECT = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1")
SOURCE = PROJECT / "inputs_v1/GSE207493/tf_regulon_v1"
OUT = PROJECT / "inputs_v1/GSE207493/tf_regulon_de_novo_v1/scenic_inputs_v1"
SEED = 26210
LINEAGES = {
    "tumor": ({"Tumor"}, 250),
    "myeloid_dc": ({"Myeloid", "Dendritic"}, 250),
    "treg": ({"Treg"}, 150),
}


def make_loom(path: Path, matrix: csc_matrix, genes: np.ndarray, meta: pd.DataFrame) -> None:
    loompy.create(
        str(path), matrix.astype(np.float32),
        row_attrs={"Gene": genes.astype(object)},
        col_attrs={
            "CellID": meta.cell_barcode.astype(str).to_numpy(object),
            "sample": meta["sample"].astype(str).to_numpy(object),
            "cell_type": meta.cell_type.astype(str).to_numpy(object),
            "F1_score": meta.f1_score.to_numpy(float),
            "TYMP_log1p": meta.TYMP_log1p.to_numpy(float),
            "S_score": meta.S_score.to_numpy(float),
            "G2M_score": meta.G2M_score.to_numpy(float),
        },
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "RUN_COMPLETE").exists():
        print("DE_NOVO_SCENIC_INPUTS_V1_REUSED")
        return
    with gzip.open(SOURCE / "raw_counts_genes_by_cells_v1.mtx.gz", "rb") as handle:
        counts = csc_matrix(mmread(handle), dtype=np.float32)
    genes = pd.read_csv(SOURCE / "genes_v1.tsv", sep="\t")["gene"].astype(str).to_numpy(dtype=str)
    cells = pd.read_csv(SOURCE / "cells_v1.tsv", sep="\t")["cell_barcode"].astype(str).to_numpy(dtype=str)
    meta = pd.read_csv(SOURCE / "cell_metadata_v1.tsv.gz", sep="\t", dtype={"cell_barcode": str})
    if counts.shape != (len(genes), len(cells)) or not np.array_equal(meta.cell_barcode, cells):
        raise RuntimeError("Frozen target-lineage source mismatch")

    rng = np.random.default_rng(SEED)
    summary_rows = []
    for lineage, (types, cap) in LINEAGES.items():
        lineage_dir = OUT / lineage
        lineage_dir.mkdir(parents=True, exist_ok=True)
        cell_mask = meta.cell_type.isin(types).to_numpy()
        one_meta = meta.loc[cell_mask].reset_index(drop=True)
        one_counts = counts[:, cell_mask].tocsc()
        detected = np.asarray((one_counts > 0).sum(axis=1)).ravel()
        totals = np.asarray(one_counts.sum(axis=1)).ravel()
        min_cells = max(20, int(np.ceil(0.01 * one_counts.shape[1])))
        upper = np.char.upper(genes)
        gene_mask = (
            (detected >= min_cells) & (totals >= 3 * min_cells)
            & (~np.char.startswith(upper, "MT-")) & (~np.char.startswith(upper, "MTRNR"))
        )
        filtered = one_counts[gene_mask, :].tocsc()
        filtered_genes = genes[gene_mask]
        selected = []
        sample_rows = []
        for sample, group in one_meta.groupby("sample", sort=True):
            indices = group.index.to_numpy()
            if len(indices) > cap:
                chosen = np.sort(rng.choice(indices, cap, replace=False))
            else:
                chosen = indices
            selected.extend(chosen.tolist())
            sample_rows.append({
                "lineage": lineage, "sample": sample, "n_all_cells": len(indices),
                "n_grn_cells": len(chosen), "max_cells_per_sample": cap,
            })
        selected = np.asarray(sorted(selected), dtype=int)
        make_loom(lineage_dir / "all_cells_filtered_counts_v1.loom", filtered, filtered_genes, one_meta)
        make_loom(
            lineage_dir / "balanced_grn_counts_v1.loom", filtered[:, selected].tocsc(), filtered_genes,
            one_meta.iloc[selected].reset_index(drop=True),
        )
        pd.DataFrame({
            "gene": genes, "detected_cells": detected, "total_counts": totals,
            "keep_for_scenic": gene_mask,
        }).to_csv(lineage_dir / "gene_qc_v1.tsv.gz", sep="\t", index=False, compression="gzip")
        pd.DataFrame(sample_rows).to_csv(lineage_dir / "sample_balancing_v1.tsv", sep="\t", index=False)
        one_meta.to_csv(lineage_dir / "cell_metadata_v1.tsv.gz", sep="\t", index=False, compression="gzip")
        summary_rows.append({
            "lineage": lineage, "cell_types": ",".join(sorted(types)), "n_all_cells": len(one_meta),
            "n_grn_cells": len(selected), "n_filtered_genes": int(gene_mask.sum()),
            "min_detected_cells_per_gene": min_cells, "max_cells_per_sample": cap,
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "de_novo_scenic_input_summary_v1.tsv", sep="\t", index=False)
    (OUT / "run_metadata_v1.json").write_text(json.dumps({
        "completed_at": datetime.now().astimezone().isoformat(), "seed": SEED,
        "network_learning": "separate by lineage; deterministic sample-balanced raw counts",
        "old_part1_regulons_used": False,
    }, indent=2) + "\n")
    (OUT / "RUN_COMPLETE").write_text(datetime.now().astimezone().isoformat() + "\n")
    print(summary.to_string(index=False))
    print("DE_NOVO_SCENIC_INPUTS_V1_PASS")


if __name__ == "__main__":
    main()
