#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


TARGETS = ["TYMP", "RRM2B"]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True, type=Path)
    p.add_argument("--genes", required=True, type=Path)
    p.add_argument("--cells", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists(): raise FileExistsError(part)
    part.mkdir(parents=True)
    genes = pd.read_csv(args.genes, sep="\t")["gene"].astype(str).tolist()
    cells = pd.read_csv(args.cells, sep="\t")["exact_cell_id"].astype(str).tolist()
    if len(genes) != 37795 or len(cells) != 11950 or len(set(genes)) != len(genes) or len(set(cells)) != len(cells):
        raise ValueError("gene/cell registry gate failed")
    if not set(TARGETS).issubset(genes): raise ValueError("target expression genes missing")
    with gzip.open(args.matrix, "rb") as handle:
        matrix = mmread(handle).tocsc()
    if matrix.shape != (len(genes), len(cells)): raise ValueError("matrix shape gate failed")
    if matrix.data.size and (not np.isfinite(matrix.data).all() or (matrix.data < 0).any() or not np.equal(matrix.data, np.floor(matrix.data)).all()):
        raise ValueError("raw count safety gate failed")
    totals = np.asarray(matrix.sum(axis=0)).ravel()
    if (totals <= 0).any(): raise ValueError("zero library size")
    out = pd.DataFrame({"exact_cell_id": cells})
    qc = []
    for gene in TARGETS:
        raw = matrix[genes.index(gene), :].toarray().ravel().astype(float)
        out[f"{gene}_raw_count"] = raw
        out[f"{gene}_log1p_cp10k"] = np.log1p(raw / totals * 10000.0)
        qc.append({"gene": gene, "detected_cells": int((raw > 0).sum()), "total_counts": float(raw.sum())})
    out.to_csv(part / "m13_tymp_rrm2b_cell_expression_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(qc).to_csv(part / "m13_expression_qc_v1.tsv", sep="\t", index=False)
    status = {"status":"PASS_EXPRESSION_PREP", "n_cells":len(cells), "n_genes":len(genes), "targets":TARGETS}
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output)


if __name__ == "__main__": main()

