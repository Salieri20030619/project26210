#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


GENES = ["CA9", "VEGFA", "ADM", "ANGPTL4", "LDHA"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--cells", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(part)
    part.mkdir(parents=True)

    gene_table = pd.read_csv(args.genes, sep="\t")
    cell_table = pd.read_csv(args.cells, sep="\t")
    if "gene" not in gene_table.columns or "exact_cell_id" not in cell_table.columns:
        raise ValueError("unexpected gene/cell registry schema")
    gene_names = gene_table["gene"].astype(str).tolist()
    cell_names = cell_table["exact_cell_id"].astype(str).tolist()
    if len(gene_names) != len(set(gene_names)) or len(cell_names) != len(set(cell_names)):
        raise ValueError("gene/cell keys are not unique")
    missing = sorted(set(GENES + ["HIF1A"]) - set(gene_names))
    if missing:
        raise ValueError(f"missing genes: {missing}")
    with gzip.open(args.matrix, "rb") as handle:
        matrix = mmread(handle).tocsc()
    if matrix.shape != (len(gene_names), len(cell_names)):
        raise ValueError(f"matrix shape mismatch: {matrix.shape}")
    if matrix.data.size and (not np.isfinite(matrix.data).all() or (matrix.data < 0).any() or not np.equal(matrix.data, np.floor(matrix.data)).all()):
        raise ValueError("raw-count safety gate failed")
    totals = np.asarray(matrix.sum(axis=0)).ravel()
    if (totals <= 0).any():
        raise ValueError("zero library size")
    positions = {gene: gene_names.index(gene) for gene in GENES + ["HIF1A"]}
    raw = matrix[[positions[g] for g in GENES], :].toarray().astype(float)
    logcp10k = np.log1p(raw / totals[None, :] * 10000.0)
    means = logcp10k.mean(axis=1)
    sds = logcp10k.std(axis=1, ddof=1)
    if (sds <= 0).any():
        raise ValueError("constant signature gene")
    z = (logcp10k - means[:, None]) / sds[:, None]
    hif_raw = matrix[positions["HIF1A"], :].toarray().ravel().astype(float)
    hif_log = np.log1p(hif_raw / totals * 10000.0)
    scores = pd.DataFrame({
        "exact_cell_id": cell_names,
        "hif1a_downstream_score": z.mean(axis=0),
        "HIF1A_raw_count": hif_raw,
        "HIF1A_log1p_cp10k": hif_log,
        "library_size": totals,
    })
    scores.to_csv(part / "cell_hif1a_downstream_score_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame({"gene": GENES, "mean_log1p_cp10k": means, "sd_log1p_cp10k": sds, "detected_cells": (raw > 0).sum(axis=1)}).to_csv(part / "signature_gene_qc_v1.tsv", sep="\t", index=False)
    pd.DataFrame([{"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in (args.matrix, args.genes, args.cells)]).to_csv(part / "INPUT_MANIFEST_v1.tsv", sep="\t", index=False)
    status = {"status": "PASS", "n_cells": len(cell_names), "n_genes": len(gene_names), "signature_genes": GENES, "HIF1A_detected_cells": int((hif_raw > 0).sum())}
    (part / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output)


if __name__ == "__main__":
    main()
