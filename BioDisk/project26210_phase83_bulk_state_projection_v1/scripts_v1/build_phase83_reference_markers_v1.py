#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_STATES = [
    "B cell", "C1QC-state Macrophage", "Dendritic cell", "Endothelial",
    "Epithelial-Normal", "Malignant", "Mast cell", "Monocyte", "NK cell",
    "Other Macrophage", "SMC/Pericyte", "SPP1-state Macrophage", "T cell",
]
EXCLUDED_EXACT = {"TYMP", "SPP1"}
EXCLUDED_PREFIX = ("MT-", "RPL", "RPS", "HBA", "HBB")
TOP_PER_STATE = 60


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", required=True, type=Path)
    parser.add_argument("--preregistration", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    signature = pd.read_csv(args.signature, sep="\t", index_col=0, compression="gzip")
    signature.index = signature.index.astype(str).str.upper()
    if list(signature.columns) != EXPECTED_STATES:
        raise RuntimeError("Unexpected Phase17 states")
    log_signature = np.log1p(signature.astype(float))
    excluded = pd.Series(
        signature.index.isin(EXCLUDED_EXACT)
        | signature.index.str.startswith(EXCLUDED_PREFIX),
        index=signature.index,
    )
    specificity = pd.DataFrame(index=signature.index, columns=signature.columns, dtype=float)
    rows = []
    for state in signature.columns:
        other = log_signature.drop(columns=state).mean(axis=1)
        values = log_signature[state] - other
        specificity[state] = values
        eligible = values.loc[~excluded & (values > 0)]
        top = eligible.nlargest(TOP_PER_STATE)
        if len(top) != TOP_PER_STATE:
            raise RuntimeError(f"Too few positive markers for {state}: {len(top)}")
        for rank, (gene, value) in enumerate(top.items(), start=1):
            rows.append({
                "state": state,
                "rank": rank,
                "gene": gene,
                "specificity": float(value),
            })
    markers = pd.DataFrame(rows)
    union = pd.Index(sorted(markers["gene"].unique()))
    reference = specificity.loc[union].clip(lower=0)
    norms = np.sqrt((reference ** 2).sum(axis=0))
    if (norms <= 0).any():
        raise RuntimeError("Zero reference column")
    reference = reference / norms
    reference.insert(0, "gene", reference.index)
    markers.to_csv(args.output / "phase83_marker_members_v1.tsv", sep="\t", index=False)
    reference.to_csv(
        args.output / "phase83_reference_specificity_v1.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    status = {
        "status": "PASS",
        "states": len(EXPECTED_STATES),
        "top_per_state": TOP_PER_STATE,
        "marker_rows": len(markers),
        "unique_marker_genes": len(union),
        "tymp_in_markers": "TYMP" in union,
        "spp1_in_markers": "SPP1" in union,
        "signature_sha256": sha256(args.signature),
        "preregistration_sha256": sha256(args.preregistration),
    }
    (args.output / "PHASE83_REFERENCE_STATUS_v1.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    if status["tymp_in_markers"] or status["spp1_in_markers"]:
        raise RuntimeError("Excluded target leaked into markers")


if __name__ == "__main__":
    main()
