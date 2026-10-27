#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
from pathlib import Path

import anndata as ad
import pandas as pd


def version(name: str) -> str:
    if importlib.util.find_spec(name) is None:
        return "NOT_INSTALLED"
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "IMPORTABLE_VERSION_UNKNOWN"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--bridge", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)

    assignments = pd.read_csv(args.assignments, sep="\t")
    bridge = pd.read_csv(args.bridge, sep="\t")
    obj = ad.read_h5ad(args.h5ad, backed="r")
    genes = set(map(str, obj.var_names))
    required_genes = ["TYMP", "HIF1A", "MST1", "MST1R", "ICAM5", "ITGAM", "ITGB2"]
    required_obs = ["curated_patient_id", "study", "total_counts_raw", "eligible_for_subclone_selection_v1"]
    malignant_states = {"F1-high Malignant", "F1-low Malignant"}
    payload = {
        "status": "PASS_C01_C04_PREFLIGHT",
        "h5ad_shape": list(obj.shape),
        "h5ad_layers": list(obj.layers.keys()),
        "h5ad_obsm": list(obj.obsm.keys()),
        "required_obs_present": {field: field in obj.obs.columns for field in required_obs},
        "required_genes_present": {gene: gene in genes for gene in required_genes},
        "assignments_rows": int(len(assignments)),
        "assignments_unique_cells": int(assignments["exact_cell_id"].nunique()),
        "assignments_patients": int(assignments["patient"].nunique()),
        "assignments_studies": int(assignments["study_v1"].nunique()),
        "malignant_assignment_cells": int(assignments["state_v2"].isin(malignant_states).sum()),
        "bridge_rows": int(len(bridge)),
        "bridge_unique_cells": int(bridge["exact_cell_id"].nunique()),
        "bridge_patients": int(bridge["patient"].nunique()),
        "bridge_tymp_high_cells": int(bridge["TYMP_high"].astype(bool).sum()),
        "bridge_cells_in_assignments": int(bridge["exact_cell_id"].isin(set(assignments["exact_cell_id"])).sum()),
        "python_packages": {name: version(name) for name in ["anndata", "scanpy", "liana", "decoupler", "omnipath"]},
    }
    failures = []
    if obj.shape != (258120, 37795): failures.append("h5ad_shape")
    if "counts" not in obj.layers: failures.append("counts_layer")
    if not all(payload["required_obs_present"].values()): failures.append("obs_fields")
    if not all(payload["required_genes_present"].values()): failures.append("target_genes")
    if payload["assignments_rows"] != payload["assignments_unique_cells"]: failures.append("assignment_duplicates")
    if payload["assignments_patients"] != 40 or payload["assignments_studies"] != 5: failures.append("assignment_scope")
    if payload["malignant_assignment_cells"] != 11921: failures.append("malignant_scope")
    if payload["bridge_rows"] != payload["bridge_unique_cells"] or payload["bridge_rows"] != 11921: failures.append("bridge_scope")
    if payload["bridge_cells_in_assignments"] != 11921: failures.append("bridge_assignment_join")
    if failures:
        payload["status"] = "FAIL_C01_C04_PREFLIGHT"
    payload["failures"] = failures
    (args.output / "status.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    assignments.groupby(["study_v1", "state_v2"], observed=True).agg(
        n_cells=("exact_cell_id", "size"), n_patients=("patient", "nunique")
    ).reset_index().to_csv(args.output / "assignment_scope_v1.tsv", sep="\t", index=False, lineterminator="\n")
    if failures:
        raise RuntimeError(f"Preflight failed: {failures}")


if __name__ == "__main__":
    main()
