#!/usr/bin/env python3
"""Audit the inputs available for Part 1 de novo cell annotation.

This script is read-only with respect to source objects.  It deliberately treats
legacy annotations as comparison fields, not as labels to copy into the rebuilt
object.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT = Path("/mnt/BioDisk/ludingyu/project26210")
CURRENT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/scvi_integration_v1/part1_scvi_integrated_v1.h5ad"
LEGACY = PROJECT / "object/01_annotation_adata_annotated.h5ad"
MEMORY_GATE = PROJECT / "check_memory_gate.sh"
DEFAULT_OUT = PROJECT / "results/part1_rebuild_v1/02_integration_annotation/annotation_input_audit_v1"


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def run_memory_gate(stage: str, records: list[dict[str, object]]) -> None:
    proc = subprocess.run([str(MEMORY_GATE)], text=True, capture_output=True)
    line = (proc.stdout or proc.stderr).strip()
    records.append({"stage": stage, "returncode": proc.returncode, "output": line})
    print(f"memory_gate[{stage}]: {line}", flush=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Memory gate failed before {stage}: {line}")


def summarize_obs(source: str, obs: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    rows: list[dict[str, object]] = []
    candidates: list[str] = []
    pattern = re.compile(r"cell|type|annot|label|conf|cluster|leiden|prediction", re.I)
    for col in obs.columns:
        series = obs[col]
        nunique = int(series.nunique(dropna=True))
        is_candidate = bool(pattern.search(str(col))) and 1 < nunique <= 500
        if is_candidate:
            candidates.append(str(col))
        examples = [str(x) for x in series.dropna().astype(str).unique()[:20]]
        rows.append(
            {
                "source": source,
                "field": str(col),
                "dtype": str(series.dtype),
                "n_rows": len(series),
                "n_missing": int(series.isna().sum()),
                "n_unique_nonmissing": nunique,
                "categorical": bool(isinstance(series.dtype, pd.CategoricalDtype)),
                "annotation_candidate": is_candidate,
                "first_levels_or_values": " | ".join(examples),
            }
        )
    return pd.DataFrame(rows), candidates


def normalized_entropy(counts: np.ndarray) -> float:
    counts = counts[counts > 0]
    if len(counts) <= 1:
        return 0.0
    prob = counts / counts.sum()
    return float(-(prob * np.log(prob)).sum() / math.log(len(prob)))


def package_rows() -> list[dict[str, object]]:
    rows = []
    for name in ["anndata", "scanpy", "scvi-tools", "scvi", "celltypist", "scikit-learn", "decoupler", "pandas", "numpy"]:
        module = {"scvi-tools": "scvi", "scikit-learn": "sklearn"}.get(name, name)
        available = importlib.util.find_spec(module) is not None
        version = ""
        if available:
            for dist in [name, module]:
                try:
                    version = importlib.metadata.version(dist)
                    break
                except importlib.metadata.PackageNotFoundError:
                    pass
        rows.append({"capability": "python_package", "name": name, "available": available, "version_or_path": version})

    model_roots = [
        Path("/home/wangmingyuan/.celltypist"),
        Path("/mnt/BioDisk/wangmingyuan/.celltypist"),
        PROJECT / "reference",
        PROJECT / "references",
    ]
    model_files: list[Path] = []
    for root in model_roots:
        if root.exists():
            model_files.extend(p for p in root.rglob("*.pkl") if p.is_file())
    if model_files:
        for path in sorted(set(model_files)):
            rows.append({"capability": "local_celltypist_model", "name": path.name, "available": True, "version_or_path": str(path)})
    else:
        rows.append({"capability": "local_celltypist_model", "name": "*.pkl", "available": False, "version_or_path": ""})
    return rows


def marker_inventory() -> tuple[pd.DataFrame, pd.DataFrame]:
    patterns = [
        "results/celltype_marker_dotplot/*.tsv",
        "results/celltype_marker_dotplot/*.csv",
        "results/celltype_marker_violin/*.csv",
        "results/*marker*/*.csv",
        "results/*marker*/*.tsv",
    ]
    files: list[Path] = []
    for pattern in patterns:
        files.extend(PROJECT.glob(pattern))
    files = sorted(set(p for p in files if p.is_file()))
    inventory: list[dict[str, object]] = []
    genes: list[dict[str, object]] = []
    gene_col_re = re.compile(r"^(gene|genes|marker|markers|feature|symbol)$", re.I)
    for path in files:
        rel = path.relative_to(PROJECT)
        inventory.append(
            {
                "path": str(rel),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
        try:
            sep = "\t" if path.suffix.lower() == ".tsv" else ","
            frame = pd.read_csv(path, sep=sep, nrows=100000)
            for col in frame.columns:
                if gene_col_re.match(str(col)):
                    for gene in frame[col].dropna().astype(str).unique():
                        genes.append({"path": str(rel), "column": str(col), "gene": gene})
        except Exception as exc:  # inventory must survive non-tabular legacy files
            inventory[-1]["parse_error"] = repr(exc)
    return pd.DataFrame(inventory), pd.DataFrame(genes)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--force", action="store_true", help="Replace only an incomplete audit directory.")
    args = parser.parse_args()
    outdir: Path = args.outdir
    complete = outdir / "RUN_COMPLETE"
    if complete.exists():
        raise SystemExit(f"Refusing to overwrite completed audit: {complete}")
    if outdir.exists() and any(outdir.iterdir()) and not args.force:
        raise SystemExit(f"Output exists and is non-empty; use --force only for an incomplete run: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    started = datetime.now().astimezone().isoformat()
    memory_records: list[dict[str, object]] = []
    run_memory_gate("before_current_h5ad", memory_records)
    current = ad.read_h5ad(CURRENT, backed="r")
    current_obs = current.obs.copy()
    current_var_names = pd.Index(current.var_names.astype(str))
    current_schema, current_candidates = summarize_obs("current_scvi", current_obs)

    run_memory_gate("before_legacy_h5ad", memory_records)
    legacy = ad.read_h5ad(LEGACY, backed="r")
    legacy_obs = legacy.obs.copy()
    legacy_var_names = pd.Index(legacy.var_names.astype(str))
    legacy_schema, legacy_candidates = summarize_obs("legacy_reference", legacy_obs)

    positions = legacy_obs.index.get_indexer(current_obs.index)
    if (positions < 0).any():
        missing = int((positions < 0).sum())
        raise RuntimeError(f"Current cells absent from legacy object: {missing}")
    if len(np.unique(positions)) != len(positions):
        raise RuntimeError("Current cells map non-uniquely into the legacy object")
    legacy_retained = legacy_obs.iloc[positions].copy()
    legacy_retained.index = current_obs.index.copy()

    pd.concat([current_schema, legacy_schema], ignore_index=True).to_csv(outdir / "h5ad_obs_schema_v1.tsv", sep="\t", index=False)

    candidate_rows = [
        {"source": "current_scvi", "field": field} for field in current_candidates
    ] + [{"source": "legacy_reference", "field": field} for field in legacy_candidates]
    pd.DataFrame(candidate_rows).to_csv(outdir / "annotation_candidate_fields_v1.tsv", sep="\t", index=False)

    count_rows: list[dict[str, object]] = []
    for field in legacy_candidates:
        for value, count in legacy_retained[field].astype("string").fillna("<NA>").value_counts(dropna=False).items():
            count_rows.append({"field": field, "value": str(value), "n_cells": int(count), "fraction": float(count / len(legacy_retained))})
    pd.DataFrame(count_rows).to_csv(outdir / "legacy_candidate_label_counts_retained_v1.tsv", sep="\t", index=False)

    leiden_fields = [c for c in current_obs.columns if re.search(r"leiden", str(c), re.I)]
    cluster_count_rows: list[dict[str, object]] = []
    for field in leiden_fields:
        counts = current_obs[field].astype("string").fillna("<NA>").value_counts(dropna=False)
        for cluster, count in counts.items():
            cluster_count_rows.append({"cluster_field": field, "cluster": str(cluster), "n_cells": int(count), "fraction": float(count / len(current_obs))})
    pd.DataFrame(cluster_count_rows).to_csv(outdir / "current_cluster_counts_v1.tsv", sep="\t", index=False)

    comparison_rows: list[dict[str, object]] = []
    concordance_rows: list[dict[str, object]] = []
    for cluster_field in leiden_fields:
        clusters = current_obs[cluster_field].astype("string").fillna("<NA>")
        for label_field in legacy_candidates:
            labels = legacy_retained[label_field].astype("string").fillna("<NA>")
            table = pd.crosstab(clusters, labels, dropna=False)
            for cluster in table.index:
                row = table.loc[cluster]
                total = int(row.sum())
                top_label = str(row.idxmax())
                top_n = int(row.max())
                concordance_rows.append(
                    {
                        "cluster_field": cluster_field,
                        "cluster": str(cluster),
                        "legacy_field": label_field,
                        "n_cells": total,
                        "top_legacy_value": top_label,
                        "top_n": top_n,
                        "top_fraction": float(top_n / total),
                        "normalized_entropy": normalized_entropy(row.to_numpy(dtype=float)),
                    }
                )
                for label, count in row[row > 0].items():
                    comparison_rows.append(
                        {
                            "cluster_field": cluster_field,
                            "cluster": str(cluster),
                            "legacy_field": label_field,
                            "legacy_value": str(label),
                            "n_cells": int(count),
                            "fraction_within_cluster": float(count / total),
                        }
                    )
    pd.DataFrame(comparison_rows).to_csv(outdir / "legacy_labels_by_current_cluster_v1.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(concordance_rows).to_csv(outdir / "cluster_legacy_concordance_v1.tsv", sep="\t", index=False)

    inventory, marker_genes = marker_inventory()
    inventory.to_csv(outdir / "legacy_marker_resource_inventory_v1.tsv", sep="\t", index=False)
    if len(marker_genes):
        marker_genes["present_current"] = marker_genes["gene"].isin(current_var_names)
        marker_genes["present_legacy"] = marker_genes["gene"].isin(legacy_var_names)
    marker_genes.to_csv(outdir / "legacy_marker_gene_presence_v1.tsv", sep="\t", index=False)
    pd.DataFrame(package_rows()).to_csv(outdir / "software_capabilities_v1.tsv", sep="\t", index=False)

    object_rows = [
        {
            "source": "current_scvi",
            "path": str(CURRENT),
            "n_obs": int(current.n_obs),
            "n_vars": int(current.n_vars),
            "obs_names_unique": bool(current.obs_names.is_unique),
            "var_names_unique": bool(current.var_names.is_unique),
            "x_backend": type(current.X).__name__,
            "layers": " | ".join(current.layers.keys()),
            "obsm": " | ".join(current.obsm.keys()),
            "obsp": " | ".join(current.obsp.keys()),
            "uns": " | ".join(current.uns.keys()),
        },
        {
            "source": "legacy_reference",
            "path": str(LEGACY),
            "n_obs": int(legacy.n_obs),
            "n_vars": int(legacy.n_vars),
            "obs_names_unique": bool(legacy.obs_names.is_unique),
            "var_names_unique": bool(legacy.var_names.is_unique),
            "x_backend": type(legacy.X).__name__,
            "layers": " | ".join(legacy.layers.keys()),
            "obsm": " | ".join(legacy.obsm.keys()),
            "obsp": " | ".join(legacy.obsp.keys()),
            "uns": " | ".join(legacy.uns.keys()),
        },
    ]
    pd.DataFrame(object_rows).to_csv(outdir / "h5ad_structure_v1.tsv", sep="\t", index=False)

    genes_equal = bool(current_var_names.equals(legacy_var_names))
    metadata = {
        "run_id": "part1_annotation_input_audit_v1",
        "started_at": started,
        "completed_at": datetime.now().astimezone().isoformat(),
        "python_executable": sys.executable,
        "command": " ".join(sys.argv),
        "current_h5ad": str(CURRENT),
        "legacy_h5ad": str(LEGACY),
        "current_shape": [int(current.n_obs), int(current.n_vars)],
        "legacy_shape": [int(legacy.n_obs), int(legacy.n_vars)],
        "current_cells_all_in_legacy": True,
        "current_cell_order_matches_legacy_subset": bool(np.array_equal(positions, np.sort(positions))),
        "gene_order_exact_match": genes_equal,
        "n_current_candidate_fields": len(current_candidates),
        "n_legacy_candidate_fields": len(legacy_candidates),
        "memory_gate": memory_records,
        "legacy_use_policy": "comparison_only_not_copied_as_rebuilt_annotation",
    }
    (outdir / "run_metadata_v1.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (outdir / "RUN_COMPLETE").write_text(metadata["completed_at"] + "\n", encoding="utf-8")
    current.file.close()
    legacy.file.close()
    print(json.dumps(metadata, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
