#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


def scan_values(dataset, block_size=1_000_000):
    finite_all = nonnegative_all = integer_all = True
    minimum = maximum = None
    if dataset.ndim == 1:
        iterator = (dataset[start : start + block_size] for start in range(0, dataset.shape[0], block_size))
    elif dataset.ndim == 2:
        rows = max(1, block_size // max(1, dataset.shape[1]))
        iterator = (dataset[start : start + rows, :] for start in range(0, dataset.shape[0], rows))
    else:
        raise ValueError(f"unsupported matrix rank: {dataset.ndim}")
    for raw_block in iterator:
        block = np.asarray(raw_block)
        finite = np.isfinite(block)
        if not bool(finite.all()):
            finite_all = False
        clean = block[finite]
        if clean.size == 0:
            continue
        if bool((clean < 0).any()):
            nonnegative_all = False
        if not np.issubdtype(clean.dtype, np.integer) and bool((clean != np.floor(clean)).any()):
            integer_all = False
        lo = float(clean.min())
        hi = float(clean.max())
        minimum = lo if minimum is None else min(minimum, lo)
        maximum = hi if maximum is None else max(maximum, hi)
    return finite_all, nonnegative_all, integer_all, minimum, maximum


def matrix_candidate(handle, path):
    node = handle[path]
    if isinstance(node, h5py.Group) and {"data", "indices", "indptr"}.issubset(node.keys()):
        shape_attr = node.attrs.get("shape")
        if shape_attr is None:
            raise ValueError(f"sparse matrix lacks shape attribute: {path}")
        shape = tuple(int(x) for x in shape_attr)
        data = node["data"]
        indices = np.asarray(node["indices"])
        indptr = np.asarray(node["indptr"])
        stored = int(data.shape[0])
        encoding = node.attrs.get("encoding-type", "sparse")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        is_csr = "csr" in str(encoding).lower()
        pointer_target = shape[0] + 1 if is_csr else shape[1] + 1
        index_limit = shape[1] if is_csr else shape[0]
        structural = (
            len(shape) == 2
            and indices.shape[0] == stored
            and indptr.shape[0] == pointer_target
            and int(indptr[0]) == 0
            and int(indptr[-1]) == stored
            and bool(np.all(indptr[1:] >= indptr[:-1]))
            and (stored == 0 or (int(indices.min()) >= 0 and int(indices.max()) < index_limit))
        )
        finite, nonnegative, integer, minimum, maximum = scan_values(data)
        storage = str(encoding)
        dtype = str(data.dtype)
    elif isinstance(node, h5py.Dataset) and node.ndim == 2:
        shape = tuple(int(x) for x in node.shape)
        stored = int(np.prod(shape))
        structural = all(value > 0 for value in shape)
        finite, nonnegative, integer, minimum, maximum = scan_values(node)
        storage = "dense"
        dtype = str(node.dtype)
    else:
        raise ValueError(f"unsupported matrix encoding at {path}")
    return {
        "matrix_path": path,
        "storage": storage,
        "shape_rows": shape[0],
        "shape_columns": shape[1],
        "stored_values": stored,
        "dtype": dtype,
        "minimum": "NA" if minimum is None else minimum,
        "maximum": "NA" if maximum is None else maximum,
        "finite_all": str(finite).lower(),
        "nonnegative_all": str(nonnegative).lower(),
        "integer_all": str(integer).lower(),
        "structural_checks": str(structural).lower(),
        "raw_count_candidate": str(structural and finite and nonnegative and integer).lower(),
    }


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--scope", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.scope.open("r", encoding="utf-8-sig", newline="") as handle:
        scope_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_libraries = sorted(
        row["library_id"]
        for row in scope_rows
        if row["study"] == "GSE242299" and row["inclusion_status"] == "include"
    )
    if len(expected_libraries) != 31 or len(set(expected_libraries)) != 31:
        raise RuntimeError("expected GSE242299 scope does not contain 31 unique libraries")
    expected_set = set(expected_libraries)

    dataset_rows = []
    matrix_paths = []
    with h5py.File(args.h5ad, "r") as handle:
        def visitor(name, node):
            if isinstance(node, h5py.Dataset):
                dataset_rows.append(
                    {
                        "path": "/" + name,
                        "shape": "x".join(str(value) for value in node.shape),
                        "dtype": str(node.dtype),
                        "compression": str(node.compression or "none"),
                    }
                )
        handle.visititems(visitor)
        if "X" in handle:
            matrix_paths.append("/X")
        if "layers" in handle:
            matrix_paths.extend(f"/layers/{key}" for key in handle["layers"].keys())
        if "raw" in handle and "X" in handle["raw"]:
            matrix_paths.append("/raw/X")
        matrix_rows = [matrix_candidate(handle, path) for path in matrix_paths]

    write_tsv(args.output / "h5ad_dataset_inventory_v2.tsv", dataset_rows, ["path", "shape", "dtype", "compression"])
    write_tsv(
        args.output / "h5ad_matrix_candidate_audit_v2.tsv",
        matrix_rows,
        [
            "matrix_path", "storage", "shape_rows", "shape_columns", "stored_values", "dtype",
            "minimum", "maximum", "finite_all", "nonnegative_all", "integer_all",
            "structural_checks", "raw_count_candidate",
        ],
    )

    adata = ad.read_h5ad(args.h5ad, backed="r")
    try:
        n_obs, n_vars = int(adata.n_obs), int(adata.n_vars)
        obs_rows = []
        exact_columns = []
        for column in adata.obs.columns:
            values = adata.obs[column].astype(str)
            unique = set(values.unique())
            exact = unique == expected_set
            if exact:
                exact_columns.append(str(column))
            if len(unique) <= 100 or any(token in str(column).lower() for token in ("library", "sample", "batch", "emulsion", "orig")):
                obs_rows.append(
                    {
                        "column": str(column),
                        "unique_values": len(unique),
                        "exact_31_library_match": str(exact).lower(),
                        "values": "|".join(sorted(unique)[:200]),
                    }
                )
        write_tsv(
            args.output / "h5ad_obs_candidate_columns_v2.tsv",
            obs_rows,
            ["column", "unique_values", "exact_31_library_match", "values"],
        )
        library_rows = []
        if exact_columns:
            library_column = exact_columns[0]
            counts = adata.obs[library_column].astype(str).value_counts()
            library_rows = [
                {"library_id": library, "cells": int(counts.get(library, 0)), "status": "PASS" if counts.get(library, 0) > 0 else "FAIL"}
                for library in expected_libraries
            ]
            write_tsv(args.output / "h5ad_library_cell_counts_v2.tsv", library_rows, ["library_id", "cells", "status"])
        else:
            library_column = "NA"
    finally:
        adata.file.close()

    raw_candidates = [
        row for row in matrix_rows
        if row["raw_count_candidate"] == "true"
        and int(row["shape_rows"]) == n_obs
        and int(row["shape_columns"]) == n_vars
    ]
    dimensions_pass = n_obs == 50236 and n_vars == 33538
    libraries_pass = bool(exact_columns) and len(library_rows) == 31 and all(row["status"] == "PASS" for row in library_rows)
    affected_pass = libraries_pass and all(
        next(row["cells"] for row in library_rows if row["library_id"] == library) > 0
        for library in ("T2_1", "N3_1")
    )
    matrix_pass = any(row["matrix_path"] == "/raw/X" for row in raw_candidates)
    overall = dimensions_pass and libraries_pass and affected_pass and matrix_pass
    metadata = {
        "h5ad": str(args.h5ad),
        "n_obs": n_obs,
        "n_vars": n_vars,
        "expected_libraries": expected_libraries,
        "exact_library_columns": exact_columns,
        "selected_library_column": library_column,
        "raw_count_candidate_paths": [row["matrix_path"] for row in raw_candidates],
        "dimensions_pass": dimensions_pass,
        "libraries_pass": libraries_pass,
        "affected_libraries_pass": affected_pass,
        "matrix_pass": matrix_pass,
        "status": "PASS" if overall else "BLOCKED",
    }
    (args.output / "h5ad_recovery_metadata_v2.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    marker = args.output / ("GSE242299_H5AD_RECOVERY_PASS_v2" if overall else "GSE242299_H5AD_RECOVERY_BLOCKED_v2")
    marker.write_text(
        f"status={metadata['status']}\n"
        f"n_obs={n_obs}\n"
        f"n_vars={n_vars}\n"
        f"library_column={library_column}\n"
        f"library_count={len(library_rows)}\n"
        f"T2_1_cells={next((row['cells'] for row in library_rows if row['library_id'] == 'T2_1'), 0)}\n"
        f"N3_1_cells={next((row['cells'] for row in library_rows if row['library_id'] == 'N3_1'), 0)}\n"
        f"raw_count_candidates={'|'.join(metadata['raw_count_candidate_paths'])}\n",
        encoding="utf-8",
    )
    if not overall:
        raise RuntimeError("official H5AD recovery gate did not pass")


if __name__ == "__main__":
    main()
