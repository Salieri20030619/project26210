#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import h5py
import numpy as np


def scan_numeric(values):
    non_integer = negative = non_finite = 0
    minimum = None
    maximum = None
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            non_finite += 1
            continue
        if not math.isfinite(number):
            non_finite += 1
            continue
        if number < 0:
            negative += 1
        if not number.is_integer():
            non_integer += 1
        minimum = number if minimum is None else min(minimum, number)
        maximum = number if maximum is None else max(maximum, number)
    return non_integer, negative, non_finite, minimum, maximum


def audit_dense_csv(path):
    features = nnz = non_integer = negative = non_finite = 0
    minimum = maximum = None
    row_width_ok = True
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        if header[0] not in {"Gene_ID", "Gene-ID"} or header[1] != "Symbol":
            raise ValueError(f"unexpected CSV annotation columns: {header[:2]}")
        cells = len(header) - 2
        if cells <= 0:
            raise ValueError("CSV header has no cell columns")
        for line in handle:
            features += 1
            parts = line.rstrip("\r\n").split(",", 2)
            if len(parts) != 3:
                row_width_ok = False
                continue
            values = np.fromstring(parts[2], dtype=float, sep=",")
            if values.size != cells:
                row_width_ok = False
                continue
            finite = np.isfinite(values)
            non_finite += int((~finite).sum())
            clean = values[finite]
            negative += int((clean < 0).sum())
            non_integer += int((clean != np.floor(clean)).sum())
            nnz += int(np.count_nonzero(clean))
            if clean.size:
                lo = float(clean.min())
                hi = float(clean.max())
                minimum = lo if minimum is None else min(minimum, lo)
                maximum = hi if maximum is None else max(maximum, hi)
    checks = row_width_ok and features > 0
    return features, cells, nnz, minimum, maximum, non_integer, negative, non_finite, checks


def scan_h5_data(dataset, chunk_size=1_000_000):
    non_integer = negative = non_finite = 0
    minimum = maximum = None
    for start in range(0, dataset.shape[0], chunk_size):
        block = np.asarray(dataset[start : start + chunk_size])
        finite = np.isfinite(block)
        non_finite += int((~finite).sum())
        clean = block[finite]
        if clean.size:
            negative += int((clean < 0).sum())
            if not np.issubdtype(clean.dtype, np.integer):
                non_integer += int((clean != np.floor(clean)).sum())
            lo = float(clean.min())
            hi = float(clean.max())
            minimum = lo if minimum is None else min(minimum, lo)
            maximum = hi if maximum is None else max(maximum, hi)
    return non_integer, negative, non_finite, minimum, maximum


def audit_h5_group(group):
    required = {"data", "indices", "indptr", "shape", "barcodes"}
    if not required.issubset(group.keys()):
        raise ValueError(f"missing H5 datasets: {sorted(required - set(group.keys()))}")
    shape = tuple(int(x) for x in np.asarray(group["shape"]))
    if len(shape) != 2:
        raise ValueError(f"invalid sparse shape: {shape}")
    features, cells = shape
    data = group["data"]
    indices = np.asarray(group["indices"])
    indptr = np.asarray(group["indptr"])
    nnz = int(data.shape[0])
    feature_len = None
    if "features" in group:
        fg = group["features"]
        for key in ("id", "name"):
            if key in fg:
                feature_len = int(fg[key].shape[0])
                break
    else:
        for key in ("genes", "gene_names"):
            if key in group:
                feature_len = int(group[key].shape[0])
                break
    structure_ok = (
        features > 0
        and cells > 0
        and indices.shape[0] == nnz
        and indptr.shape[0] == cells + 1
        and int(indptr[0]) == 0
        and int(indptr[-1]) == nnz
        and bool(np.all(indptr[1:] >= indptr[:-1]))
        and (nnz == 0 or (int(indices.min()) >= 0 and int(indices.max()) < features))
        and int(group["barcodes"].shape[0]) == cells
        and feature_len == features
    )
    ni, neg, nf, lo, hi = scan_h5_data(data)
    return features, cells, nnz, lo, hi, ni, neg, nf, structure_ok


def audit_h5(path):
    with h5py.File(path, "r") as handle:
        if "matrix" in handle:
            return audit_h5_group(handle["matrix"])
        groups = [handle[key] for key in handle.keys() if isinstance(handle[key], h5py.Group)]
        groups = [group for group in groups if "data" in group]
        if not groups:
            raise ValueError("no 10x sparse matrix group found")
        results = [audit_h5_group(group) for group in groups]
        cells = {result[1] for result in results}
        if len(cells) != 1:
            raise ValueError("legacy H5 genome groups have inconsistent barcode counts")
        return (
            sum(result[0] for result in results),
            results[0][1],
            sum(result[2] for result in results),
            min(result[3] for result in results if result[3] is not None),
            max(result[4] for result in results if result[4] is not None),
            sum(result[5] for result in results),
            sum(result[6] for result in results),
            sum(result[7] for result in results),
            all(result[8] for result in results),
        )


def count_gzip_lines(path):
    count = 0
    with gzip.open(path, "rb") as handle:
        for _ in handle:
            count += 1
    return count


def audit_mtx(matrix_path, features_path, barcodes_path):
    dimensions = None
    observed = non_integer = negative = non_finite = 0
    minimum = maximum = None
    coordinates_ok = True
    with gzip.open(matrix_path, "rt", encoding="utf-8") as handle:
        banner = handle.readline().strip()
        if not banner.lower().startswith("%%matrixmarket matrix coordinate"):
            raise ValueError("unsupported MatrixMarket banner")
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.split()
            if dimensions is None:
                if len(parts) != 3:
                    raise ValueError("invalid MatrixMarket dimensions")
                dimensions = tuple(int(x) for x in parts)
                continue
            if len(parts) < 3:
                coordinates_ok = False
                continue
            row_index, column_index = int(parts[0]), int(parts[1])
            observed += 1
            if not (1 <= row_index <= dimensions[0] and 1 <= column_index <= dimensions[1]):
                coordinates_ok = False
            ni, neg, nf, lo, hi = scan_numeric([parts[2]])
            non_integer += ni
            negative += neg
            non_finite += nf
            if lo is not None:
                minimum = lo if minimum is None else min(minimum, lo)
                maximum = hi if maximum is None else max(maximum, hi)
    if dimensions is None:
        raise ValueError("MatrixMarket dimensions not found")
    features, cells, declared = dimensions
    structure_ok = (
        features > 0
        and cells > 0
        and declared == observed
        and count_gzip_lines(features_path) == features
        and count_gzip_lines(barcodes_path) == cells
        and coordinates_ok
    )
    return features, cells, observed, minimum, maximum, non_integer, negative, non_finite, structure_ok


def format_number(value):
    if value is None:
        return "NA"
    return str(int(value)) if float(value).is_integer() else repr(float(value))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract-root", type=Path, required=True)
    parser.add_argument("--scope", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with args.scope.open("r", encoding="utf-8-sig", newline="") as handle:
        scope = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in scope if row["study"] in {"GSE156632", "GSE159115", "GSE242299"}]
    if len(selected) != 57:
        raise RuntimeError(f"expected 57 public archive libraries, observed {len(selected)}")

    audit_rows = []
    for row in selected:
        study = row["study"]
        members = row["archive_members"].split("|")
        paths = [args.extract_root / study / member for member in members]
        if not all(path.is_file() for path in paths):
            missing = [str(path) for path in paths if not path.is_file()]
            raise FileNotFoundError(";".join(missing))
        error = ""
        try:
            if study == "GSE156632":
                result = audit_dense_csv(paths[0])
                source_format = "dense_csv_gzip"
            elif study == "GSE159115":
                result = audit_h5(paths[0])
                source_format = "10x_h5"
            else:
                by_suffix = {path.name.rsplit("_", 1)[-1]: path for path in paths}
                matrix_path = next(path for path in paths if path.name.endswith("_matrix.mtx.gz"))
                features_path = next(path for path in paths if path.name.endswith("_features.tsv.gz"))
                barcodes_path = next(path for path in paths if path.name.endswith("_barcodes.tsv.gz"))
                result = audit_mtx(matrix_path, features_path, barcodes_path)
                source_format = "10x_mtx_gzip"
            features, cells, nnz, lo, hi, ni, neg, nf, structure_ok = result
            status = "PASS" if structure_ok and ni == 0 and neg == 0 and nf == 0 else "FAIL"
        except Exception as exc:
            features = cells = nnz = lo = hi = "NA"
            ni = neg = nf = "NA"
            structure_ok = False
            status = "FAIL"
            error = f"{type(exc).__name__}:{exc}"
            source_format = "unknown"
        audit_rows.append(
            {
                "study": study,
                "geo_accession": row["geo_accession"],
                "library_id": row["library_id"],
                "inclusion_status": row["inclusion_status"],
                "source_format": source_format,
                "source_files": "|".join(path.name for path in paths),
                "features": features,
                "cells": cells,
                "nnz": nnz,
                "min_count": format_number(lo) if lo != "NA" else "NA",
                "max_count": format_number(hi) if hi != "NA" else "NA",
                "non_integer_values": ni,
                "negative_values": neg,
                "non_finite_values": nf,
                "structural_checks": str(bool(structure_ok)).lower(),
                "status": status,
                "error": error,
            }
        )

    args.output.mkdir(parents=True, exist_ok=True)
    audit_path = args.output / "raw_count_structure_audit_v3.tsv"
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)

    summary_rows = []
    for study in ("GSE156632", "GSE159115", "GSE242299"):
        rows = [row for row in audit_rows if row["study"] == study]
        summary_rows.append(
            {
                "study": study,
                "audited_libraries": len(rows),
                "included_libraries": sum(row["inclusion_status"] == "include" for row in rows),
                "excluded_libraries": sum(row["inclusion_status"] == "exclude" for row in rows),
                "pass_libraries": sum(row["status"] == "PASS" for row in rows),
                "fail_libraries": sum(row["status"] == "FAIL" for row in rows),
                "status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
            }
        )
    summary_path = args.output / "raw_count_structure_summary_v3.tsv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    if len(audit_rows) != 57 or sum(row["inclusion_status"] == "include" for row in audit_rows) != 54:
        raise RuntimeError("scope count mismatch after audit")
    expected_failures = {("GSE242299", "GSM7757833"), ("GSE242299", "GSM7757855")}
    observed_failures = {
        (row["study"], row["geo_accession"])
        for row in audit_rows
        if row["status"] == "FAIL"
    }
    if observed_failures != expected_failures:
        raise RuntimeError(
            f"unexpected failure set: observed={sorted(observed_failures)} expected={sorted(expected_failures)}"
        )
    status_path = args.output / "PUBLIC_GEO_RAW_COUNT_STRUCTURE_BLOCKED_v3"
    status_path.write_text(
        "status=BLOCKED_SOURCE_DEFECTS\n"
        "audited_public_libraries=57\n"
        "pass_libraries=55\n"
        "fail_libraries=2\n"
        "failed_libraries=GSM7757833:T2_1|GSM7757855:N3_1\n"
        "downstream_started=false\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
