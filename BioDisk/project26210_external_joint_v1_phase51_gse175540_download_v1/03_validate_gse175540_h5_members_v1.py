#!/usr/bin/env python3
"""Read-only structural and count-safety preflight for GSE175540 10x H5 members."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import shutil
import tarfile
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np


EXPECTED_SAMPLES = 24
EXPECTED_H5_MEMBERS = 48
H5_SUFFIXES = ("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")


def check_numeric_dataset(ds: h5py.Dataset, block: int = 1_000_000) -> tuple[bool, str, int]:
    total = 0
    for start in range(0, int(ds.shape[0]), block):
        arr = np.asarray(ds[start : start + block])
        if not np.issubdtype(arr.dtype, np.number):
            return False, f"non_numeric_dtype:{arr.dtype}", total
        if not np.all(np.isfinite(arr)):
            return False, "non_finite_count", total
        if np.any(arr < 0):
            return False, "negative_count", total
        if np.issubdtype(arr.dtype, np.floating) and np.any(arr != np.floor(arr)):
            return False, "non_integer_count", total
        total += int(arr.astype(np.int64, copy=False).sum(dtype=np.int64))
    return True, "PASS", total


def validate_h5(path: Path, member_name: str) -> dict[str, object]:
    accession = Path(member_name).name.split("_", 1)[0]
    matrix_type = "filtered" if member_name.endswith(H5_SUFFIXES[0]) else "raw"
    row: dict[str, object] = {
        "sample_accession": accession,
        "matrix_type": matrix_type,
        "archive_member": member_name,
        "status": "FAIL",
        "reason": "not_checked",
        "genes": "",
        "barcodes": "",
        "nnz": "",
        "total_counts": "",
        "count_dtype": "",
    }
    try:
        with h5py.File(path, "r") as handle:
            if "matrix" not in handle:
                row["reason"] = "missing_matrix_group"
                return row
            matrix = handle["matrix"]
            required = ("data", "indices", "indptr", "shape", "barcodes", "features")
            missing = [key for key in required if key not in matrix]
            if missing:
                row["reason"] = "missing:" + ",".join(missing)
                return row
            shape = np.asarray(matrix["shape"][:], dtype=np.int64)
            if shape.size != 2 or np.any(shape <= 0):
                row["reason"] = f"invalid_shape:{shape.tolist()}"
                return row
            genes, barcodes = map(int, shape)
            data = matrix["data"]
            indices = matrix["indices"]
            indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
            nnz = int(data.shape[0])
            if int(indices.shape[0]) != nnz:
                row["reason"] = "data_indices_length_mismatch"
                return row
            if indptr.size != barcodes + 1 or indptr[0] != 0 or indptr[-1] != nnz:
                row["reason"] = "invalid_indptr_boundaries"
                return row
            if np.any(np.diff(indptr) < 0):
                row["reason"] = "non_monotonic_indptr"
                return row
            if int(matrix["barcodes"].shape[0]) != barcodes:
                row["reason"] = "barcode_shape_mismatch"
                return row
            features = matrix["features"]
            if "id" not in features or "name" not in features:
                row["reason"] = "missing_feature_id_or_name"
                return row
            if int(features["id"].shape[0]) != genes or int(features["name"].shape[0]) != genes:
                row["reason"] = "feature_shape_mismatch"
                return row
            for start in range(0, nnz, 1_000_000):
                idx = np.asarray(indices[start : start + 1_000_000], dtype=np.int64)
                if np.any(idx < 0) or np.any(idx >= genes):
                    row["reason"] = "index_out_of_bounds"
                    return row
            safe, reason, total_counts = check_numeric_dataset(data)
            row.update(
                {
                    "genes": genes,
                    "barcodes": barcodes,
                    "nnz": nnz,
                    "total_counts": total_counts,
                    "count_dtype": str(data.dtype),
                    "status": "PASS" if safe else "FAIL",
                    "reason": reason,
                }
            )
            return row
    except Exception as exc:  # keep a complete per-member audit instead of stopping early
        row["reason"] = f"{type(exc).__name__}:{exc}"
        return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()

    archive = args.archive.resolve()
    output_root = args.output_root.resolve()
    final_dir = output_root / "h5_preflight_v1"
    part_dir = output_root / "h5_preflight_v1.part"
    if final_dir.exists() or part_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {final_dir} or {part_dir}")
    if not archive.is_file():
        raise SystemExit(f"Archive not found: {archive}")
    part_dir.mkdir(parents=True)

    rows: list[dict[str, object]] = []
    with tarfile.open(archive, "r") as tf:
        members = [m for m in tf.getmembers() if m.isfile() and m.name.endswith(H5_SUFFIXES)]
        members.sort(key=lambda item: item.name)
        for member in members:
            source = tf.extractfile(member)
            if source is None:
                rows.append({"archive_member": member.name, "status": "FAIL", "reason": "extractfile_none"})
                continue
            with tempfile.NamedTemporaryFile(dir=part_dir, suffix=".h5.part", delete=False) as tmp:
                tmp_path = Path(tmp.name)
                shutil.copyfileobj(source, tmp, length=1024 * 1024)
            try:
                rows.append(validate_h5(tmp_path, member.name))
            finally:
                tmp_path.unlink(missing_ok=True)

    fields = [
        "sample_accession", "matrix_type", "archive_member", "status", "reason",
        "genes", "barcodes", "nnz", "total_counts", "count_dtype",
    ]
    with (part_dir / "gse175540_h5_matrix_validation_v1.tsv").open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    sample_ids = sorted({str(row.get("sample_accession", "")) for row in rows})
    failures = [row for row in rows if row.get("status") != "PASS"]
    matrix_types = {kind: sum(row.get("matrix_type") == kind for row in rows) for kind in ("raw", "filtered")}
    passed = (
        len(rows) == EXPECTED_H5_MEMBERS
        and len(sample_ids) == EXPECTED_SAMPLES
        and matrix_types == {"raw": EXPECTED_SAMPLES, "filtered": EXPECTED_SAMPLES}
        and not failures
    )
    status = {
        "schema_version": "phase51_gse175540_h5_preflight_v1",
        "status": "PASS" if passed else "FAIL",
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "archive": str(archive),
        "python": platform.python_version(),
        "h5py": h5py.__version__,
        "numpy": np.__version__,
        "h5_members": len(rows),
        "samples": len(sample_ids),
        "matrix_type_counts": matrix_types,
        "failed_members": len(failures),
        "checks": [
            "10x matrix group and required datasets exist",
            "matrix dimensions agree with features and barcodes",
            "CSC pointers and indices are internally consistent",
            "counts are finite, non-negative, and integer-valued",
        ],
    }
    (part_dir / "GSE175540_H5_PREFLIGHT_STATUS_v1.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (part_dir / ("GSE175540_H5_PREFLIGHT_PASS_v1" if passed else "GSE175540_H5_PREFLIGHT_FAIL_v1")).write_text(
        f"status={status['status']}\nh5_members={len(rows)}\nsamples={len(sample_ids)}\nfailed_members={len(failures)}\n",
        encoding="utf-8",
    )
    if not passed:
        return 1
    os.replace(part_dir, final_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
