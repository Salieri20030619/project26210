#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
from pathlib import Path

BASE = Path("/mnt/BioDisk/wangmingyuan/project26210_immunotherapy_scrna_download_precheck_v1")
OUT = BASE / "precheck" / "tenx_sample_precheck_v1.tsv"
SUMMARY = BASE / "precheck" / "tenx_sample_precheck_summary_v1.tsv"


def count_lines(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for _ in handle)


def matrix_dimensions(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("%"):
                rows, cols, nnz = line.split()[:3]
                return int(rows), int(cols), int(nnz)
    raise ValueError(f"No MatrixMarket dimension line found: {path}")


def sample_rows(dataset: str, folder: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for matrix in sorted(folder.rglob("*.matrix.mtx.gz")):
        stem = matrix.name.removesuffix(".matrix.mtx.gz")
        barcode = folder / f"{stem}.barcodes.tsv.gz"
        feature = folder / f"{stem}.features.tsv.gz"
        complete = barcode.exists() and feature.exists()
        gsm = stem.split("_", 1)[0]
        if complete:
            barcode_count = count_lines(barcode)
            feature_count = count_lines(feature)
            matrix_rows, matrix_cols, nnz = matrix_dimensions(matrix)
            dimension_match = matrix_rows == feature_count and matrix_cols == barcode_count
        else:
            barcode_count = feature_count = matrix_rows = matrix_cols = nnz = ""
            dimension_match = False
        rows.append(
            {
                "dataset": dataset,
                "sample_key": stem,
                "gsm": gsm,
                "files_complete": complete,
                "barcodes": barcode_count,
                "features": feature_count,
                "matrix_rows": matrix_rows,
                "matrix_cols": matrix_cols,
                "matrix_nnz": nnz,
                "dimension_match": dimension_match,
            }
        )
    return rows


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    all_rows = []
    all_rows.extend(sample_rows("GSE264586", BASE / "downloads" / "GSE264586" / "extracted"))
    all_rows.extend(sample_rows("GSE294109", BASE / "downloads" / "GSE294109" / "files"))

    fields = [
        "dataset",
        "sample_key",
        "gsm",
        "files_complete",
        "barcodes",
        "features",
        "matrix_rows",
        "matrix_cols",
        "matrix_nnz",
        "dimension_match",
    ]
    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    with SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["dataset", "libraries", "total_barcodes", "complete_libraries", "dimension_match_libraries"])
        for dataset in ("GSE264586", "GSE294109"):
            subset = [row for row in all_rows if row["dataset"] == dataset]
            writer.writerow(
                [
                    dataset,
                    len(subset),
                    sum(int(row["barcodes"]) for row in subset if row["barcodes"] != ""),
                    sum(row["files_complete"] is True for row in subset),
                    sum(row["dimension_match"] is True for row in subset),
                ]
            )


if __name__ == "__main__":
    main()
