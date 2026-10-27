#!/usr/bin/env python3

import csv
import gzip
import math
import os
import sys
from pathlib import Path

input_dir = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/"
    "external_inputs_v1/GSE178481"
)

output_file = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/"
    "phase3_input_audit_20260808_215808/"
    "gse178481_content_audit.tsv"
)

files = sorted(input_dir.glob("*.count.csv.gz"))

if len(files) != 26:
    raise RuntimeError(
        f"Expected 26 count CSV files, found {len(files)}"
    )

reference_genes = None

header = [
    "file",
    "genes",
    "cells",
    "gene_order_matches_reference",
    "invalid_values",
    "negative_values",
    "non_integer_values",
    "nonzero_values",
    "total_values",
    "nonzero_fraction",
    "total_umi",
    "max_count",
    "status",
    "reason",
]

with output_file.open("w", encoding="utf-8", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(header)

    for file_index, path in enumerate(files, start=1):
        print(
            f"[{file_index}/{len(files)}] {path.name}",
            flush=True
        )

        genes = []
        invalid_values = 0
        negative_values = 0
        non_integer_values = 0
        nonzero_values = 0
        total_values = 0
        total_umi = 0.0
        max_count = 0.0
        width_mismatch = False

        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)

            try:
                csv_header = next(reader)
            except StopIteration:
                writer.writerow([
                    path.name, 0, 0, False, 0, 0, 0, 0, 0,
                    0, 0, 0, "FAIL", "EMPTY_FILE"
                ])
                continue

            cells = len(csv_header) - 1
            expected_width = cells + 1

            for row_number, row in enumerate(reader, start=2):
                if len(row) != expected_width:
                    width_mismatch = True
                    continue

                genes.append(row[0])

                for raw_value in row[1:]:
                    total_values += 1

                    try:
                        value = float(raw_value)
                    except (TypeError, ValueError):
                        invalid_values += 1
                        continue

                    if not math.isfinite(value):
                        invalid_values += 1
                        continue

                    if value < 0:
                        negative_values += 1

                    if not value.is_integer():
                        non_integer_values += 1

                    if value != 0:
                        nonzero_values += 1

                    total_umi += value

                    if value > max_count:
                        max_count = value

        if reference_genes is None:
            reference_genes = genes
            gene_order_matches = True
        else:
            gene_order_matches = genes == reference_genes

        reasons = []

        if width_mismatch:
            reasons.append("WIDTH_MISMATCH")
        if not gene_order_matches:
            reasons.append("GENE_ORDER_MISMATCH")
        if invalid_values > 0:
            reasons.append("INVALID_VALUE")
        if negative_values > 0:
            reasons.append("NEGATIVE_VALUE")
        if non_integer_values > 0:
            reasons.append("NON_INTEGER_VALUE")

        status = "PASS" if not reasons else "FAIL"
        reason = "OK" if not reasons else ";".join(reasons)

        nonzero_fraction = (
            nonzero_values / total_values
            if total_values > 0
            else float("nan")
        )

        writer.writerow([
            path.name,
            len(genes),
            cells,
            gene_order_matches,
            invalid_values,
            negative_values,
            non_integer_values,
            nonzero_values,
            total_values,
            f"{nonzero_fraction:.10f}",
            f"{total_umi:.0f}",
            f"{max_count:.0f}",
            status,
            reason,
        ])

print(f"OUTPUT={output_file}")
