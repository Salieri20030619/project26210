#!/usr/bin/env python3

import csv
import gzip
import os
import re
import sys
import time
from pathlib import Path


ROOT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1"
)

INPUT_DIR = ROOT / "external_inputs_v1" / "GSE178481"

AUDIT_FILE = (
    ROOT
    / "phase3_input_audit_20260808_215808"
    / "gse178481_content_audit.tsv"
)

OUTPUT_ROOT = (
    ROOT
    / "phase4c_gse178481_seurat_v1"
    / "sparse_inputs"
)

MANIFEST_FILE = (
    ROOT
    / "phase4c_gse178481_seurat_v1"
    / "sparse_conversion_manifest.tsv"
)


def load_audit(path):
    records = {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "file",
            "genes",
            "cells",
            "nonzero_values",
            "total_umi",
            "status",
        }

        missing = required.difference(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                "Audit file missing columns: "
                + ",".join(sorted(missing))
            )

        for row in reader:
            if row["status"] != "PASS":
                raise RuntimeError(
                    f"Audit status is not PASS: {row['file']}"
                )

            records[row["file"]] = {
                "genes": int(row["genes"]),
                "cells": int(row["cells"]),
                "nonzero_values": int(row["nonzero_values"]),
                "total_umi": int(row["total_umi"]),
            }

    return records


def sample_id_from_filename(filename):
    suffix = ".count.csv.gz"

    if not filename.endswith(suffix):
        raise ValueError(f"Unexpected filename: {filename}")

    return filename[:-len(suffix)]


def parse_sample_metadata(sample_id):
    match = re.match(
        r"^GSM(?P<gsm>[0-9]+)_RCC-(?P<patient>[^-]+)-"
        r"(?P<tissue>Normal|PTumor)(?P<region>[0-9]*)$",
        sample_id,
    )

    if match is None:
        raise ValueError(
            f"Cannot parse sample metadata from: {sample_id}"
        )

    gsm = "GSM" + match.group("gsm")
    patient = match.group("patient")
    tissue_raw = match.group("tissue")
    region_raw = match.group("region")

    tissue = (
        "Normal"
        if tissue_raw == "Normal"
        else "PrimaryTumor"
    )

    region = (
        region_raw
        if region_raw
        else "1"
    )

    return {
        "gsm": gsm,
        "patient": patient,
        "tissue": tissue,
        "region": region,
    }


def write_text_lines_gzip(path, values):
    temp_path = Path(str(path) + ".tmp")

    with gzip.open(
        temp_path,
        "wt",
        encoding="utf-8",
        newline="\n",
        compresslevel=6,
    ) as handle:
        for value in values:
            handle.write(str(value))
            handle.write("\n")

    os.replace(temp_path, path)


def convert_one(input_path, expected, output_dir):
    output_dir.mkdir(parents=True, exist_ok=False)

    matrix_path = output_dir / "matrix.mtx.gz"
    genes_path = output_dir / "features.tsv.gz"
    barcodes_path = output_dir / "barcodes.tsv.gz"
    metadata_path = output_dir / "sample_metadata.tsv"
    complete_path = output_dir / ".conversion_complete"

    matrix_temp = Path(str(matrix_path) + ".tmp")
    genes_temp = Path(str(genes_path) + ".tmp")
    barcodes_temp = Path(str(barcodes_path) + ".tmp")
    metadata_temp = Path(str(metadata_path) + ".tmp")

    sample_id = sample_id_from_filename(input_path.name)
    metadata = parse_sample_metadata(sample_id)

    observed_genes = 0
    observed_cells = 0
    observed_nonzero = 0
    observed_total_umi = 0
    max_count = 0

    started = time.time()

    try:
        with gzip.open(
            input_path,
            "rt",
            encoding="utf-8-sig",
            newline="",
        ) as input_handle:
            reader = csv.reader(input_handle)

            try:
                header = next(reader)
            except StopIteration:
                raise RuntimeError(
                    f"Input file is empty: {input_path}"
                )

            if len(header) < 2:
                raise RuntimeError(
                    f"Invalid CSV header: {input_path}"
                )

            barcodes = header[1:]
            observed_cells = len(barcodes)

            if observed_cells != expected["cells"]:
                raise RuntimeError(
                    f"Cell mismatch for {input_path.name}: "
                    f"observed={observed_cells}, "
                    f"expected={expected['cells']}"
                )

            write_text_lines_gzip(barcodes_temp, barcodes)

            with gzip.open(
                genes_temp,
                "wt",
                encoding="utf-8",
                newline="\n",
                compresslevel=6,
            ) as genes_handle, gzip.open(
                matrix_temp,
                "wt",
                encoding="utf-8",
                newline="\n",
                compresslevel=6,
            ) as matrix_handle:
                matrix_handle.write(
                    "%%MatrixMarket matrix coordinate integer general\n"
                )
                matrix_handle.write(
                    "% Generated from GSE178481 dense count CSV\n"
                )
                matrix_handle.write(
                    f"{expected['genes']} "
                    f"{expected['cells']} "
                    f"{expected['nonzero_values']}\n"
                )

                expected_width = observed_cells + 1

                for gene_index, row in enumerate(reader, start=1):
                    if len(row) != expected_width:
                        raise RuntimeError(
                            f"Width mismatch in {input_path.name}, "
                            f"CSV row={gene_index + 1}, "
                            f"observed={len(row)}, "
                            f"expected={expected_width}"
                        )

                    gene = row[0]
                    genes_handle.write(gene)
                    genes_handle.write("\n")
                    observed_genes += 1

                    for cell_index, raw_value in enumerate(
                        row[1:],
                        start=1,
                    ):
                        try:
                            value = int(raw_value)
                        except ValueError:
                            numeric_value = float(raw_value)

                            if not numeric_value.is_integer():
                                raise RuntimeError(
                                    f"Non-integer value encountered in "
                                    f"{input_path.name}: {raw_value}"
                                )

                            value = int(numeric_value)

                        if value == 0:
                            continue

                        if value < 0:
                            raise RuntimeError(
                                f"Negative value encountered in "
                                f"{input_path.name}: {value}"
                            )

                        matrix_handle.write(
                            f"{gene_index} {cell_index} {value}\n"
                        )

                        observed_nonzero += 1
                        observed_total_umi += value

                        if value > max_count:
                            max_count = value

        if observed_genes != expected["genes"]:
            raise RuntimeError(
                f"Gene mismatch for {input_path.name}: "
                f"observed={observed_genes}, "
                f"expected={expected['genes']}"
            )

        if observed_nonzero != expected["nonzero_values"]:
            raise RuntimeError(
                f"Nonzero mismatch for {input_path.name}: "
                f"observed={observed_nonzero}, "
                f"expected={expected['nonzero_values']}"
            )

        if observed_total_umi != expected["total_umi"]:
            raise RuntimeError(
                f"UMI mismatch for {input_path.name}: "
                f"observed={observed_total_umi}, "
                f"expected={expected['total_umi']}"
            )

        with metadata_temp.open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.writer(
                handle,
                delimiter="\t",
                lineterminator="\n",
            )

            writer.writerow([
                "sample_id",
                "gsm",
                "patient_id",
                "tissue",
                "region",
                "source_dataset",
                "input_file",
            ])

            writer.writerow([
                sample_id,
                metadata["gsm"],
                metadata["patient"],
                metadata["tissue"],
                metadata["region"],
                "GSE178481",
                input_path.name,
            ])

        os.replace(matrix_temp, matrix_path)
        os.replace(genes_temp, genes_path)
        os.replace(barcodes_temp, barcodes_path)
        os.replace(metadata_temp, metadata_path)

        elapsed = time.time() - started

        with complete_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            handle.write("status=PASS\n")
            handle.write(f"sample_id={sample_id}\n")
            handle.write(f"genes={observed_genes}\n")
            handle.write(f"cells={observed_cells}\n")
            handle.write(f"nonzero_values={observed_nonzero}\n")
            handle.write(f"total_umi={observed_total_umi}\n")
            handle.write(f"max_count={max_count}\n")
            handle.write(f"elapsed_seconds={elapsed:.3f}\n")

        return {
            "sample_id": sample_id,
            "input_file": input_path.name,
            "genes": observed_genes,
            "cells": observed_cells,
            "nonzero_values": observed_nonzero,
            "total_umi": observed_total_umi,
            "max_count": max_count,
            "elapsed_seconds": elapsed,
            "output_dir": str(output_dir),
            "status": "PASS",
            "reason": "OK",
        }

    except Exception:
        for temp_path in [
            matrix_temp,
            genes_temp,
            barcodes_temp,
            metadata_temp,
        ]:
            if temp_path.exists():
                temp_path.unlink()

        raise


def main():
    audit = load_audit(AUDIT_FILE)
    files = sorted(INPUT_DIR.glob("*.count.csv.gz"))

    if len(files) != 26:
        raise RuntimeError(
            f"Expected 26 input files, found {len(files)}"
        )

    if set(path.name for path in files) != set(audit):
        missing_input = sorted(set(audit) - set(path.name for path in files))
        missing_audit = sorted(
            set(path.name for path in files) - set(audit)
        )

        raise RuntimeError(
            f"Input/audit mismatch. "
            f"Missing input={missing_input}; "
            f"missing audit={missing_audit}"
        )

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    MANIFEST_FILE.parent.mkdir(parents=True, exist_ok=True)

    if MANIFEST_FILE.exists():
        raise RuntimeError(
            f"Manifest already exists; refusing to overwrite: "
            f"{MANIFEST_FILE}"
        )

    manifest_temp = Path(str(MANIFEST_FILE) + ".tmp")

    fieldnames = [
        "sample_id",
        "input_file",
        "genes",
        "cells",
        "nonzero_values",
        "total_umi",
        "max_count",
        "elapsed_seconds",
        "output_dir",
        "status",
        "reason",
    ]

    total_cells = 0
    total_nonzero = 0
    total_umi = 0

    with manifest_temp.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for index, input_path in enumerate(files, start=1):
            sample_id = sample_id_from_filename(input_path.name)
            output_dir = OUTPUT_ROOT / sample_id

            print(
                f"[{index}/{len(files)}] START {input_path.name}",
                flush=True,
            )

            if output_dir.exists():
                raise RuntimeError(
                    f"Output directory already exists; "
                    f"refusing to overwrite: {output_dir}"
                )

            result = convert_one(
                input_path=input_path,
                expected=audit[input_path.name],
                output_dir=output_dir,
            )

            writer.writerow(result)
            handle.flush()

            total_cells += result["cells"]
            total_nonzero += result["nonzero_values"]
            total_umi += result["total_umi"]

            print(
                f"[{index}/{len(files)}] PASS "
                f"cells={result['cells']} "
                f"nnz={result['nonzero_values']} "
                f"umi={result['total_umi']} "
                f"seconds={result['elapsed_seconds']:.1f}",
                flush=True,
            )

    os.replace(manifest_temp, MANIFEST_FILE)

    print("===== COMPLETION =====")
    print("status=PASS")
    print(f"files={len(files)}")
    print(f"cells={total_cells}")
    print(f"nonzero_values={total_nonzero}")
    print(f"total_umi={total_umi}")
    print(f"manifest={MANIFEST_FILE}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(
            f"FATAL: {type(exc).__name__}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)
