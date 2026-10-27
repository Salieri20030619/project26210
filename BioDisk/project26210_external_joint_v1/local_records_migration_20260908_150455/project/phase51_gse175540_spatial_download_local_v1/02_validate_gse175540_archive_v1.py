#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_ARCHIVE_BYTES = 1_122_877_440
EXPECTED_SAMPLES = 24
EXPECTED_FFPE = 12
EXPECTED_FROZEN = 12
REQUIRED_SUFFIXES = {
    "TLS_annotation.csv.gz": 23,
    "aligned_fiducials.jpg.gz": 24,
    "detected_tissue_image.jpg.gz": 24,
    "filtered_feature_bc_matrix.h5": 24,
    "raw_feature_bc_matrix.h5": 24,
    "scalefactors_json.json.gz": 24,
    "tissue_hires_image.png.gz": 24,
    "tissue_lowres_image.png.gz": 24,
    "tissue_positions_list.csv.gz": 24,
}
EXPECTED_TLS_MISSING_SAMPLE = "GSM5924034"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, value: str) -> None:
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(value, encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python 02_validate_gse175540_archive_v1.py <download_dir.part>")
    root = Path(sys.argv[1]).resolve()
    downloads = root / "downloads"
    archive = downloads / "GSE175540_RAW.tar"
    filelist = downloads / "filelist.txt"
    metadata = downloads / "GSE175540_family.xml.tgz"

    if archive.stat().st_size != EXPECTED_ARCHIVE_BYTES:
        raise RuntimeError("Official archive byte count does not match")
    if filelist.stat().st_size == 0 or metadata.stat().st_size == 0:
        raise RuntimeError("Downloaded metadata file is empty")

    official_sizes: dict[str, int] = {}
    with filelist.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["#Archive/File"] == "File":
                official_sizes[row["Name"]] = int(row["Size"])

    with tarfile.open(archive, mode="r:") as handle:
        members = [member for member in handle.getmembers() if member.isfile()]
    observed_sizes = {Path(member.name).name: member.size for member in members}
    if len(observed_sizes) != len(members):
        raise RuntimeError("Archive contains duplicate basenames")
    if observed_sizes != official_sizes:
        missing = sorted(set(official_sizes) - set(observed_sizes))
        extra = sorted(set(observed_sizes) - set(official_sizes))
        size_bad = sorted(name for name in set(observed_sizes) & set(official_sizes) if observed_sizes[name] != official_sizes[name])
        raise RuntimeError(f"Archive/filelist mismatch missing={missing} extra={extra} size_bad={size_bad}")

    sample_ids = sorted({match.group(0) for name in observed_sizes if (match := re.match(r"GSM\d+", name))})
    if len(sample_ids) != EXPECTED_SAMPLES:
        raise RuntimeError(f"Expected 24 GSM samples, observed {len(sample_ids)}")

    suffix_counts = Counter()
    for name in observed_sizes:
        for suffix in REQUIRED_SUFFIXES:
            if name.endswith(suffix):
                suffix_counts[suffix] += 1
    if dict(suffix_counts) != REQUIRED_SUFFIXES:
        raise RuntimeError(f"Required file-type counts differ: {dict(suffix_counts)}")

    tls_samples = {
        sample
        for sample in sample_ids
        if any(name.startswith(sample + "_") and name.endswith("TLS_annotation.csv.gz") for name in observed_sizes)
    }
    tls_missing = sorted(set(sample_ids) - tls_samples)
    if tls_missing != [EXPECTED_TLS_MISSING_SAMPLE]:
        raise RuntimeError(f"Unexpected TLS annotation coverage: missing={tls_missing}")

    ffpe = {sample for sample in sample_ids if any(name.startswith(sample + "_ffpe_") for name in observed_sizes)}
    frozen = {sample for sample in sample_ids if any(name.startswith(sample + "_frozen_") for name in observed_sizes)}
    if len(ffpe) != EXPECTED_FFPE or len(frozen) != EXPECTED_FROZEN or ffpe & frozen:
        raise RuntimeError(f"Preparation split invalid: FFPE={len(ffpe)} frozen={len(frozen)}")

    archive_sha = sha256(archive)
    sample_table = root / "gse175540_sample_inventory_v1.tsv"
    rows = ["sample_id\tpreparation\tmember_count\tofficial_member_bytes"]
    for sample in sample_ids:
        names = [name for name in observed_sizes if name.startswith(sample + "_")]
        preparation = "FFPE" if sample in ffpe else "frozen"
        rows.append(f"{sample}\t{preparation}\t{len(names)}\t{sum(observed_sizes[name] for name in names)}")
    atomic_text(sample_table, "\n".join(rows) + "\n")

    status = {
        "schema_version": "phase51_gse175540_archive_validation_v1",
        "status": "PASS",
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "official_accession": "GSE175540",
        "archive_bytes": archive.stat().st_size,
        "archive_sha256": archive_sha,
        "archive_file_members": len(members),
        "official_filelist_rows": len(official_sizes),
        "samples": len(sample_ids),
        "ffpe_samples": len(ffpe),
        "frozen_samples": len(frozen),
        "required_suffix_counts": REQUIRED_SUFFIXES,
        "tls_annotation_samples": len(tls_samples),
        "tls_annotation_missing_sample": EXPECTED_TLS_MISSING_SAMPLE,
        "raw_reads_public": False,
        "processed_spatial_package_public": True,
    }
    atomic_text(root / "GSE175540_DOWNLOAD_VALIDATION_STATUS_v1.json", json.dumps(status, indent=2) + "\n")
    atomic_text(
        root / "GSE175540_DOWNLOAD_PASS_v1",
        "\n".join(
            [
                "status=PASS",
                f"archive_bytes={status['archive_bytes']}",
                f"archive_sha256={archive_sha}",
                f"archive_file_members={status['archive_file_members']}",
                f"samples={status['samples']}",
                f"ffpe_samples={status['ffpe_samples']}",
                f"frozen_samples={status['frozen_samples']}",
            ]
        )
        + "\n",
    )


if __name__ == "__main__":
    main()
