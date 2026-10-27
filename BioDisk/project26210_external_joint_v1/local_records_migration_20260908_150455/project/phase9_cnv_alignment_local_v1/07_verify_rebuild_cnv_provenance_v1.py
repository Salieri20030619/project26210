#!/usr/bin/env python3
"""Read-only provenance check proving the Phase 9 CNV source belongs to Part 1 rebuild."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def serializable(value):
    if isinstance(value, dict):
        return {str(k): serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value) if not isinstance(value, (str, int, float, bool, type(None))) else value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnv-h5ad", required=True)
    parser.add_argument("--source-h5ad", required=True)
    parser.add_argument("--build-metadata", required=True)
    args = parser.parse_args()

    cnv_path = Path(args.cnv_h5ad).resolve()
    source_path = Path(args.source_h5ad).resolve()
    metadata_path = Path(args.build_metadata).resolve()
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    checks = {
        "cnv_path_is_rebuild": "/results/part1_rebuild_v1/" in str(cnv_path),
        "source_path_is_rebuild": "/results/part1_rebuild_v1/" in str(source_path),
        "metadata_source_path_exact": Path(metadata["source_h5ad"]).resolve() == source_path,
        "metadata_output_path_exact": Path(metadata["output_h5ad"]).resolve() == cnv_path,
        "source_sha256_exact": sha256(source_path) == metadata["source_h5ad_sha256"],
        "cnv_sha256_exact": sha256(cnv_path) == metadata["output_h5ad_sha256"],
    }

    source = ad.read_h5ad(source_path, backed="r")
    cnv = ad.read_h5ad(cnv_path, backed="r")
    checks.update({
        "shape_exact": source.shape == cnv.shape == tuple(metadata["shape"]),
        "obs_names_exact": source.obs_names.equals(cnv.obs_names),
        "var_names_exact": source.var_names.equals(cnv.var_names),
        "source_has_part1_rebuild_uns": "part1_rebuild" in source.uns,
        "source_has_part1_rebuild_v1_uns": "part1_rebuild_v1" in source.uns,
        "cnv_has_part1_rebuild_uns": "part1_rebuild" in cnv.uns,
        "cnv_has_part1_rebuild_v1_uns": "part1_rebuild_v1" in cnv.uns,
        "cnv_has_malignancy_cnv_uns": "part1_malignancy_cnv_v1" in cnv.uns,
    })
    status = "PASS" if all(checks.values()) else "FAIL"
    report = {
        "status": status,
        "checks": checks,
        "source_h5ad": str(source_path),
        "source_sha256": metadata["source_h5ad_sha256"],
        "cnv_h5ad": str(cnv_path),
        "cnv_sha256": metadata["output_h5ad_sha256"],
        "source_uns_part1_rebuild": serializable(source.uns.get("part1_rebuild")),
        "source_uns_part1_rebuild_v1": serializable(source.uns.get("part1_rebuild_v1")),
        "cnv_uns_part1_malignancy_cnv_v1": serializable(cnv.uns.get("part1_malignancy_cnv_v1")),
    }
    print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
