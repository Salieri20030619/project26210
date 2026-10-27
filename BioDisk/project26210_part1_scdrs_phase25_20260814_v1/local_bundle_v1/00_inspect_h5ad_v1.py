#!/usr/bin/env python
"""Read-only structural inspection of an H5AD used for Phase 25."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    args = parser.parse_args()
    obj = ad.read_h5ad(args.h5ad, backed="r")
    payload = {
        "path": str(args.h5ad),
        "shape": list(obj.shape),
        "x_type": type(obj.X).__name__,
        "layers": list(obj.layers.keys()),
        "raw_present": obj.raw is not None,
        "obs_columns": list(obj.obs.columns),
        "obsm": list(obj.obsm.keys()),
        "var_columns": list(obj.var.columns),
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
