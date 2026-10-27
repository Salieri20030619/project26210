#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import anndata as ad


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    args = parser.parse_args()
    obj = ad.read_h5ad(args.h5ad, backed="r")
    try:
        print(json.dumps({
            "shape": list(obj.shape),
            "obs_columns": list(map(str, obj.obs.columns)),
            "obsm_keys": list(map(str, obj.obsm.keys())),
            "uns_keys": list(map(str, obj.uns.keys())),
            "factor_names": list(map(str, obj.uns.get("mod", {}).get("factor_names", []))),
        }, indent=2))
    finally:
        obj.file.close()


if __name__ == "__main__":
    main()
