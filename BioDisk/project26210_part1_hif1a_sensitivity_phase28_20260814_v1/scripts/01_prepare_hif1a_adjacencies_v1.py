#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genie3", required=True, type=Path)
    parser.add_argument("--grnboost2", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    part = args.output.with_name(args.output.name + ".part")
    if part.exists():
        raise FileExistsError(part)
    part.mkdir(parents=True)

    rows = []
    inputs = []
    for method, path in (("GENIE3", args.genie3), ("GRNBoost2", args.grnboost2)):
        if not path.is_file():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path, sep="\t")
        if list(frame.columns) != ["TF", "target", "importance"]:
            raise ValueError(f"unexpected adjacency schema: {path}")
        if frame[["TF", "target"]].isna().any().any() or frame["importance"].isna().any():
            raise ValueError(f"missing adjacency values: {path}")
        hif = frame.loc[frame["TF"].astype(str).eq("HIF1A")].copy()
        if len(hif) < 20 or hif["target"].duplicated().any() or (hif["importance"] <= 0).any():
            raise ValueError(f"HIF1A edge safety gate failed for {method}")
        name = method.lower() + "_hif1a_adjacencies_v1.tsv"
        hif.to_csv(part / name, sep="\t", index=False)
        rows.append({
            "method": method,
            "n_all_edges": len(frame),
            "n_hif1a_edges": len(hif),
            "n_unique_targets": hif["target"].nunique(),
            "max_importance": float(hif["importance"].max()),
            "min_importance": float(hif["importance"].min()),
            "output": name,
        })
        inputs.append({"method": method, "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)})
    pd.DataFrame(rows).to_csv(part / "hif1a_adjacency_summary_v1.tsv", sep="\t", index=False)
    pd.DataFrame(inputs).to_csv(part / "INPUT_MANIFEST_v1.tsv", sep="\t", index=False)
    (part / "status.json").write_text(json.dumps({"status": "PASS", "methods": 2, "tf": "HIF1A"}, indent=2) + "\n", encoding="utf-8")
    os.replace(part, args.output)


if __name__ == "__main__":
    main()

