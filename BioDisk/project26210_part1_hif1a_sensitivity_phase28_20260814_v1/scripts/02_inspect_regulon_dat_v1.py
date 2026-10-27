#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dat", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--present-marker", required=True, type=Path)
    args = parser.parse_args()
    for path in (args.summary, args.targets, args.present_marker):
        if path.exists():
            raise FileExistsError(path)
    with args.dat.open("rb") as handle:
        regulons = pickle.load(handle)
    rows = []
    targets = []
    for regulon in regulons:
        name = str(regulon.name)
        gene2weight = dict(regulon.gene2weight)
        rows.append({"regulon": name, "n_targets": len(gene2weight), "context": ";".join(sorted(map(str, regulon.context)))})
        for gene, weight in sorted(gene2weight.items()):
            targets.append({"regulon": name, "target": str(gene), "weight": float(weight)})
    pd.DataFrame(rows, columns=["regulon", "n_targets", "context"]).to_csv(args.summary, sep="\t", index=False)
    pd.DataFrame(targets, columns=["regulon", "target", "weight"]).to_csv(args.targets, sep="\t", index=False)
    positive = [x for x in rows if x["regulon"] == "HIF1A(+)"]
    if positive:
        args.present_marker.write_text(json.dumps({"regulon": "HIF1A(+)", "n_targets": positive[0]["n_targets"]}) + "\n", encoding="utf-8")
    print(json.dumps({"n_regulons": len(rows), "hif1a_positive_present": bool(positive)}))


if __name__ == "__main__":
    main()

