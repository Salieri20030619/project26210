#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--regulons", type=Path, required=True); parser.add_argument("--output", type=Path, required=True); args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    with args.regulons.open("rb") as handle: regulons = pickle.load(handle)
    rows = []
    for regulon in regulons:
        name = str(regulon.name)
        rows.append({"regulon":name, "tf":name.split("(")[0], "n_targets":len(regulon.gene2weight), "targets":";".join(sorted(map(str, regulon.gene2weight.keys())))})
    pd.DataFrame(rows).sort_values("regulon").to_csv(args.output, sep="\t", index=False, lineterminator="\n")


if __name__ == "__main__": main()
