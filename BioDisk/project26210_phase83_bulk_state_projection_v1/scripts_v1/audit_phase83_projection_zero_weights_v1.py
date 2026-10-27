#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--weights", required=True, type=Path)
p.add_argument("--diagnostics", required=True, type=Path)
p.add_argument("--output", required=True, type=Path)
a = p.parse_args()
if a.output.exists():
    raise FileExistsError(a.output)
a.output.mkdir(parents=True)
w = pd.read_csv(a.weights, sep="\t")
d = pd.read_csv(a.diagnostics, sep="\t")
summary = (
    w.groupby(["dataset", "state"])["weight"]
    .agg(n="size", zero_count=lambda x: int((x <= 1e-12).sum()),
         minimum="min", q25=lambda x: x.quantile(.25), median="median",
         q75=lambda x: x.quantile(.75), maximum="max")
    .reset_index()
)
summary["zero_fraction"] = summary["zero_count"] / summary["n"]
summary.to_csv(a.output / "phase83_state_weight_distribution_v1.tsv", sep="\t", index=False)
diag = d.groupby("dataset").agg(
    n=("sample_id", "size"),
    median_relative_reconstruction_error=("relative_reconstruction_error", "median"),
    q25_relative_reconstruction_error=("relative_reconstruction_error", lambda x: x.quantile(.25)),
    q75_relative_reconstruction_error=("relative_reconstruction_error", lambda x: x.quantile(.75)),
    median_active_states=("active_states", "median"),
).reset_index()
diag.to_csv(a.output / "phase83_projection_diagnostic_summary_v1.tsv", sep="\t", index=False)
print(summary[summary["state"].str.contains("Macrophage")].to_string(index=False))
