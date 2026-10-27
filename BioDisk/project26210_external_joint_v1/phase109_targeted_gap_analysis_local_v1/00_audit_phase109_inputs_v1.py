from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", type=Path, required=True)
parser.add_argument("--assignments", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)

assign = pd.read_csv(args.assignments, sep="\t")
adata = ad.read_h5ad(args.h5ad, backed="r")
obs_names = pd.Index(adata.obs_names.astype(str))
assignment_ids = pd.Index(assign["exact_cell_id"].astype(str))
indexer = obs_names.get_indexer(assignment_ids)

candidate_tokens = [
    "patient",
    "study",
    "sample",
    "dataset",
    "cell",
    "annot",
    "malig",
    "cnv",
    "f1",
    "tymp",
    "mta",
    "polyamine",
    "metab",
]
column_rows = []
for column in adata.obs.columns:
    values = adata.obs[column]
    column_rows.append(
        {
            "column": str(column),
            "dtype": str(values.dtype),
            "nonmissing": int(values.notna().sum()),
            "unique_nonmissing": int(values.nunique(dropna=True)),
            "candidate": any(token in str(column).lower() for token in candidate_tokens),
            "examples": ";".join(map(str, values.dropna().astype(str).unique()[:8])),
        }
    )
pd.DataFrame(column_rows).to_csv(args.output / "h5ad_obs_column_audit_v1.tsv", sep="\t", index=False)

state_counts = (
    assign.groupby(["study_v1", "state_v2"], dropna=False)
    .agg(cells=("exact_cell_id", "size"), patients=("patient", "nunique"))
    .reset_index()
)
state_counts.to_csv(args.output / "assignment_state_counts_v1.tsv", sep="\t", index=False)

payload = {
    "status": "PASS_INPUT_AUDIT" if (indexer >= 0).all() else "FAIL_INPUT_AUDIT",
    "h5ad_path": str(args.h5ad),
    "h5ad_sha256": sha256(args.h5ad),
    "h5ad_shape": [int(adata.n_obs), int(adata.n_vars)],
    "h5ad_obs_columns": list(map(str, adata.obs.columns)),
    "h5ad_layers": list(map(str, adata.layers.keys())),
    "h5ad_obsm": list(map(str, adata.obsm.keys())),
    "assignments_path": str(args.assignments),
    "assignments_sha256": sha256(args.assignments),
    "assignment_rows": int(len(assign)),
    "assignment_unique_ids": int(assign.exact_cell_id.nunique()),
    "assignment_patients": int(assign.patient.nunique()),
    "assignment_studies": sorted(assign.study_v1.unique().tolist()),
    "assignment_states": sorted(assign.state_v2.unique().tolist()),
    "assignment_ids_missing_from_h5ad": int((indexer < 0).sum()),
    "candidate_obs_columns": [row["column"] for row in column_rows if row["candidate"]],
}
(args.output / "PHASE109_INPUT_AUDIT_STATUS_v1.json").write_text(
    json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(payload, ensure_ascii=False, indent=2))
raise SystemExit(0 if payload["status"] == "PASS_INPUT_AUDIT" else 1)
