#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import pandas as pd


OUTPUT_DIR = Path("/mnt/BioDisk/wangmingyuan/project26210_singlecell_inventory_v1")
INPUTS = {
    "part1_cnv_annotated": Path(
        "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/"
        "cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
    ),
    "legacy_annotated": Path(
        "/mnt/BioDisk/ludingyu/project26210/object/01_annotation_adata_annotated.h5ad"
    ),
    "part4_cell2location_reference": Path(
        "/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/ref/ref_signatures.h5ad"
    ),
}
KEY_GENES = [
    "TYMP", "SPP1", "APOE", "C1QA", "C1QB", "C1QC", "TREM2", "TYROBP",
    "PECAM1", "VWF", "KDR", "EMCN", "ENG", "RGS5", "CSPG4", "ACTA2",
    "CD3D", "CD8A", "CD4", "FOXP3", "PDCD1", "LAG3", "HAVCR2", "CXCL13",
    "VEGFA", "HIF1A", "POSTN", "COL1A1", "CCL2", "CCR2",
]
INTEREST_HINTS = (
    "cell", "type", "label", "annot", "tissue", "disease", "study", "sample",
    "patient", "malign", "cnv", "tumor", "normal", "source", "major", "minor",
)


def scalar(value):
    if pd.isna(value):
        return None
    if hasattr(value, "item"):
        try:
            return value.item()
        except (ValueError, AttributeError):
            pass
    return str(value) if not isinstance(value, (str, int, float, bool)) else value


def audit_one(name: str, path: Path) -> tuple[dict, list[dict]]:
    if not path.exists():
        raise FileNotFoundError(path)
    obj = ad.read_h5ad(path, backed="r")
    obs = obj.obs
    var_names = set(map(str, obj.var_names))
    summary = {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "n_obs": int(obj.n_obs),
        "n_vars": int(obj.n_vars),
        "obs_columns": list(map(str, obs.columns)),
        "layers": list(map(str, obj.layers.keys())),
        "obsm": list(map(str, obj.obsm.keys())),
        "has_raw": obj.raw is not None,
        "key_gene_presence": {gene: gene in var_names for gene in KEY_GENES},
        "selected_value_counts": {},
    }
    rows: list[dict] = []
    for column in obs.columns:
        series = obs[column]
        nunique = int(series.nunique(dropna=False))
        column_lower = str(column).lower()
        selected = nunique <= 100 and any(hint in column_lower for hint in INTEREST_HINTS)
        if not selected:
            continue
        counts = series.value_counts(dropna=False)
        summary["selected_value_counts"][str(column)] = {
            str(scalar(level)): int(count) for level, count in counts.items()
        }
        for level, count in counts.items():
            rows.append(
                {
                    "object": name,
                    "column": str(column),
                    "level": scalar(level),
                    "n": int(count),
                    "n_unique_column": nunique,
                }
            )
    if getattr(obj, "file", None) is not None:
        obj.file.close()
    return summary, rows


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "complete",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "read_mode": "anndata backed='r'; metadata only; expression matrix not materialized",
        "objects": {},
    }
    rows: list[dict] = []
    for name, path in INPUTS.items():
        summary, object_rows = audit_one(name, path)
        result["objects"][name] = summary
        rows.extend(object_rows)
    json_path = OUTPUT_DIR / "singlecell_h5ad_metadata_inventory_v1.json"
    tsv_path = OUTPUT_DIR / "singlecell_h5ad_selected_value_counts_v1.tsv"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame(rows).to_csv(tsv_path, sep="\t", index=False)
    print(json.dumps({
        "status": "complete",
        "json": str(json_path),
        "tsv": str(tsv_path),
        "objects": {k: {"n_obs": v["n_obs"], "n_vars": v["n_vars"]} for k, v in result["objects"].items()},
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
