#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

H5AD = Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad")
META67 = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2")
META207 = Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1/results_v1/GSE207493/gse207493_cell_metadata_scores_v1.tsv.gz")
PANEL = ["TYMP","MTAP","ADORA2B","STAT3","SPP1","ARG1","VEGFA","IL10","TIMP1","GPNMB","LGALS3","FABP5","CD9","LPL","C1QA","C1QB","C1QC","FCN1","VCAN","S100A8","S100A9","CTSS","LILRB1","TGFB1","VSIG4","MARCO","CD274","IDO1","MAT2A","MAT2B","AMD1","SRM","SMS"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    obj = ad.read_h5ad(H5AD, backed="r")
    obs = pd.Index(obj.obs_names.astype(str))
    pos = pd.Series(np.arange(len(obs), dtype=int), index=obs)
    genes = pd.Index(obj.var_names.astype(str).str.upper())
    rows = []
    for cohort in ("GSE264586", "GSE294109"):
        path = META67 / f"{cohort}_phase67_metadata_v1.tsv.gz"
        meta = pd.read_csv(path, sep="\t", compression="gzip", dtype=str)
        gsm = meta.cell_id.str.extract(r"(GSM\d+)", expand=False)
        barcode = meta.cell_id.str.rsplit("__", n=1).str[-1]
        mapped = gsm + "_" + barcode
        myeloid = meta.lineage_v1.eq("Myeloid_DC")
        eligible = meta.loc[myeloid, ["patient_key", "timepoint"]].value_counts()
        rows.append({
            "cohort": cohort, "metadata_rows": len(meta), "mapped_rows": int(mapped.isin(obs).sum()),
            "missing_rows": int((~mapped.isin(obs)).sum()), "myeloid_cells": int(myeloid.sum()),
            "myeloid_units_ge50": int((eligible >= 50).sum()), "metadata_sha256": sha256(path),
        })
    g207 = pd.read_csv(META207, sep="\t", compression="gzip", dtype=str)
    g207_ids = pd.Index(g207.cell_barcode.astype(str))
    myeloid207 = g207.cell_type.eq("Myeloid_DC")
    units207 = g207.loc[myeloid207, "sample"].value_counts()
    rows.append({
        "cohort": "GSE207493", "metadata_rows": len(g207), "mapped_rows": int(g207_ids.isin(obs).sum()),
        "missing_rows": int((~g207_ids.isin(obs)).sum()), "myeloid_cells": int(myeloid207.sum()),
        "myeloid_units_ge50": int((units207 >= 50).sum()), "metadata_sha256": sha256(META207),
    })
    matrix_sample = np.asarray(obj.X[: min(1000, obj.n_obs), : min(100, obj.n_vars)].toarray() if hasattr(obj.X[: min(1000, obj.n_obs), : min(100, obj.n_vars)], "toarray") else obj.X[: min(1000, obj.n_obs), : min(100, obj.n_vars)])
    payload = {
        "h5ad": str(H5AD), "shape": [int(obj.n_obs), int(obj.n_vars)], "h5ad_sha256": sha256(H5AD),
        "x_sample_min": float(np.nanmin(matrix_sample)), "x_sample_max": float(np.nanmax(matrix_sample)),
        "x_sample_noninteger_fraction": float(np.mean(np.abs(matrix_sample - np.rint(matrix_sample)) > 1e-8)),
        "panel_present": [g for g in PANEL if g in set(genes)], "panel_missing": [g for g in PANEL if g not in set(genes)],
        "mapping": rows,
    }
    (out / "phase71_unified_input_audit_v1.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    with (out / "phase71_unified_mapping_summary_v1.tsv").open("w", encoding="utf-8", newline="") as fh:
        fields = ["cohort","metadata_rows","mapped_rows","missing_rows","myeloid_cells","myeloid_units_ge50","metadata_sha256"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(json.dumps({"status":"PASS" if all(r["missing_rows"] == 0 for r in rows) else "FAIL", "shape":payload["shape"], "mapping":rows}))


if __name__ == "__main__":
    main()
