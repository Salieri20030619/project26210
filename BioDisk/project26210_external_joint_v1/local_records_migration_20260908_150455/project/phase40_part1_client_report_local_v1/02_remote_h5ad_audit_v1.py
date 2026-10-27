#!/usr/bin/env python3
"""Read-only audit of the authoritative Part 1 H5AD objects on the project server."""

import json
import os

import anndata as ad


PATHS = {
    "cell_annotated": "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/02_integration_annotation/cell_annotation_v1/part1_scvi_integrated_cell_annotated_v1.h5ad",
    "cnv_annotated": "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad",
}


def counts(obs, candidates):
    for key in candidates:
        if key in obs.columns:
            vc = obs[key].astype(str).value_counts(dropna=False)
            return key, {str(k): int(v) for k, v in vc.items()}
    return None, {}


records = []
for label, path in PATHS.items():
    obj = ad.read_h5ad(path, backed="r")
    study_key, study_counts = counts(obj.obs, ["study", "dataset", "dataset_id", "src__study"])
    patient_key, patient_counts = counts(obj.obs, ["patient_id", "patient", "Patient", "src__patient"])
    celltype_key, celltype_counts = counts(
        obj.obs,
        ["cell_type", "celltype", "cell_type_broad", "annotation", "cell_annotation"],
    )
    cnv_key, cnv_counts = counts(
        obj.obs,
        ["cnv_consensus_label", "cnv_status", "malignancy_cnv", "malignancy_label"],
    )
    eligible_key, eligible_counts = counts(
        obj.obs,
        ["eligible_for_subclone_selection_v1", "eligible_for_subclone_selection"],
    )
    records.append(
        {
            "label": label,
            "path": path,
            "size_bytes": int(os.path.getsize(path)),
            "n_obs": int(obj.n_obs),
            "n_vars": int(obj.n_vars),
            "obs_columns": [str(x) for x in obj.obs.columns],
            "study_key": study_key,
            "study_counts": study_counts,
            "patient_key": patient_key,
            "n_patients": len(patient_counts),
            "patient_counts": patient_counts,
            "celltype_key": celltype_key,
            "celltype_counts": celltype_counts,
            "cnv_key": cnv_key,
            "cnv_counts": cnv_counts,
            "eligible_key": eligible_key,
            "eligible_counts": eligible_counts,
        }
    )
    obj.file.close()

print(json.dumps({"status": "PASS", "records": records}, ensure_ascii=False, sort_keys=True))
