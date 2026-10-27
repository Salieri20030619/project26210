#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


INPUT = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/"
    "cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
)
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_singlecell_inventory_v1")
KEY_GENES = [
    "TYMP", "SPP1", "APOE", "C1QA", "C1QB", "C1QC", "TREM2", "TYROBP",
    "PECAM1", "VWF", "KDR", "EMCN", "ENG", "RGS5", "CSPG4", "ACTA2",
    "CD3D", "CD8A", "CD4", "FOXP3", "PDCD1", "LAG3", "HAVCR2", "CXCL13",
    "VEGFA", "HIF1A", "POSTN", "COL1A1", "CCL2", "CCR2",
]
FOCUS_GROUPS = {
    "Macrophage", "Monocyte", "Endothelial", "T cell", "Pericyte/SMC",
    "Conservative malignant",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    obj = ad.read_h5ad(INPUT, backed="r")
    obs = obj.obs.copy()
    obs["joint_group_v1"] = obs["cell_type_broad_v1"].astype(str)
    malignant = obs["malignant_conservative_v1"].astype(bool)
    obs.loc[malignant, "joint_group_v1"] = "Conservative malignant"
    obs["sample_key"] = obs["study"].astype(str) + "||" + obs["sample_id"].astype(str)

    patient_group = (
        obs.groupby(["curated_patient_id", "joint_group_v1"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    rows = []
    for group, frame in obs.groupby("joint_group_v1", observed=True):
        pcounts = patient_group.loc[patient_group["joint_group_v1"].eq(group), "n_cells"]
        rows.append({
            "joint_group_v1": group,
            "n_cells": int(len(frame)),
            "n_studies": int(frame["study"].nunique()),
            "n_patients": int(frame["curated_patient_id"].nunique()),
            "n_patients_ge20_cells": int((pcounts >= 20).sum()),
            "n_samples": int(frame["sample_key"].nunique()),
            "n_tumor_cells": int(frame["source_tissue"].eq("Tumor").sum()),
            "n_normal_cells": int(frame["source_tissue"].eq("Normal").sum()),
        })
    coverage = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    coverage_path = OUT / "current_joint_cell_pool_coverage_v1.tsv"
    coverage.to_csv(coverage_path, sep="\t", index=False)

    present = [gene for gene in KEY_GENES if gene in obj.var_names]
    indices = [int(obj.var_names.get_loc(gene)) for gene in present]
    matrix = obj.layers["counts"][:, indices]
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)

    expression_rows = []
    study_rows = []
    groups = obs["joint_group_v1"].astype(str).to_numpy()
    studies = obs["study"].astype(str).to_numpy()
    patients = obs["curated_patient_id"].astype(str).to_numpy()
    for group in sorted(FOCUS_GROUPS):
        mask = groups == group
        if not mask.any():
            continue
        for j, gene in enumerate(present):
            values = matrix[mask, j]
            detected = values > 0
            patient_detection = pd.DataFrame({
                "patient": patients[mask], "detected": detected,
            }).groupby("patient", observed=True)["detected"].mean()
            expression_rows.append({
                "joint_group_v1": group,
                "gene": gene,
                "n_cells": int(mask.sum()),
                "detected_cells": int(detected.sum()),
                "detection_fraction": float(detected.mean()),
                "mean_raw_count": float(values.mean()),
                "mean_raw_count_among_detected": float(values[detected].mean()) if detected.any() else 0.0,
                "n_patients": int(patient_detection.size),
                "patients_with_any_detection": int((patient_detection > 0).sum()),
                "median_patient_detection_fraction": float(patient_detection.median()),
            })
            for study in sorted(pd.unique(studies[mask])):
                submask = mask & (studies == study)
                subvalues = matrix[submask, j]
                subdetected = subvalues > 0
                study_rows.append({
                    "joint_group_v1": group,
                    "study": study,
                    "gene": gene,
                    "n_cells": int(submask.sum()),
                    "detected_cells": int(subdetected.sum()),
                    "detection_fraction": float(subdetected.mean()),
                    "mean_raw_count": float(subvalues.mean()),
                })
    expression_path = OUT / "current_joint_cell_pool_key_gene_detection_v1.tsv"
    study_path = OUT / "current_joint_cell_pool_key_gene_by_study_v1.tsv"
    pd.DataFrame(expression_rows).to_csv(expression_path, sep="\t", index=False)
    pd.DataFrame(study_rows).to_csv(study_path, sep="\t", index=False)

    metadata = {
        "status": "complete",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input": str(INPUT),
        "n_cells": int(obj.n_obs),
        "n_genes": int(obj.n_vars),
        "outputs": {
            "coverage": str(coverage_path),
            "expression": str(expression_path),
            "expression_by_study": str(study_path),
        },
    }
    metadata_path = OUT / "current_joint_cell_pool_metadata_v1.json"
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    if getattr(obj, "file", None) is not None:
        obj.file.close()
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
