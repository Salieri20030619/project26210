#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
OUT = Path("/mnt/BioDisk/wangmingyuan/project26210_singlecell_inventory_v1")
REFERENCE = ROOT / "results/part4_spatial/ref/ref_signatures.h5ad"
MAPPING = ROOT / "results/part1_rebuild_v1/00_input_audit/geo_sample_provenance_v1.tsv"
KEY_GENES = [
    "TYMP", "SPP1", "APOE", "C1QA", "C1QB", "C1QC", "TREM2", "TYROBP",
    "PECAM1", "VWF", "KDR", "EMCN", "ENG", "RGS5", "CSPG4", "ACTA2",
    "CD3D", "CD8A", "CD4", "FOXP3", "PDCD1", "LAG3", "HAVCR2", "CXCL13",
    "VEGFA", "HIF1A", "POSTN", "COL1A1", "CCL2", "CCR2",
]


def coverage_table(obs: pd.DataFrame, group_col: str) -> pd.DataFrame:
    patient_group = (
        obs.groupby(["curated_patient_id", group_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    sample_group = (
        obs.groupby(["sample_key", group_col], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    rows = []
    for group, frame in obs.groupby(group_col, observed=True):
        patient_counts = patient_group.loc[patient_group[group_col].eq(group), "n_cells"]
        sample_counts = sample_group.loc[sample_group[group_col].eq(group), "n_cells"]
        rows.append({
            group_col: str(group),
            "n_cells": int(len(frame)),
            "n_studies": int(frame["study"].nunique()),
            "n_patients": int(frame["curated_patient_id"].nunique()),
            "n_patients_ge20_cells": int((patient_counts >= 20).sum()),
            "n_samples": int(frame["sample_key"].nunique()),
            "n_samples_ge20_cells": int((sample_counts >= 20).sum()),
            "n_tumor_cells": int(frame["Disease"].eq("ccRCC").sum()),
            "n_normal_cells": int(frame["Disease"].eq("Normal").sum()),
        })
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def expression_summary(obj: ad.AnnData, obs: pd.DataFrame, group_col: str) -> pd.DataFrame:
    present_genes = [gene for gene in KEY_GENES if gene in obj.var_names]
    indices = [int(obj.var_names.get_loc(gene)) for gene in present_genes]
    matrix = obj.layers["counts"][:, indices]
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix)
    rows = []
    groups = obs[group_col].astype(str).to_numpy()
    patients = obs["curated_patient_id"].astype(str).to_numpy()
    for group in sorted(pd.unique(groups)):
        mask = groups == group
        group_patients = patients[mask]
        for j, gene in enumerate(present_genes):
            values = matrix[mask, j]
            detected = values > 0
            patient_detection = pd.DataFrame({
                "patient": group_patients,
                "detected": detected,
            }).groupby("patient", observed=True)["detected"].mean()
            rows.append({
                group_col: group,
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
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(MAPPING, sep="\t")
    mapping = mapping[["study", "local_sample_id", "curated_patient_id", "source_tissue"]].copy()
    mapping["sample_key"] = mapping["study"].astype(str) + "||" + mapping["local_sample_id"].astype(str)
    if mapping["sample_key"].duplicated().any():
        raise RuntimeError("Non-unique study/local_sample_id mapping")

    obj = ad.read_h5ad(REFERENCE, backed="r")
    obs = obj.obs.copy()
    obs["sample_key"] = obs["study"].astype(str) + "||" + obs["sample_id"].astype(str)
    obs = obs.merge(
        mapping[["sample_key", "curated_patient_id", "source_tissue"]],
        on="sample_key",
        how="left",
        validate="many_to_one",
    )
    if obs["curated_patient_id"].isna().any():
        missing = sorted(obs.loc[obs["curated_patient_id"].isna(), "sample_key"].unique())
        raise RuntimeError(f"Reference samples missing patient mapping: {missing}")

    outputs = {}
    for group_col in ["c2l_celltype", "macro_subtype"]:
        coverage = coverage_table(obs, group_col)
        coverage_path = OUT / f"{group_col}_patient_sample_coverage_v1.tsv"
        coverage.to_csv(coverage_path, sep="\t", index=False)
        expression = expression_summary(obj, obs, group_col)
        expression_path = OUT / f"{group_col}_key_gene_detection_v1.tsv"
        expression.to_csv(expression_path, sep="\t", index=False)
        study_path = OUT / f"{group_col}_by_study_v1.tsv"
        pd.crosstab(obs[group_col], obs["study"]).to_csv(study_path, sep="\t")
        outputs[group_col] = {
            "coverage": str(coverage_path),
            "expression": str(expression_path),
            "by_study": str(study_path),
        }

    metadata = {
        "status": "complete",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "reference": str(REFERENCE),
        "mapping": str(MAPPING),
        "n_cells": int(len(obs)),
        "n_patients": int(obs["curated_patient_id"].nunique()),
        "n_samples": int(obs["sample_key"].nunique()),
        "n_studies": int(obs["study"].nunique()),
        "unmapped_cells": int(obs["curated_patient_id"].isna().sum()),
        "key_genes": KEY_GENES,
        "outputs": outputs,
    }
    metadata_path = OUT / "joint_reference_coverage_metadata_v1.json"
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    if getattr(obj, "file", None) is not None:
        obj.file.close()
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
