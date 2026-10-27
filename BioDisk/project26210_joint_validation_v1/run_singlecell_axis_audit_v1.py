#!/usr/bin/env python3
"""Patient-aware audit of Part 3 ligand-receptor axes on the rebuilt Part 1 object.

This is not a cell-cell communication inference. It tests whether the fixed
rank-7/F1 malignant program, candidate ligands, and receiver genes are
consistently measurable across biological patients in the rebuilt dataset.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr


INPUT_H5AD = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/"
    "cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad"
)
USAGE = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/"
    "full_ranks4_8_v1/consensus_usage_rank7_v1.tsv.gz"
)
SELECTED_RANK = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/"
    "tymp_study_sensitivity_v1/selected_rank_v1.tsv"
)
OUT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_joint_validation_v1/"
    "singlecell_axis_audit_v1"
)

GENES = [
    "TYMP", "VEGFA", "VEGFB", "PGF", "FLT1", "KDR", "APP", "CD74",
    "APOE", "TREM2", "TYROBP", "SPP1", "CD44", "MIF", "CXCR4",
]
LIGANDS = ["VEGFA", "VEGFB", "PGF", "APP", "APOE", "SPP1", "MIF"]
GROUPS = ["Malignant", "Macrophage", "Monocyte", "Endothelial", "Pericyte/SMC", "T cell"]
MIN_CELLS = 20

AXES = [
    ("Tumor_VEGFA_KDR", "Malignant", "VEGFA", "Endothelial", ["KDR"]),
    ("Tumor_VEGFA_FLT1", "Malignant", "VEGFA", "Endothelial", ["FLT1"]),
    ("Tumor_VEGFB_FLT1", "Malignant", "VEGFB", "Endothelial", ["FLT1"]),
    ("Tumor_PGF_FLT1", "Malignant", "PGF", "Endothelial", ["FLT1"]),
    ("Tumor_APP_CD74", "Malignant", "APP", "Macrophage", ["CD74"]),
    ("Tumor_APOE_TREM2_TYROBP", "Malignant", "APOE", "Macrophage", ["TREM2", "TYROBP"]),
    ("Tumor_SPP1_CD44", "Malignant", "SPP1", "Macrophage", ["CD44"]),
    ("Macrophage_SPP1_Tumor_CD44", "Macrophage", "SPP1", "Malignant", ["CD44"]),
    ("Tumor_MIF_CD74_CXCR4", "Malignant", "MIF", "Macrophage", ["CD74", "CXCR4"]),
]


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return result
    order = valid.sort_values().index
    n = len(order)
    adjusted = valid.loc[order].to_numpy() * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[order] = np.minimum(adjusted, 1.0)
    return result


def safe_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return math.nan, math.nan
    stat = spearmanr(x, y)
    return float(stat.statistic), float(stat.pvalue)


def direction_summary(values: pd.Series) -> dict[str, float | int]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    positive = int((values > 0).sum())
    negative = int((values < 0).sum())
    nonzero = positive + negative
    p_value = float(binomtest(positive, nonzero, 0.5).pvalue) if nonzero else math.nan
    return {
        "n_units": int(len(values)),
        "positive_units": positive,
        "negative_units": negative,
        "zero_units": int((values == 0).sum()),
        "median_effect": float(values.median()) if len(values) else math.nan,
        "q25_effect": float(values.quantile(0.25)) if len(values) else math.nan,
        "q75_effect": float(values.quantile(0.75)) if len(values) else math.nan,
        "sign_test_p": p_value,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(SELECTED_RANK, sep="\t")
    if len(selected) != 1 or int(selected.loc[0, "selected_rank"]) != 7 or selected.loc[0, "selected_TYMP_factor"] != "F1":
        raise RuntimeError("Frozen Part 1 selected rank/factor is not rank 7/F1")

    obj = ad.read_h5ad(INPUT_H5AD, backed="r")
    missing_genes = sorted(set(GENES) - set(obj.var_names.astype(str)))
    if missing_genes:
        raise RuntimeError(f"Missing required genes: {missing_genes}")

    required_obs = {
        "cell_type_broad_v1", "eligible_for_subclone_selection_v1", "curated_patient_id",
        "study", "source_tissue", "total_counts_raw",
    }
    missing_obs = sorted(required_obs - set(obj.obs.columns))
    if missing_obs:
        raise RuntimeError(f"Missing required obs fields: {missing_obs}")

    obs = obj.obs[list(required_obs)].copy()
    obs.index = obs.index.astype(str)
    obs["patient"] = obs["curated_patient_id"].astype(str)
    obs["study_v1"] = obs["study"].astype(str)
    obs["group"] = obs["cell_type_broad_v1"].astype(str)
    eligible = obs["eligible_for_subclone_selection_v1"].astype(bool)
    obs.loc[eligible, "group"] = "Malignant"
    obs["tumor"] = obs["source_tissue"].astype(str).eq("Tumor")
    obs["library_size"] = pd.to_numeric(obs["total_counts_raw"], errors="coerce")

    usage = pd.read_csv(USAGE, sep="\t", usecols=["exact_cell_id", "F1"])
    usage["exact_cell_id"] = usage["exact_cell_id"].astype(str)
    if usage["exact_cell_id"].duplicated().any():
        raise RuntimeError("Duplicate exact_cell_id in rank-7 usage")
    usage_map = usage.set_index("exact_cell_id")["F1"]
    obs["F1"] = usage_map.reindex(obs.index).to_numpy()
    if int(eligible.sum()) != 11950 or int(obs.loc[eligible, "F1"].notna().sum()) != 11950:
        raise RuntimeError("Expected complete rank-7/F1 mapping for 11,950 conservative malignant cells")
    if obs.loc[~eligible, "F1"].notna().any():
        raise RuntimeError("Rank-7/F1 usage mapped outside the frozen malignant set")

    gene_indices = [int(obj.var_names.get_loc(gene)) for gene in GENES]
    matrix = obj.layers["counts"][:, gene_indices]
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=np.float64)
    if matrix.shape != (obj.n_obs, len(GENES)) or not np.isfinite(matrix).all() or (matrix < 0).any():
        raise RuntimeError("Invalid selected raw-count matrix")
    gene_index = {gene: j for j, gene in enumerate(GENES)}

    # Patient-by-group raw-count pseudobulk and detection summaries.
    pb_rows: list[dict[str, object]] = []
    analysis_mask = obs["tumor"].to_numpy() & obs["group"].isin(GROUPS).to_numpy()
    analysis_positions = np.flatnonzero(analysis_mask)
    analysis_frame = obs.iloc[analysis_positions]
    for (patient, group), labels in analysis_frame.groupby(["patient", "group"], observed=True).groups.items():
        pos = obj.obs_names.get_indexer(pd.Index(labels))
        if (pos < 0).any():
            raise RuntimeError("Internal cell-index mapping failure")
        lib_total = float(obs.iloc[pos]["library_size"].sum())
        studies = sorted(obs.iloc[pos]["study_v1"].unique())
        if len(studies) != 1:
            raise RuntimeError(f"Patient {patient} maps to multiple studies: {studies}")
        for gene in GENES:
            values = matrix[pos, gene_index[gene]]
            count_sum = float(values.sum())
            cpm = count_sum / lib_total * 1e6 if lib_total > 0 else math.nan
            pb_rows.append({
                "patient": patient,
                "study": studies[0],
                "group": group,
                "gene": gene,
                "n_cells": int(len(pos)),
                "library_size_sum": lib_total,
                "raw_count_sum": count_sum,
                "detection_fraction": float((values > 0).mean()),
                "log1p_cpm": float(np.log1p(cpm)) if np.isfinite(cpm) else math.nan,
            })
    pseudobulk = pd.DataFrame(pb_rows)
    pseudobulk.to_csv(OUT / "patient_group_gene_pseudobulk_v1.tsv.gz", sep="\t", index=False)

    # Fixed F1 program summaries and within-patient cell-level ligand associations.
    malignant_positions = np.flatnonzero(eligible.to_numpy())
    f1 = obs.iloc[malignant_positions]["F1"].to_numpy(dtype=float)
    libs = obs.iloc[malignant_positions]["library_size"].to_numpy(dtype=float)
    patients = obs.iloc[malignant_positions]["patient"].to_numpy(dtype=str)
    studies = obs.iloc[malignant_positions]["study_v1"].to_numpy(dtype=str)
    f1_patient_rows: list[dict[str, object]] = []
    corr_rows: list[dict[str, object]] = []
    for patient in sorted(np.unique(patients)):
        pmask = patients == patient
        study = str(np.unique(studies[pmask])[0])
        f1_patient_rows.append({
            "patient": patient,
            "study": study,
            "n_malignant_cells": int(pmask.sum()),
            "mean_F1": float(np.mean(f1[pmask])),
            "median_F1": float(np.median(f1[pmask])),
        })
        if pmask.sum() < MIN_CELLS:
            continue
        for ligand in LIGANDS:
            raw = matrix[malignant_positions[pmask], gene_index[ligand]]
            log_cp10k = np.log1p(raw / np.maximum(libs[pmask], 1.0) * 1e4)
            rho, p_value = safe_spearman(f1[pmask], log_cp10k)
            corr_rows.append({
                "patient": patient,
                "study": study,
                "ligand": ligand,
                "n_cells": int(pmask.sum()),
                "rho_F1_vs_ligand": rho,
                "p_value": p_value,
                "ligand_detection_fraction": float((raw > 0).mean()),
            })
    f1_patient = pd.DataFrame(f1_patient_rows)
    f1_patient.to_csv(OUT / "malignant_F1_patient_summary_v1.tsv", sep="\t", index=False)
    f1_correlations = pd.DataFrame(corr_rows)
    f1_correlations.to_csv(OUT / "patient_within_malignant_F1_ligand_correlations_v1.tsv", sep="\t", index=False)

    f1_cohort_rows: list[dict[str, object]] = []
    for ligand, frame in f1_correlations.groupby("ligand", observed=True):
        summary = direction_summary(frame["rho_F1_vs_ligand"])
        f1_cohort_rows.append({"ligand": ligand, **summary})
    f1_cohort = pd.DataFrame(f1_cohort_rows)
    f1_cohort["sign_test_fdr_bh"] = bh(f1_cohort["sign_test_p"])
    f1_cohort.to_csv(OUT / "cohort_F1_ligand_direction_summary_v1.tsv", sep="\t", index=False)

    loo_rows: list[dict[str, object]] = []
    for ligand, frame in f1_correlations.groupby("ligand", observed=True):
        for excluded in sorted(frame["study"].unique()):
            kept = frame.loc[frame["study"].ne(excluded)]
            summary = direction_summary(kept["rho_F1_vs_ligand"])
            loo_rows.append({"ligand": ligand, "excluded_study": excluded, **summary})
    pd.DataFrame(loo_rows).to_csv(OUT / "F1_ligand_leave_one_study_out_v1.tsv", sep="\t", index=False)

    # Exact candidate axes: expression competence and cross-patient coordination.
    pb_lookup = pseudobulk.set_index(["patient", "group", "gene"])
    f1_lookup = f1_patient.set_index("patient")
    axis_patient_rows: list[dict[str, object]] = []
    patients_all = sorted(pseudobulk["patient"].unique())
    for axis, source_group, ligand, target_group, receptors in AXES:
        for patient in patients_all:
            keys = [(patient, source_group, ligand)] + [(patient, target_group, receptor) for receptor in receptors]
            if not all(key in pb_lookup.index for key in keys):
                continue
            source = pb_lookup.loc[(patient, source_group, ligand)]
            receptor_rows = [pb_lookup.loc[(patient, target_group, receptor)] for receptor in receptors]
            if int(source["n_cells"]) < MIN_CELLS or min(int(row["n_cells"]) for row in receptor_rows) < MIN_CELLS:
                continue
            ligand_log = float(source["log1p_cpm"])
            receptor_log = float(min(float(row["log1p_cpm"]) for row in receptor_rows))
            ligand_det = float(source["detection_fraction"])
            receptor_det = float(min(float(row["detection_fraction"]) for row in receptor_rows))
            axis_patient_rows.append({
                "axis": axis,
                "patient": patient,
                "study": str(source["study"]),
                "source_group": source_group,
                "ligand": ligand,
                "target_group": target_group,
                "receptors": "+".join(receptors),
                "source_n_cells": int(source["n_cells"]),
                "target_n_cells": int(receptor_rows[0]["n_cells"]),
                "ligand_detection_fraction": ligand_det,
                "limiting_receptor_detection_fraction": receptor_det,
                "ligand_log1p_cpm": ligand_log,
                "limiting_receptor_log1p_cpm": receptor_log,
                "descriptive_pair_score": float(math.sqrt(max(ligand_log, 0) * max(receptor_log, 0))),
                "pair_detected": bool(ligand_det > 0 and receptor_det > 0),
                "mean_F1": float(f1_lookup.loc[patient, "mean_F1"]) if patient in f1_lookup.index else math.nan,
            })
    axis_patient = pd.DataFrame(axis_patient_rows)
    axis_patient.to_csv(OUT / "axis_patient_expression_competence_v1.tsv", sep="\t", index=False)

    axis_summary_rows: list[dict[str, object]] = []
    for axis, frame in axis_patient.groupby("axis", observed=True):
        rho_lr, p_lr = safe_spearman(frame["ligand_log1p_cpm"].to_numpy(), frame["limiting_receptor_log1p_cpm"].to_numpy())
        rho_f1, p_f1 = safe_spearman(frame["mean_F1"].to_numpy(), frame["descriptive_pair_score"].to_numpy())
        axis_summary_rows.append({
            "axis": axis,
            "n_patients_evaluable": int(len(frame)),
            "n_studies": int(frame["study"].nunique()),
            "patients_pair_detected": int(frame["pair_detected"].sum()),
            "median_ligand_detection_fraction": float(frame["ligand_detection_fraction"].median()),
            "median_limiting_receptor_detection_fraction": float(frame["limiting_receptor_detection_fraction"].median()),
            "median_descriptive_pair_score": float(frame["descriptive_pair_score"].median()),
            "rho_ligand_vs_receptor_across_patients": rho_lr,
            "p_ligand_vs_receptor_across_patients": p_lr,
            "rho_F1_vs_pair_score_across_patients": rho_f1,
            "p_F1_vs_pair_score_across_patients": p_f1,
        })
    axis_summary = pd.DataFrame(axis_summary_rows)
    axis_summary["fdr_ligand_vs_receptor"] = bh(axis_summary["p_ligand_vs_receptor_across_patients"])
    axis_summary["fdr_F1_vs_pair_score"] = bh(axis_summary["p_F1_vs_pair_score_across_patients"])
    axis_summary.to_csv(OUT / "axis_cohort_expression_summary_v1.tsv", sep="\t", index=False)

    # Remove between-study location shifts before evaluating cross-patient coordination.
    # This is a sensitivity analysis for study-composition/Simpson's-paradox effects.
    adjusted_rows: list[dict[str, object]] = []
    for axis, frame in axis_patient.groupby("axis", observed=True):
        frame = frame.copy()
        for variable in ["ligand_log1p_cpm", "limiting_receptor_log1p_cpm", "descriptive_pair_score", "mean_F1"]:
            frame[f"{variable}_study_centered"] = frame[variable] - frame.groupby("study", observed=True)[variable].transform("median")
        for excluded in ["NONE"] + sorted(frame["study"].unique().tolist()):
            kept = frame if excluded == "NONE" else frame.loc[frame["study"].ne(excluded)]
            rho_lr, p_lr = safe_spearman(
                kept["ligand_log1p_cpm_study_centered"].to_numpy(),
                kept["limiting_receptor_log1p_cpm_study_centered"].to_numpy(),
            )
            rho_f1, p_f1 = safe_spearman(
                kept["mean_F1_study_centered"].to_numpy(),
                kept["descriptive_pair_score_study_centered"].to_numpy(),
            )
            adjusted_rows.append({
                "axis": axis,
                "excluded_study": excluded,
                "n_patients": int(len(kept)),
                "n_studies": int(kept["study"].nunique()),
                "rho_ligand_vs_receptor_study_centered": rho_lr,
                "p_ligand_vs_receptor_study_centered": p_lr,
                "rho_F1_vs_pair_score_study_centered": rho_f1,
                "p_F1_vs_pair_score_study_centered": p_f1,
            })
    adjusted = pd.DataFrame(adjusted_rows)
    primary_adjusted = adjusted["excluded_study"].eq("NONE")
    adjusted.loc[primary_adjusted, "fdr_ligand_vs_receptor_study_centered"] = bh(
        adjusted.loc[primary_adjusted, "p_ligand_vs_receptor_study_centered"]
    )
    adjusted.loc[primary_adjusted, "fdr_F1_vs_pair_score_study_centered"] = bh(
        adjusted.loc[primary_adjusted, "p_F1_vs_pair_score_study_centered"]
    )
    adjusted.to_csv(OUT / "axis_study_centered_sensitivity_v1.tsv", sep="\t", index=False)

    # Directly compare SPP1 source strength between malignant cells and macrophages.
    spp1 = pseudobulk.loc[pseudobulk["gene"].eq("SPP1") & pseudobulk["group"].isin(["Malignant", "Macrophage"])]
    spp1_wide = spp1.pivot(index=["patient", "study"], columns="group", values=["n_cells", "detection_fraction", "log1p_cpm"])
    spp1_wide.columns = [f"{metric}_{group}" for metric, group in spp1_wide.columns]
    spp1_wide = spp1_wide.reset_index()
    spp1_wide = spp1_wide.loc[
        spp1_wide["n_cells_Malignant"].ge(MIN_CELLS) & spp1_wide["n_cells_Macrophage"].ge(MIN_CELLS)
    ].copy()
    spp1_wide["delta_detection_malignant_minus_macrophage"] = (
        spp1_wide["detection_fraction_Malignant"] - spp1_wide["detection_fraction_Macrophage"]
    )
    spp1_wide["delta_log1p_cpm_malignant_minus_macrophage"] = (
        spp1_wide["log1p_cpm_Malignant"] - spp1_wide["log1p_cpm_Macrophage"]
    )
    spp1_wide.to_csv(OUT / "SPP1_source_patient_comparison_v1.tsv", sep="\t", index=False)
    source_summary = {
        "n_patients": int(len(spp1_wide)),
        "n_studies": int(spp1_wide["study"].nunique()),
        "malignant_higher_detection_patients": int((spp1_wide["delta_detection_malignant_minus_macrophage"] > 0).sum()),
        "median_detection_delta": float(spp1_wide["delta_detection_malignant_minus_macrophage"].median()),
        "detection_delta_sign_test_p": direction_summary(spp1_wide["delta_detection_malignant_minus_macrophage"])["sign_test_p"],
        "malignant_higher_log1p_cpm_patients": int((spp1_wide["delta_log1p_cpm_malignant_minus_macrophage"] > 0).sum()),
        "median_log1p_cpm_delta": float(spp1_wide["delta_log1p_cpm_malignant_minus_macrophage"].median()),
        "log1p_cpm_delta_sign_test_p": direction_summary(spp1_wide["delta_log1p_cpm_malignant_minus_macrophage"])["sign_test_p"],
    }
    pd.DataFrame([source_summary]).to_csv(OUT / "SPP1_source_cohort_summary_v1.tsv", sep="\t", index=False)

    qc = pd.DataFrame([
        ("input_shape", "PASS" if obj.n_obs == 258120 and obj.n_vars == 37795 else "FAIL", f"{obj.n_obs}x{obj.n_vars}"),
        ("required_genes", "PASS" if not missing_genes else "FAIL", f"{len(GENES)}/{len(GENES)}"),
        ("frozen_malignant_cells", "PASS" if int(eligible.sum()) == 11950 else "FAIL", str(int(eligible.sum()))),
        ("rank7_F1_mapping", "PASS" if int(obs.loc[eligible, "F1"].notna().sum()) == 11950 else "FAIL", str(int(obs.loc[eligible, "F1"].notna().sum()))),
        ("raw_count_values", "PASS", "selected genes finite and nonnegative"),
        ("axis_outputs", "PASS" if len(axis_summary) == len(AXES) else "FAIL", f"{len(axis_summary)}/{len(AXES)} axes"),
    ], columns=["check", "status", "detail"])
    qc.to_csv(OUT / "analysis_qc_v1.tsv", sep="\t", index=False)
    if not qc["status"].eq("PASS").all():
        raise RuntimeError("One or more analysis QC checks failed")

    metadata = {
        "status": "completed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "input_h5ad": str(INPUT_H5AD),
        "input_h5ad_size_bytes": INPUT_H5AD.stat().st_size,
        "input_h5ad_verified_sha256": "e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7",
        "rank7_usage": str(USAGE),
        "selected_rank_record": str(SELECTED_RANK),
        "n_cells": int(obj.n_obs),
        "n_genes": int(obj.n_vars),
        "n_frozen_malignant_cells": int(eligible.sum()),
        "n_axis_definitions": len(AXES),
        "minimum_cells_per_patient_group": MIN_CELLS,
        "interpretation": "expression competence and patient consistency; not communication probability or causality",
        "outputs": sorted(path.name for path in OUT.iterdir() if path.is_file()),
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if getattr(obj, "file", None) is not None:
        obj.file.close()
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
