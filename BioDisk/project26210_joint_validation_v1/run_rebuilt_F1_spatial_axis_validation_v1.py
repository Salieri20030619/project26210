#!/usr/bin/env python3
"""Spatially validate rebuilt Part 1 F1-linked candidate axes in 12 TJ-RCC cases.

The analysis projects an axis-gene-independent F1 signature into Visium spots,
uses broad receiver abundance only, and treats patients as biological replicates.
It tests spatial association, not direct communication or causality.
"""

from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.spatial import cKDTree
from scipy.stats import binomtest, rankdata, spearmanr


MANIFEST = Path(
    "/mnt/BioDisk/ludingyu/project26210/external/ccRCC_spatial/TJ_RCC_figshare_24599295/"
    "analysis_cell2location_full_v1/input_manifest.tsv"
)
F1_TOP = Path(
    "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/05_cnmf/"
    "cnv_clone_program_function_v1/rank7_program_top100_genes_v1.tsv"
)
OUT = Path(
    "/mnt/BioDisk/wangmingyuan/project26210_joint_validation_v1/"
    "rebuilt_F1_spatial_axis_validation_v1"
)
MEMORY_GATE = Path("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")

SEED = 20260805
N_NEIGHBORS = 6
N_PERM = 9999
GRIDS = (4, 5)

AXES = [
    ("VEGFA_KDR", "VEGFA", ["KDR"], "Endothelial"),
    ("VEGFA_FLT1", "VEGFA", ["FLT1"], "Endothelial"),
    ("VEGFB_FLT1", "VEGFB", ["FLT1"], "Endothelial"),
    ("PGF_FLT1", "PGF", ["FLT1"], "Endothelial"),
    ("APP_CD74", "APP", ["CD74"], "Macrophage"),
    ("APOE_TREM2_TYROBP", "APOE", ["TREM2", "TYROBP"], "Macrophage"),
    ("SPP1_CD44", "SPP1", ["CD44"], "Macrophage"),
    ("MIF_CD74_CXCR4", "MIF", ["CD74", "CXCR4"], "Macrophage"),
]
AXIS_GENES = sorted({ligand for _, ligand, _, _ in AXES} | {g for _, _, genes, _ in AXES for g in genes} | {"TYMP"})
TUMOR_CONTEXT = ["CA9", "EPCAM", "KRT8", "KRT18", "KRT19", "NDUFA4L2"]
ABUNDANCE_PREFIXES = ["q05cell_abundance_w_sf_means_per_cluster_mu_fg_", "q05cell_abundance_w_sf_"]


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


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = float(np.nanstd(values))
    if not np.isfinite(sd) or sd <= 0:
        return np.zeros_like(values)
    return (values - np.nanmean(values)) / sd


def percentile(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average", pct=True).to_numpy(dtype=float)


def residualize(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    design = [np.ones(len(y))] + [zscore(value) for value in covariates]
    x = np.column_stack(design)
    beta = np.linalg.lstsq(x, np.asarray(y, dtype=float), rcond=None)[0]
    return zscore(np.asarray(y, dtype=float) - x @ beta)


def memory_gate(stage: str, rows: list[dict[str, object]]) -> None:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip()
    rows.append({"stage": stage, "exit_code": result.returncode, "output": line})
    if result.returncode != 0:
        raise RuntimeError(f"Memory gate failed at {stage}: {line}")


def clean_abundance(adata: ad.AnnData) -> pd.DataFrame:
    frame = adata.obsm["q05_cell_abundance_w_sf"]
    if not isinstance(frame, pd.DataFrame):
        names = adata.uns.get("mod", {}).get("factor_names")
        frame = pd.DataFrame(frame, index=adata.obs_names, columns=names)
    frame = frame.copy()
    cleaned = []
    for column in frame.columns:
        value = str(column)
        for prefix in ABUNDANCE_PREFIXES:
            value = value.replace(prefix, "")
        cleaned.append(value)
    frame.columns = cleaned
    required = {"TYMP-High Malignant", "TYMP-Low Malignant", "Endothelial", "SPP1+ TAM", "APOE+C1QC+ TAM", "Other Macrophage"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"Missing q05 abundance columns: {missing}")
    return frame


def block_labels(coords: np.ndarray, grid: int) -> np.ndarray:
    bins = []
    for axis in range(2):
        ranks = pd.Series(coords[:, axis]).rank(method="first")
        bins.append(pd.qcut(ranks, q=grid, labels=False).to_numpy())
    return bins[0] * grid + bins[1]


def block_medians(x: np.ndarray, y: np.ndarray, labels: np.ndarray, mask: np.ndarray) -> pd.DataFrame:
    frame = pd.DataFrame({"block": labels[mask], "x": x[mask], "y": y[mask]})
    return frame.groupby("block", as_index=False).agg(x=("x", "median"), y=("y", "median"), n=("x", "size"))


def permutation_spearman(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float]:
    if len(x) < 6 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return math.nan, math.nan
    rho = float(spearmanr(x, y).statistic)
    xr = rankdata(x) - np.mean(rankdata(x))
    yr = rankdata(y) - np.mean(rankdata(y))
    denominator = float(np.sqrt(np.sum(xr**2) * np.sum(yr**2)))
    if denominator <= 0:
        return math.nan, math.nan
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(N_PERM):
        null = float(np.dot(xr, rng.permutation(yr)) / denominator)
        extreme += int(abs(null) >= abs(rho))
    return rho, float((1 + extreme) / (N_PERM + 1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(MANIFEST, sep="\t")
    if len(manifest) != 12 or manifest["sample"].duplicated().any():
        raise RuntimeError("Expected 12 unique TJ-RCC samples")

    top = pd.read_csv(F1_TOP, sep="\t")
    top = top.loc[top["factor"].eq("F1")].sort_values("gene_rank")
    top = top.loc[~top["gene"].isin(AXIS_GENES) & ~top["gene"].str.startswith("HLA-")]
    f1_signature = top.head(30)[["gene", "gene_rank", "loading"]].copy()
    if len(f1_signature) != 30:
        raise RuntimeError("Could not construct 30-gene independent F1 signature")
    f1_signature.to_csv(OUT / "F1_projection_signature_v1.tsv", sep="\t", index=False)

    requested = sorted(set(f1_signature["gene"]) | set(AXIS_GENES) | set(TUMOR_CONTEXT))
    patient_rows: list[dict[str, object]] = []
    score_rows: list[pd.DataFrame] = []
    gene_qc_rows: list[dict[str, object]] = []
    memory_rows: list[dict[str, object]] = []

    for sample_number, item in manifest.sort_values("sample").reset_index(drop=True).iterrows():
        sample = str(item["sample"])
        path = Path(item["path"])
        memory_gate(f"before_{sample}", memory_rows)
        adata = ad.read_h5ad(path)
        present = [gene for gene in requested if gene in adata.var_names]
        missing = sorted(set(requested) - set(present))
        if any(gene in missing for gene in AXIS_GENES):
            raise RuntimeError(f"{sample}: missing axis genes {sorted(set(AXIS_GENES) & set(missing))}")
        used_f1 = [gene for gene in f1_signature["gene"] if gene in present]
        if len(used_f1) < 24:
            raise RuntimeError(f"{sample}: only {len(used_f1)}/30 F1 signature genes available")

        matrix = adata[:, present].X
        if issparse(matrix):
            matrix = matrix.toarray()
        matrix = np.asarray(matrix, dtype=float)
        totals = np.asarray(adata.X.sum(axis=1)).reshape(-1).astype(float)
        if np.any(totals <= 0) or not np.isfinite(matrix).all():
            raise RuntimeError(f"{sample}: invalid expression matrix")
        log_cp10k = np.log1p(matrix / totals[:, None] * 1e4)
        expression = pd.DataFrame(log_cp10k, index=adata.obs_names, columns=present)
        abundance = clean_abundance(adata)
        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        if coords.shape != (adata.n_obs, 2):
            raise RuntimeError(f"{sample}: invalid spatial coordinates")

        weights = f1_signature.set_index("gene").loc[used_f1, "loading"].to_numpy(dtype=float)
        weights = weights / weights.sum()
        f1_score = np.column_stack([zscore(expression[gene].to_numpy()) for gene in used_f1]) @ weights
        tumor_context_genes = [gene for gene in TUMOR_CONTEXT if gene in expression]
        tumor_context = np.column_stack([zscore(expression[gene].to_numpy()) for gene in tumor_context_genes]).mean(axis=1)
        tumor_abundance = abundance["TYMP-High Malignant"].to_numpy() + abundance["TYMP-Low Malignant"].to_numpy()
        macrophage_abundance = abundance[["SPP1+ TAM", "APOE+C1QC+ TAM", "Other Macrophage"]].sum(axis=1).to_numpy()
        endothelial_abundance = abundance["Endothelial"].to_numpy()
        n_genes = adata.obs["n_genes_by_counts"].to_numpy(dtype=float)
        x, y = coords[:, 0], coords[:, 1]
        spatial_covariates = [totals, n_genes, x, y, x**2, y**2, x * y]
        f1_residual = residualize(f1_score, [tumor_context, tumor_abundance, *spatial_covariates])
        neighbors = np.asarray(cKDTree(coords).query(coords, k=N_NEIGHBORS + 1)[1][:, 1:], dtype=int)
        subsets = {
            "all_spots": np.ones(adata.n_obs, dtype=bool),
            "tumor_enriched": percentile(tumor_abundance) >= 0.5,
        }

        sample_scores = pd.DataFrame({
            "sample": sample,
            "spot": adata.obs_names.astype(str),
            "F1_projection_score": f1_score,
            "F1_projection_residual": f1_residual,
            "tumor_context_score": tumor_context,
            "tumor_abundance": tumor_abundance,
            "macrophage_abundance": macrophage_abundance,
            "endothelial_abundance": endothelial_abundance,
        })

        for axis_index, (axis, ligand, receptors, target) in enumerate(AXES):
            ligand_raw = expression[ligand].to_numpy(dtype=float)
            ligand_residual = residualize(ligand_raw, [tumor_abundance, *spatial_covariates])
            receptor_raw = np.column_stack([zscore(expression[gene].to_numpy()) for gene in receptors]).mean(axis=1)
            receptor_residual = residualize(receptor_raw, spatial_covariates)
            target_abundance = endothelial_abundance if target == "Endothelial" else macrophage_abundance
            receiver_competence = np.sqrt(percentile(target_abundance) * percentile(receptor_residual))
            neighbor_receiver = receiver_competence[neighbors].mean(axis=1)
            sample_scores[f"{axis}__ligand_residual"] = ligand_residual
            sample_scores[f"{axis}__neighbor_receiver_competence"] = neighbor_receiver

            for subset, mask in subsets.items():
                for grid in GRIDS:
                    labels = block_labels(coords, grid)
                    tests = {
                        "F1_to_neighbor_receiver": (f1_residual, neighbor_receiver),
                        "ligand_to_neighbor_receiver": (ligand_residual, neighbor_receiver),
                    }
                    for test_index, (test, (left, right)) in enumerate(tests.items()):
                        blocks = block_medians(left, right, labels, mask)
                        rho, p_value = permutation_spearman(
                            blocks["x"].to_numpy(), blocks["y"].to_numpy(),
                            SEED + sample_number * 10000 + axis_index * 100 + test_index * 10 + grid,
                        )
                        patient_rows.append({
                            "sample": sample,
                            "axis": axis,
                            "ligand": ligand,
                            "receptors": "+".join(receptors),
                            "target": target,
                            "subset": subset,
                            "grid": f"{grid}x{grid}",
                            "test": test,
                            "n_spots": int(mask.sum()),
                            "n_blocks": int(len(blocks)),
                            "rho": rho,
                            "permutation_p": p_value,
                        })
            gene_qc_rows.append({
                "sample": sample,
                "axis": axis,
                "ligand_detection_fraction": float((matrix[:, present.index(ligand)] > 0).mean()),
                "minimum_receptor_detection_fraction": float(min((matrix[:, present.index(gene)] > 0).mean() for gene in receptors)),
                "n_F1_signature_genes": len(used_f1),
            })
        score_rows.append(sample_scores)

    patient = pd.DataFrame(patient_rows)
    patient["within_patient_fdr_bh"] = patient.groupby("sample", observed=True)["permutation_p"].transform(bh)
    patient.to_csv(OUT / "patient_spatial_axis_tests_v1.tsv", sep="\t", index=False)
    pd.concat(score_rows, ignore_index=True).to_csv(OUT / "patient_spot_scores_v1.tsv.gz", sep="\t", index=False)
    pd.DataFrame(gene_qc_rows).to_csv(OUT / "spatial_gene_detection_qc_v1.tsv", sep="\t", index=False)
    pd.DataFrame(memory_rows).to_csv(OUT / "memory_gates_v1.tsv", sep="\t", index=False)

    cohort_rows: list[dict[str, object]] = []
    keys = ["axis", "ligand", "receptors", "target", "subset", "grid", "test"]
    for values, frame in patient.groupby(keys, observed=True):
        rhos = frame["rho"].dropna()
        positive = int((rhos > 0).sum())
        negative = int((rhos < 0).sum())
        nonzero = positive + negative
        cohort_rows.append({
            **dict(zip(keys, values)),
            "n_patients": int(len(rhos)),
            "positive_patients": positive,
            "negative_patients": negative,
            "median_rho": float(rhos.median()),
            "q25_rho": float(rhos.quantile(0.25)),
            "q75_rho": float(rhos.quantile(0.75)),
            "sign_test_p": float(binomtest(positive, nonzero, 0.5).pvalue) if nonzero else math.nan,
            "patients_within_fdr_lt_0_05": int((frame["within_patient_fdr_bh"] < 0.05).sum()),
        })
    cohort = pd.DataFrame(cohort_rows)
    cohort["sign_test_fdr_bh_global"] = bh(cohort["sign_test_p"])
    cohort.to_csv(OUT / "cohort_spatial_axis_direction_summary_v1.tsv", sep="\t", index=False)

    primary = cohort.loc[cohort["subset"].eq("tumor_enriched")].copy()
    decisions = []
    for axis, frame in primary.groupby("axis", observed=True):
        f1 = frame.loc[frame["test"].eq("F1_to_neighbor_receiver")]
        ligand = frame.loc[frame["test"].eq("ligand_to_neighbor_receiver")]
        f1_global = int((f1["positive_patients"].ge(9) & f1["sign_test_fdr_bh_global"].lt(0.05)).sum())
        ligand_global = int((ligand["positive_patients"].ge(9) & ligand["sign_test_fdr_bh_global"].lt(0.05)).sum())
        f1_nominal = int((f1["positive_patients"].ge(9) & f1["sign_test_p"].lt(0.05)).sum())
        ligand_nominal = int((ligand["positive_patients"].ge(9) & ligand["sign_test_p"].lt(0.05)).sum())
        if f1_global == 2 and ligand_global == 2:
            decision = "exact_axis_supported_both_grids"
        elif f1_global > 0 and ligand_global == 0:
            decision = "receiver_context_only_not_ligand_specific"
        elif ligand_global > 0 and f1_global == 0:
            decision = "ligand_context_only_not_F1_specific"
        elif f1_global > 0 or ligand_global > 0:
            decision = "component_support_not_exact_axis"
        elif f1_nominal > 0 and ligand_nominal > 0:
            decision = "suggestive_but_not_global_FDR"
        else:
            decision = "exact_axis_not_spatially_supported"
        decisions.append({
            "axis": axis,
            "F1_global_FDR_supported_grids": f1_global,
            "ligand_global_FDR_supported_grids": ligand_global,
            "F1_nominal_supported_grids": f1_nominal,
            "ligand_nominal_supported_grids": ligand_nominal,
            "decision": decision,
        })
    decision_table = pd.DataFrame(decisions).sort_values("axis")
    decision_table.to_csv(OUT / "primary_axis_decision_v1.tsv", sep="\t", index=False)

    metadata = {
        "status": "completed",
        "created_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_patients": 12,
        "n_axes": len(AXES),
        "n_neighbors": N_NEIGHBORS,
        "block_grids": list(GRIDS),
        "permutations": N_PERM,
        "F1_signature": "top 30 rank-7/F1 genes after excluding all tested axis genes, TYMP, and HLA genes",
        "primary_test": "tumor-enriched spots: rebuilt F1 residual versus spatial-lagged receiver competence",
        "inference_boundary": "spatial association and expression competence only; not direct communication or causality",
        "input_manifest": str(MANIFEST),
        "F1_top_gene_input": str(F1_TOP),
    }
    (OUT / "run_metadata_v1.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
