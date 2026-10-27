#!/usr/bin/env python3
"""Sample-level single-cell support for the POSTN-CAF/APOE-TAM mechanism."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import binomtest, rankdata, wilcoxon


ROOT = Path("/mnt/BioDisk/ludingyu/project26210")
INPUT = ROOT / "results/part4_spatial/ref/ref_signatures.h5ad"
SCRIPT = ROOT / "results/part4_spatial/analyze_postn_apoe_singlecell_support_v1.py"
OUT = ROOT / "results/part4_spatial/postn_apoe_singlecell_support_v1"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
QC = OUT / "qc"
MEMORY_GATE = ROOT / "check_memory_gate.sh"
MIN_CELLS = 10

SIGNATURES = {
    "POSTN_ECM_proxy": [
        "POSTN", "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2",
        "COL5A1", "COL5A2", "FN1", "FAP", "PDGFRA", "LOX", "LOXL2",
        "MMP2", "MMP14", "THBS1", "TGFBI",
    ],
    "CAF_recruitment_ligands": ["POSTN", "CCL2", "CSF1"],
    "TGFb_receptors": ["TGFBR1", "TGFBR2"],
    "APOE_TAM_expression_proxy": [
        "APOE", "C1QA", "C1QB", "C1QC", "CD68", "CD163", "MRC1",
        "LGALS3", "CTSB", "LPL",
    ],
    "TAM_recruitment_receptors": ["ITGAV", "ITGB3", "CCR2", "CSF1R"],
    "TAM_feedback_ligands": ["TGFB1", "VEGFA", "PGF"],
}

GENES = sorted(set(
    [gene for values in SIGNATURES.values() for gene in values]
    + ["TGFB1", "VEGFA", "POSTN", "CCL2", "CSF1", "ITGAV", "ITGB3", "CCR2", "CSF1R"]
))

PRIMARY_COMPARISONS = [
    ("C1_fibroblast_POSTN_ECM", "fibroblast_candidate", "other_SMC_Pericyte", "POSTN_ECM_proxy", 1),
    ("C2_fibroblast_recruitment_ligands", "fibroblast_candidate", "other_SMC_Pericyte", "CAF_recruitment_ligands", 1),
    ("C3_fibroblast_TGFb_receptors", "fibroblast_candidate", "other_SMC_Pericyte", "TGFb_receptors", 1),
    ("C4_APOE_TAM_identity", "APOE_TAM", "other_macrophage", "APOE_TAM_expression_proxy", 1),
    ("C5_APOE_TAM_recruitment_receptors", "APOE_TAM", "other_macrophage", "TAM_recruitment_receptors", 1),
    ("C6_APOE_TAM_feedback_ligands", "APOE_TAM", "other_macrophage", "TAM_feedback_ligands", 1),
    ("C7_malignant_TGFB1_source", "malignant", "normal_epithelial", "TGFB1", 1),
]

GENE_COMPARISONS = [
    ("fibroblast_candidate", "other_SMC_Pericyte", gene)
    for gene in ["POSTN", "CCL2", "CSF1", "TGFBR1", "TGFBR2", "COL1A1", "FN1", "FAP"]
] + [
    ("APOE_TAM", "other_macrophage", gene)
    for gene in ["APOE", "C1QC", "ITGAV", "ITGB3", "CCR2", "CSF1R", "TGFB1", "VEGFA", "PGF"]
]


def memory_gate(stage: str) -> str:
    result = subprocess.run([str(MEMORY_GATE)], capture_output=True, text=True, check=False)
    line = (result.stdout or result.stderr).strip()
    print(f"[memory_gate] {stage} {line}", flush=True)
    if result.returncode != 0:
        raise RuntimeError(line)
    return line


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    answer = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    if not valid.any():
        return answer
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    corrected = np.empty_like(ranked)
    corrected[order] = np.minimum(ranked, 1.0)
    answer[valid] = corrected
    return answer


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = values.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(values)
    return (values - values.mean()) / sd


def donor_id(study: str, sample_id: str) -> str | None:
    if study == "GSE156632":
        match = re.match(r"RCC(\d+)[tn]$", sample_id)
        return f"{study}_RCC{match.group(1)}" if match else None
    if study == "GSE242299":
        match = re.match(r"[NT](\d+)_", sample_id)
        return f"{study}_{match.group(1)}" if match else None
    return None


def paired_test(deltas: np.ndarray) -> dict:
    deltas = np.asarray(deltas, dtype=float)
    deltas = deltas[np.isfinite(deltas)]
    nonzero = deltas[deltas != 0]
    if len(nonzero) == 0:
        return {"median_delta": 0.0, "wilcoxon_p": 1.0, "sign_p": 1.0, "n_positive": 0, "n_units": len(deltas)}
    p_wilcoxon = float(wilcoxon(nonzero, alternative="two-sided", method="auto").pvalue)
    n_positive = int(np.sum(nonzero > 0))
    p_sign = float(binomtest(n_positive, len(nonzero), 0.5, alternative="two-sided").pvalue)
    return {
        "median_delta": float(np.median(deltas)),
        "wilcoxon_p": p_wilcoxon,
        "sign_p": p_sign,
        "n_positive": n_positive,
        "n_units": len(deltas),
    }


def group_masks(obs: pd.DataFrame) -> dict[str, np.ndarray]:
    cell_type = obs["cell_type"].astype(str)
    raw = obs["celltypist_label_raw"].astype(str)
    c2l = obs["c2l_celltype"].astype(str)
    cellchat = obs["cellchat_group"].astype(str)
    fibroblast = (cell_type == "SMC/Pericyte") & (raw == "Fibroblasts")
    apoe = c2l == "APOE+C1QC+ TAM"
    return {
        "fibroblast_candidate": fibroblast.to_numpy(),
        "other_SMC_Pericyte": ((cell_type == "SMC/Pericyte") & ~fibroblast).to_numpy(),
        "APOE_TAM": apoe.to_numpy(),
        "other_macrophage": ((cell_type == "Macrophage") & ~apoe).to_numpy(),
        "malignant": cellchat.isin(["TYMP-High Malignant", "TYMP-Low Malignant"]).to_numpy(),
        "normal_epithelial": (cellchat == "Epithelial-Normal").to_numpy(),
    }


def build_pseudobulk(obs: pd.DataFrame, values: pd.DataFrame, masks: dict[str, np.ndarray]) -> pd.DataFrame:
    rows = []
    sample_ids = obs["sample_id"].astype(str).to_numpy()
    studies = obs["study"].astype(str).to_numpy()
    for group, group_mask in masks.items():
        subset_samples = np.unique(sample_ids[group_mask])
        for sample_id in subset_samples:
            mask = group_mask & (sample_ids == sample_id)
            if mask.sum() < MIN_CELLS:
                continue
            study_values = np.unique(studies[mask])
            if len(study_values) != 1:
                raise ValueError(f"Multiple studies for {sample_id}")
            record = {
                "study": study_values[0],
                "sample_id": sample_id,
                "group": group,
                "n_cells": int(mask.sum()),
            }
            record.update({column: float(values.loc[mask, column].mean()) for column in values.columns})
            rows.append(record)
    return pd.DataFrame(rows)


def paired_values(pseudobulk: pd.DataFrame, group_a: str, group_b: str, variable: str) -> pd.DataFrame:
    frame = pseudobulk[pseudobulk["group"].isin([group_a, group_b])]
    wide = frame.pivot(index=["study", "sample_id"], columns="group", values=variable).dropna()
    wide = wide.reset_index()
    wide["delta"] = wide[group_a] - wide[group_b]
    return wide


def malignant_matched_donors(pseudobulk: pd.DataFrame, variable: str) -> pd.DataFrame:
    frame = pseudobulk[pseudobulk["group"].isin(["malignant", "normal_epithelial"])].copy()
    frame["donor_id"] = [donor_id(study, sample) for study, sample in zip(frame["study"], frame["sample_id"])]
    frame = frame.dropna(subset=["donor_id"])
    donor_values = frame.groupby(["study", "donor_id", "group"], as_index=False)[variable].mean()
    wide = donor_values.pivot(index=["study", "donor_id"], columns="group", values=variable).dropna().reset_index()
    wide["delta"] = wide["malignant"] - wide["normal_epithelial"]
    return wide


def main() -> None:
    for directory in [TABLES, FIGURES, QC]:
        directory.mkdir(parents=True, exist_ok=True)
    memory_lines = [memory_gate("start")]

    adata = ad.read_h5ad(INPUT, backed="r")
    obs = adata.obs[["study", "sample_id", "Disease", "cell_type", "celltypist_label_raw", "c2l_celltype", "cellchat_group", "nCount_RNA"]].copy()
    present = [gene for gene in GENES if gene in adata.var_names]
    matrix = adata[:, present].to_memory().X
    if issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    totals = obs["nCount_RNA"].to_numpy(float)
    if np.any(totals <= 0):
        raise ValueError("Non-positive nCount_RNA")
    log_cpm = np.log1p(matrix / totals[:, None] * 1e4)
    expression = pd.DataFrame(log_cpm, columns=present)
    adata.file.close()

    values = pd.DataFrame(index=obs.index)
    for gene in present:
        values[gene] = zscore(expression[gene].to_numpy())
    signature_qc = []
    for signature, requested in SIGNATURES.items():
        genes = [gene for gene in requested if gene in values.columns]
        if len(genes) < 2:
            raise ValueError(f"Too few genes for {signature}: {genes}")
        values[signature] = values[genes].mean(axis=1).to_numpy()
        signature_qc.append(
            {
                "signature": signature,
                "requested": ",".join(requested),
                "present": ",".join(genes),
                "missing": ",".join(sorted(set(requested) - set(genes))),
                "n_present": len(genes),
            }
        )

    masks = group_masks(obs)
    group_qc = []
    for group, mask in masks.items():
        group_qc.append(
            {
                "group": group,
                "n_cells": int(mask.sum()),
                "n_samples": int(obs.loc[mask, "sample_id"].nunique()),
                "n_studies": int(obs.loc[mask, "study"].nunique()),
                "study_counts": ";".join(f"{key}:{value}" for key, value in obs.loc[mask, "study"].value_counts().items()),
            }
        )

    pseudobulk = build_pseudobulk(obs, values, masks)
    primary_rows = []
    paired_detail = []
    for comparison, group_a, group_b, variable, expected in PRIMARY_COMPARISONS:
        if comparison == "C7_malignant_TGFB1_source":
            wide = malignant_matched_donors(pseudobulk, variable)
            unit = "matched_donor"
            unit_id = "donor_id"
        else:
            wide = paired_values(pseudobulk, group_a, group_b, variable)
            unit = "sample"
            unit_id = "sample_id"
        test = paired_test(wide["delta"].to_numpy())
        primary_rows.append(
            {
                "comparison": comparison,
                "group_a": group_a,
                "group_b": group_b,
                "variable": variable,
                "unit": unit,
                "expected_sign": expected,
                **test,
            }
        )
        for row in wide.to_dict("records"):
            paired_detail.append(
                {
                    "comparison": comparison,
                    "study": row["study"],
                    "unit_id": row[unit_id],
                    "delta": row["delta"],
                }
            )

    primary = pd.DataFrame(primary_rows)
    primary["bh_q"] = bh(primary["wilcoxon_p"].to_numpy())
    primary["direction_matches"] = np.sign(primary["median_delta"]) == primary["expected_sign"]
    primary["support_class"] = np.where(
        primary["direction_matches"]
        & (primary["bh_q"] < 0.05)
        & (primary["n_positive"] >= np.ceil(2 * primary["n_units"] / 3)),
        "sample_level_supported",
        np.where(primary["direction_matches"], "direction_only", "not_supported"),
    )
    paired_detail = pd.DataFrame(paired_detail)

    gene_rows = []
    for group_a, group_b, gene in GENE_COMPARISONS:
        if gene not in pseudobulk.columns:
            gene_rows.append(
                {
                    "group_a": group_a,
                    "group_b": group_b,
                    "gene": gene,
                    "median_delta": np.nan,
                    "wilcoxon_p": np.nan,
                    "sign_p": np.nan,
                    "n_positive": 0,
                    "n_units": 0,
                    "availability": "gene_absent_from_matrix",
                }
            )
            continue
        wide = paired_values(pseudobulk, group_a, group_b, gene)
        test = paired_test(wide["delta"].to_numpy())
        gene_rows.append({"group_a": group_a, "group_b": group_b, "gene": gene, **test, "availability": "available"})
    gene_results = pd.DataFrame(gene_rows)
    gene_results["bh_q"] = bh(gene_results["wilcoxon_p"].to_numpy())

    study_rows = []
    for comparison, frame in paired_detail.groupby("comparison"):
        for study, study_frame in frame.groupby("study"):
            if len(study_frame) < 3:
                continue
            result = paired_test(study_frame["delta"].to_numpy())
            study_rows.append({"comparison": comparison, "study": study, **result})
    study_results = pd.DataFrame(study_rows)

    pseudobulk.to_csv(TABLES / "sample_group_pseudobulk.tsv", sep="\t", index=False)
    primary.to_csv(TABLES / "primary_sample_level_tests.tsv", sep="\t", index=False)
    paired_detail.to_csv(TABLES / "paired_effects.tsv", sep="\t", index=False)
    gene_results.to_csv(TABLES / "secondary_gene_level_tests.tsv", sep="\t", index=False)
    study_results.to_csv(TABLES / "study_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(signature_qc).to_csv(TABLES / "signature_gene_qc.tsv", sep="\t", index=False)
    pd.DataFrame(group_qc).to_csv(TABLES / "group_qc.tsv", sep="\t", index=False)
    (QC / "memory_gate.txt").write_text("\n".join(memory_lines) + "\n", encoding="utf-8")

    order = [item[0] for item in PRIMARY_COMPARISONS]
    data = [paired_detail.loc[paired_detail["comparison"] == comparison, "delta"].to_numpy() for comparison in order]
    fig, ax = plt.subplots(figsize=(11, 6))
    box = ax.boxplot(data, vert=False, labels=[value.split("_", 1)[0] for value in order], showfliers=False, patch_artist=True)
    for patch in box["boxes"]:
        patch.set_facecolor("#d7e8f5")
    rng = np.random.default_rng(2621005)
    for y, values_array in enumerate(data, start=1):
        ax.scatter(values_array, y + rng.uniform(-0.12, 0.12, len(values_array)), s=13, alpha=0.55, color="#2c5f8a")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Within-unit mean signature/expression difference (group A minus group B)")
    ax.set_title("Single-cell reference: sample/donor-level paired mechanism tests")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig01_primary_paired_effects.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig01_primary_paired_effects.pdf", bbox_inches="tight")
    plt.close(fig)

    selected_groups = ["fibroblast_candidate", "other_SMC_Pericyte", "APOE_TAM", "other_macrophage", "malignant", "normal_epithelial"]
    selected_genes = [
        gene for gene in
        ["POSTN", "COL1A1", "CCL2", "CSF1", "TGFBR1", "TGFBR2", "APOE", "C1QC", "ITGAV", "ITGB3", "CCR2", "CSF1R", "TGFB1", "VEGFA"]
        if gene in values.columns
    ]
    group_means = []
    for group in selected_groups:
        mask = masks[group]
        group_means.append([float(values.loc[mask, gene].mean()) for gene in selected_genes])
    group_means = np.asarray(group_means)
    fig, ax = plt.subplots(figsize=(13, 5))
    image = ax.imshow(group_means, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(selected_genes)), selected_genes, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(selected_groups)), selected_groups)
    ax.set_title("Cell-level standardized expression by analysis group (descriptive)")
    fig.colorbar(image, ax=ax, label="Mean gene z-score")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig02_group_gene_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / "fig02_group_gene_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    n_supported = int((primary["support_class"] == "sample_level_supported").sum())
    summary = {
        "input": str(INPUT),
        "input_size_bytes": INPUT.stat().st_size,
        "input_sha256": sha256(INPUT),
        "script": str(SCRIPT),
        "script_sha256": sha256(SCRIPT),
        "n_cells": len(obs),
        "n_studies": int(obs["study"].nunique()),
        "minimum_cells_per_group_per_sample": MIN_CELLS,
        "n_primary_comparisons": len(primary),
        "n_sample_level_supported": n_supported,
        "limitations": [
            "Fibroblast candidates use the raw CellTypist fibroblast label inside the final SMC/Pericyte compartment; they are not a newly clustered POSTN+ CAF state.",
            "APOE+C1QC+ TAM cells are strongly concentrated in GSE242299, so study composition is reported and paired within-sample effects are primary.",
            "Expression competence and group enrichment do not reproduce intervention-based causal evidence.",
        ],
    }
    (OUT / "run_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    readme = f"""# POSTN-CAF/APOE-TAM 单细胞支持分析 v1

本分析只读使用 cell2location 参考对象中的 230,618 个细胞，以样本或匹配供者为独立单位检验论文机制所需的细胞表达能力。7 个预设比较中有 {n_supported} 个达到样本级支持标准；数值见 `tables/primary_sample_level_tests.tsv`。

纤维细胞候选定义为最终 SMC/Pericyte 区室中、原始 CellTypist 标签为 Fibroblasts 的细胞；这不是重新聚类得到的 POSTN⁺CAF 亚群。APOE+C1QC+ TAM 明显集中于 GSE242299，因此采用样本内配对比较并保留研究分层敏感性表。该结果仍是观察性表达支持，不能替代论文中的中和抗体、共培养或小鼠干预证据。
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(primary.to_string(index=False), flush=True)
    print(f"[done] {OUT}", flush=True)


if __name__ == "__main__":
    main()
