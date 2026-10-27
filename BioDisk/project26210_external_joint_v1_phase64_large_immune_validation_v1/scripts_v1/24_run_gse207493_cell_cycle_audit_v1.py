#!/usr/bin/env python3
"""Patient-aware cell-cycle confounding audit for GSE207493 F1/TYMP findings."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

for variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "8")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


SEED = 26210
EXPECTED_SAMPLES = [
    "RCC81", "RCC84", "RCC86", "RCC87", "RCC94", "RCC96", "RCC99",
    "RCC100", "RCC101", "RCC103", "RCC104", "RCC106", "RCC112",
    "RCC113", "RCC114", "RCC115", "RCC116", "RCC119", "RCC120",
]
MIN_CELLS = {
    "RCC81": 8, "RCC84": 7, "RCC86": 9, "RCC87": 9, "RCC94": 8,
    "RCC96": 8, "RCC99": 9, "RCC100": 10, "RCC101": 10,
    "RCC103": 10, "RCC104": 10, "RCC106": 10, "RCC112": 10,
    "RCC113": 9, "RCC114": 7, "RCC115": 9, "RCC116": 6,
    "RCC119": 10, "RCC120": 10,
}
MAX_FEATURES = {
    "RCC81": 3500, "RCC84": 4700, "RCC86": 5000, "RCC87": 2900,
    "RCC94": 3700, "RCC96": 3300, "RCC99": 2800, "RCC100": 3800,
    "RCC101": 2400, "RCC103": 3300, "RCC104": 5000,
    "RCC106": 4000, "RCC112": 4000, "RCC113": 4500,
    "RCC114": 3200, "RCC115": 2400, "RCC116": 4500,
    "RCC119": 4000, "RCC120": 3500,
}
T_TYPES = {"Treg", "Exhausted_T", "Proliferating_T", "CD8_T", "CD4_T"}


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def bh_adjust(values: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.notna() & np.isfinite(values)
    p = values.loc[finite].astype(float)
    if p.empty:
        return output
    order = np.argsort(p.to_numpy())
    ranked = p.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    output.loc[p.index] = restored
    return output


def find_matrices(root: Path) -> dict[str, tuple[Path, str]]:
    candidates: dict[str, list[tuple[Path, str]]] = {}
    for matrix in root.rglob("*matrix.mtx*"):
        if not matrix.is_file():
            continue
        match = re.search(r"RCC\d+", str(matrix), flags=re.IGNORECASE)
        if match is None:
            continue
        sample = match.group(0).upper()
        prefix = matrix.name.split("matrix.mtx", 1)[0]
        candidates.setdefault(sample, []).append((matrix.parent, prefix))
    chosen = {}
    for sample, entries in candidates.items():
        entries.sort(key=lambda item: (
            "filtered_feature_bc_matrix" not in str(item[0]), len(str(item[0]))
        ))
        chosen[sample] = entries[0]
    return chosen


def reconstruct(raw_root: Path) -> ad.AnnData:
    matrices = find_matrices(raw_root)
    missing = sorted(set(EXPECTED_SAMPLES) - set(matrices))
    if missing:
        raise RuntimeError(f"missing matrices: {missing}")
    objects = []
    for sample in EXPECTED_SAMPLES:
        directory, prefix = matrices[sample]
        one = sc.read_10x_mtx(directory, var_names="gene_symbols", make_unique=True, prefix=prefix)
        sc.pp.filter_genes(one, min_cells=MIN_CELLS[sample])
        one.var["mt"] = one.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(one, qc_vars=["mt"], inplace=True, percent_top=None)
        keep = (
            (one.obs["n_genes_by_counts"] > 500)
            & (one.obs["n_genes_by_counts"] < MAX_FEATURES[sample])
            & (one.obs["pct_counts_mt"] < 10)
        )
        one = one[keep].copy()
        one.obs_names = [f"{sample}_{barcode}" for barcode in one.obs_names]
        one.obs["sample"] = sample
        objects.append(one)
    merged = ad.concat(objects, join="outer", merge="same", fill_value=0)
    merged.obs["sample"] = pd.Categorical(merged.obs["sample"], categories=EXPECTED_SAMPLES)
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    return merged


def signed_rank(values: pd.Series) -> tuple[float, float]:
    values = values.dropna().astype(float)
    if len(values) < 5 or not np.any(values != 0):
        return np.nan, np.nan
    result = stats.wilcoxon(values, alternative="two-sided")
    return float(result.statistic), float(result.pvalue)


def residualize(y: np.ndarray, s_score: np.ndarray, g2m_score: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(y)), s_score, g2m_score])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    return y - design @ coefficients


def contrast_masks(cell_type: pd.Series, contrast: str) -> tuple[np.ndarray, np.ndarray]:
    if contrast == "F1_Tumor_vs_nonTumor":
        return (cell_type == "Tumor").to_numpy(), (cell_type != "Tumor").to_numpy()
    if contrast == "TYMP_MyeloidDC_vs_other":
        return cell_type.isin(["Myeloid", "Dendritic"]).to_numpy(), (
            ~cell_type.isin(["Myeloid", "Dendritic"])
        ).to_numpy()
    if contrast in {"TYMP_Treg_vs_otherT", "F1_Treg_vs_otherT"}:
        return (cell_type == "Treg").to_numpy(), (
            cell_type.isin(T_TYPES) & (cell_type != "Treg")
        ).to_numpy()
    raise ValueError(contrast)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_root", type=Path)
    parser.add_argument("cycle_genes_tsv", type=Path)
    parser.add_argument("prior_cell_metadata_tsv_gz", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 8
    np.random.seed(SEED)

    merged = reconstruct(args.raw_root)
    prior = pd.read_csv(args.prior_cell_metadata_tsv_gz, sep="\t")
    prior = prior.set_index("cell_barcode", verify_integrity=True)
    missing_metadata = merged.obs_names.difference(prior.index)
    extra_metadata = prior.index.difference(merged.obs_names)
    if len(missing_metadata) or len(extra_metadata):
        raise RuntimeError(
            f"cell metadata mismatch: missing={len(missing_metadata)}, extra={len(extra_metadata)}"
        )
    prior = prior.loc[merged.obs_names]
    for column in ["leiden", "cell_type", "f1_score", "TYMP_log1p"]:
        merged.obs[column] = prior[column].to_numpy()

    cycle = pd.read_csv(args.cycle_genes_tsv, sep="\t")
    s_all = cycle.loc[cycle.phase == "S", "gene"].astype(str).tolist()
    g2m_all = cycle.loc[cycle.phase == "G2M", "gene"].astype(str).tolist()
    s_present = [gene for gene in s_all if gene in merged.var_names]
    g2m_present = [gene for gene in g2m_all if gene in merged.var_names]
    sc.tl.score_genes_cell_cycle(
        merged, s_genes=s_present, g2m_genes=g2m_present, random_state=SEED
    )
    merged.obs["cycling_score"] = merged.obs[["S_score", "G2M_score"]].max(axis=1)

    phase = merged.obs.groupby(["cell_type", "phase"], observed=True).size().rename("n_cells").reset_index()
    phase["fraction_within_celltype"] = phase["n_cells"] / phase.groupby(
        "cell_type", observed=True
    )["n_cells"].transform("sum")
    write_tsv(phase, args.output_dir / "gse207493_cell_cycle_phase_by_celltype_v1.tsv")

    phase_sample = merged.obs.groupby(
        ["sample", "cell_type", "phase"], observed=True
    ).size().rename("n_cells").reset_index()
    phase_sample["fraction_within_sample_celltype"] = phase_sample["n_cells"] / phase_sample.groupby(
        ["sample", "cell_type"], observed=True
    )["n_cells"].transform("sum")
    write_tsv(phase_sample, args.output_dir / "gse207493_cell_cycle_phase_by_sample_celltype_v1.tsv")

    cycle_summary = merged.obs.groupby("cell_type", observed=True).agg(
        n_cells=("cell_type", "size"), n_samples=("sample", "nunique"),
        S_score_mean=("S_score", "mean"), G2M_score_mean=("G2M_score", "mean"),
        cycling_score_mean=("cycling_score", "mean"),
        cycling_fraction=("phase", lambda x: float(np.mean(x != "G1"))),
    ).reset_index()
    write_tsv(cycle_summary, args.output_dir / "gse207493_cell_cycle_summary_v1.tsv")

    contrast_features = {
        "F1_Tumor_vs_nonTumor": "f1_score",
        "TYMP_MyeloidDC_vs_other": "TYMP_log1p",
        "TYMP_Treg_vs_otherT": "TYMP_log1p",
        "F1_Treg_vs_otherT": "f1_score",
    }
    effect_rows = []
    correlation_rows = []
    for sample in EXPECTED_SAMPLES:
        sample_frame = merged.obs.loc[merged.obs["sample"] == sample].copy()
        s_values = sample_frame["S_score"].to_numpy(float)
        g2m_values = sample_frame["G2M_score"].to_numpy(float)
        for contrast, feature in contrast_features.items():
            target, reference = contrast_masks(sample_frame["cell_type"], contrast)
            if target.sum() < 20 or reference.sum() < 20:
                continue
            y = sample_frame[feature].to_numpy(float)
            residual = residualize(y, s_values, g2m_values)
            effect_rows.append({
                "sample": sample, "contrast": contrast, "feature": feature,
                "n_target": int(target.sum()), "n_reference": int(reference.sum()),
                "raw_mean_difference": float(y[target].mean() - y[reference].mean()),
                "cycle_adjusted_mean_difference": float(
                    residual[target].mean() - residual[reference].mean()
                ),
            })

        correlation_groups = {
            "Tumor": sample_frame.cell_type == "Tumor",
            "Myeloid_DC": sample_frame.cell_type.isin(["Myeloid", "Dendritic"]),
            "Treg": sample_frame.cell_type == "Treg",
            "other_T": sample_frame.cell_type.isin(T_TYPES) & (sample_frame.cell_type != "Treg"),
        }
        for group, mask in correlation_groups.items():
            if mask.sum() < 20:
                continue
            for feature in ["f1_score", "TYMP_log1p"]:
                for cycle_feature in ["S_score", "G2M_score", "cycling_score"]:
                    rho, p_value = stats.spearmanr(
                        sample_frame.loc[mask, feature], sample_frame.loc[mask, cycle_feature]
                    )
                    correlation_rows.append({
                        "sample": sample, "cell_group": group, "feature": feature,
                        "cycle_feature": cycle_feature, "n_cells": int(mask.sum()),
                        "rho": rho, "p_value_cell_level_descriptive": p_value,
                    })

    effects = pd.DataFrame(effect_rows)
    write_tsv(effects, args.output_dir / "gse207493_cell_cycle_adjusted_effects_by_sample_v1.tsv")
    effect_summary_rows = []
    for contrast, group in effects.groupby("contrast", observed=True):
        raw_stat, raw_p = signed_rank(group["raw_mean_difference"])
        adj_stat, adj_p = signed_rank(group["cycle_adjusted_mean_difference"])
        effect_summary_rows.append({
            "contrast": contrast, "feature": group["feature"].iloc[0],
            "n_samples": len(group),
            "median_raw_mean_difference": group["raw_mean_difference"].median(),
            "raw_wilcoxon_statistic": raw_stat, "raw_p_value": raw_p,
            "median_cycle_adjusted_mean_difference": group[
                "cycle_adjusted_mean_difference"
            ].median(),
            "adjusted_wilcoxon_statistic": adj_stat, "adjusted_p_value": adj_p,
        })
    effect_summary = pd.DataFrame(effect_summary_rows)
    effect_summary["raw_q_value"] = bh_adjust(effect_summary["raw_p_value"])
    effect_summary["adjusted_q_value"] = bh_adjust(effect_summary["adjusted_p_value"])
    write_tsv(effect_summary, args.output_dir / "gse207493_cell_cycle_adjusted_effect_summary_v1.tsv")

    correlations = pd.DataFrame(correlation_rows)
    write_tsv(correlations, args.output_dir / "gse207493_feature_cycle_correlations_by_sample_v1.tsv")
    correlation_summary_rows = []
    for keys, group in correlations.groupby(
        ["cell_group", "feature", "cycle_feature"], observed=True
    ):
        statistic, p_value = signed_rank(group["rho"])
        correlation_summary_rows.append({
            "cell_group": keys[0], "feature": keys[1], "cycle_feature": keys[2],
            "n_samples": len(group), "median_spearman_rho": group["rho"].median(),
            "wilcoxon_statistic": statistic, "p_value": p_value,
        })
    correlation_summary = pd.DataFrame(correlation_summary_rows)
    correlation_summary["q_value"] = bh_adjust(correlation_summary["p_value"])
    write_tsv(
        correlation_summary,
        args.output_dir / "gse207493_feature_cycle_correlation_summary_v1.tsv",
    )

    cell_output = merged.obs.reset_index().rename(columns={"index": "cell_barcode"})[
        ["cell_barcode", "sample", "leiden", "cell_type", "phase", "S_score",
         "G2M_score", "cycling_score", "f1_score", "TYMP_log1p"]
    ]
    write_tsv(cell_output, args.output_dir / "gse207493_cell_cycle_scores_v1.tsv.gz")

    gates = pd.DataFrame([
        {"gate": "qc_cells", "observed": merged.n_obs, "target": 102723,
         "status": "PASS" if merged.n_obs == 102723 else "FAIL",
         "note": "exact author-reported QC total"},
        {"gate": "cell_metadata_matched", "observed": merged.n_obs - len(missing_metadata),
         "target": merged.n_obs, "status": "PASS" if not len(missing_metadata) else "FAIL",
         "note": "one-to-one barcode match to frozen cell annotations"},
        {"gate": "S_genes_present", "observed": len(s_present), "target": len(s_all),
         "status": "PASS" if len(s_present) / len(s_all) >= 0.8 else "FAIL",
         "note": "Seurat::cc.genes.updated.2019"},
        {"gate": "G2M_genes_present", "observed": len(g2m_present), "target": len(g2m_all),
         "status": "PASS" if len(g2m_present) / len(g2m_all) >= 0.8 else "FAIL",
         "note": "Seurat::cc.genes.updated.2019"},
        {"gate": "phase_assigned", "observed": int(merged.obs.phase.notna().sum()),
         "target": merged.n_obs,
         "status": "PASS" if merged.obs.phase.notna().all() else "FAIL",
         "note": "Scanpy score_genes_cell_cycle"},
    ])
    write_tsv(gates, args.output_dir / "gse207493_cell_cycle_data_gate_v1.tsv")
    print(gates.to_string(index=False))
    print(cycle_summary.sort_values("cycling_fraction", ascending=False).to_string(index=False))
    print(effect_summary.to_string(index=False))


if __name__ == "__main__":
    main()
