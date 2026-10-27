#!/usr/bin/env python3
"""Patient-aware metabolic heterogeneity analysis for GSE207493."""

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
from scipy import stats


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

PATHWAY_PANEL = {
    "R-HSA-500753": "nucleotide",
    "R-HSA-73621": "nucleotide",
    "R-HSA-73614": "nucleotide",
    "R-HSA-83936": "nucleotide",
    "R-HSA-74259": "nucleotide",
    "R-HSA-73817": "nucleotide",
    "R-HSA-74217": "nucleotide",
    "R-HSA-70171": "glucose_energy",
    "R-HSA-70263": "glucose_energy",
    "R-HSA-70268": "glucose_energy",
    "R-HSA-204174": "glucose_energy",
    "R-HSA-71403": "glucose_energy",
    "R-HSA-1428517": "glucose_energy",
    "R-HSA-611105": "glucose_energy",
    "R-HSA-71336": "glucose_energy",
    "R-HSA-8978868": "lipid",
    "R-HSA-77289": "lipid",
    "R-HSA-191273": "lipid",
    "R-HSA-1655829": "lipid",
    "R-HSA-174403": "redox_hypoxia",
    "R-HSA-3299685": "redox_hypoxia",
    "R-HSA-1234174": "redox_hypoxia",
    "R-HSA-1234158": "redox_hypoxia",
    "R-HSA-8964539": "amino_acid_cofactor",
    "R-HSA-71240": "amino_acid_cofactor",
    "R-HSA-977347": "amino_acid_cofactor",
    "R-HSA-1614635": "amino_acid_cofactor",
    "R-HSA-196757": "amino_acid_cofactor",
    "R-HSA-189451": "amino_acid_cofactor",
    "R-HSA-917937": "amino_acid_cofactor",
    "R-HSA-71291": "amino_acid_cofactor",
    "R-HSA-70895": "amino_acid_cofactor",
}


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


def signed_rank(values: pd.Series) -> tuple[float, float]:
    values = values.dropna().astype(float)
    if len(values) < 5 or not np.any(values != 0):
        return np.nan, np.nan
    result = stats.wilcoxon(values, alternative="two-sided")
    return float(result.statistic), float(result.pvalue)


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
        print(f"loaded {sample}: {one.n_obs} cells", flush=True)
    merged = ad.concat(objects, join="outer", merge="same", fill_value=0)
    merged.obs["sample"] = pd.Categorical(merged.obs["sample"], categories=EXPECTED_SAMPLES)
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    print(f"merged: {merged.n_obs} cells x {merged.n_vars} genes", flush=True)
    return merged


def parse_gmt(path: Path) -> dict[str, tuple[str, list[str]]]:
    output = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                output[fields[1]] = (fields[0], list(dict.fromkeys(fields[2:])))
    return output


def residualize_matrix(values: np.ndarray, s_score: np.ndarray, g2m_score: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(values)), s_score, g2m_score])
    coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def quantile_masks(values: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    low_cut, high_cut = np.quantile(values, [0.25, 0.75])
    if not np.isfinite(low_cut) or not np.isfinite(high_cut) or high_cut <= low_cut:
        return None
    low = values <= low_cut
    high = values >= high_cut
    if low.sum() < 5 or high.sum() < 5:
        return None
    return high, low


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_root", type=Path)
    parser.add_argument("reactome_gmt", type=Path)
    parser.add_argument("cell_cycle_metadata_tsv_gz", type=Path)
    parser.add_argument("signature_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 8
    np.random.seed(SEED)

    merged = reconstruct(args.raw_root)
    metadata = pd.read_csv(args.cell_cycle_metadata_tsv_gz, sep="\t")
    metadata = metadata.set_index("cell_barcode", verify_integrity=True)
    missing_metadata = merged.obs_names.difference(metadata.index)
    extra_metadata = metadata.index.difference(merged.obs_names)
    if len(missing_metadata) or len(extra_metadata):
        raise RuntimeError(
            f"cell metadata mismatch: missing={len(missing_metadata)}, extra={len(extra_metadata)}"
        )
    metadata = metadata.loc[merged.obs_names]
    for column in [
        "leiden", "cell_type", "phase", "S_score", "G2M_score",
        "cycling_score", "f1_score", "TYMP_log1p",
    ]:
        merged.obs[column] = metadata[column].to_numpy()

    signature = pd.read_csv(args.signature_tsv, sep="\t")
    f1_genes = set(
        signature.loc[
            signature.signature_id == "f1_top100_weighted", "gene"
        ].astype(str).str.upper()
    )
    reactome = parse_gmt(args.reactome_gmt)
    panel_rows = []
    score_columns = []
    f1_excluded_score_columns = []
    tymp_excluded_score_columns = []
    score_to_pathway = {}
    for index, (stable_id, category) in enumerate(PATHWAY_PANEL.items(), start=1):
        if stable_id not in reactome:
            panel_rows.append({
                "stable_id": stable_id, "pathway": "", "category": category,
                "genes_total": 0, "genes_present": 0, "coverage": 0,
                "score_column": "", "status": "FAIL_NOT_IN_GMT",
            })
            continue
        name, genes = reactome[stable_id]
        present = [gene for gene in genes if gene in merged.var_names]
        f1_excluded = [gene for gene in present if gene.upper() not in f1_genes]
        tymp_excluded = [gene for gene in present if gene.upper() != "TYMP"]
        score_column = f"METAB_{index:02d}"
        f1_excluded_column = f"{score_column}_F1X"
        tymp_excluded_column = f"{score_column}_TYMPX"
        status = "PASS" if min(len(present), len(f1_excluded), len(tymp_excluded)) >= 5 else "FAIL_LOW_COVERAGE"
        panel_rows.append({
            "stable_id": stable_id, "pathway": name, "category": category,
            "genes_total": len(genes), "genes_present": len(present),
            "genes_after_F1_exclusion": len(f1_excluded),
            "genes_after_TYMP_exclusion": len(tymp_excluded),
            "coverage": len(present) / len(genes) if genes else 0,
            "score_column": score_column if status == "PASS" else "",
            "F1_excluded_score_column": f1_excluded_column if status == "PASS" else "",
            "TYMP_excluded_score_column": tymp_excluded_column if status == "PASS" else "",
            "status": status,
        })
        if status != "PASS":
            continue
        for variant_column, variant_genes in [
            (score_column, present),
            (f1_excluded_column, f1_excluded),
            (tymp_excluded_column, tymp_excluded),
        ]:
            sc.tl.score_genes(
                merged, gene_list=variant_genes,
                ctrl_size=max(20, min(50, len(variant_genes))),
                n_bins=25, score_name=variant_column, random_state=SEED,
                ctrl_as_ref=False, use_raw=False,
            )
            score_to_pathway[variant_column] = {
                "stable_id": stable_id, "pathway": name, "category": category,
            }
        score_columns.append(score_column)
        f1_excluded_score_columns.append(f1_excluded_column)
        tymp_excluded_score_columns.append(tymp_excluded_column)
        print(f"scored {stable_id}: {name} ({len(present)}/{len(genes)} genes)", flush=True)

    panel = pd.DataFrame(panel_rows)
    write_tsv(panel, args.output_dir / "gse207493_metabolic_pathway_panel_v1.tsv")
    if len(score_columns) < 20:
        raise RuntimeError(f"only {len(score_columns)} metabolic pathways passed coverage")

    celltype_rows = []
    heterogeneity_rows = []
    for (sample, cell_type), group in merged.obs.groupby(
        ["sample", "cell_type"], observed=True
    ):
        if len(group) < 20:
            continue
        for score_column in score_columns:
            values = group[score_column].to_numpy(float)
            info = score_to_pathway[score_column]
            row = {
                "sample": sample, "cell_type": cell_type,
                "stable_id": info["stable_id"], "pathway": info["pathway"],
                "category": info["category"], "n_cells": len(values),
                "mean_score": np.mean(values), "median_score": np.median(values),
                "sd_score": np.std(values, ddof=1),
                "iqr_score": np.quantile(values, 0.75) - np.quantile(values, 0.25),
                "mad_score": stats.median_abs_deviation(values, scale=1.0),
            }
            heterogeneity_rows.append(row)
    heterogeneity = pd.DataFrame(heterogeneity_rows)
    write_tsv(
        heterogeneity,
        args.output_dir / "gse207493_metabolic_heterogeneity_by_sample_celltype_v1.tsv.gz",
    )
    heterogeneity_summary = heterogeneity.groupby(
        ["cell_type", "stable_id", "pathway", "category"], observed=True
    ).agg(
        n_samples=("sample", "nunique"), median_cells=("n_cells", "median"),
        median_mean_score=("mean_score", "median"),
        median_iqr_score=("iqr_score", "median"),
        median_mad_score=("mad_score", "median"),
    ).reset_index()
    write_tsv(
        heterogeneity_summary,
        args.output_dir / "gse207493_metabolic_heterogeneity_summary_v1.tsv",
    )

    contrasts = {
        "F1_high_vs_low_Tumor": ("Tumor", "f1_score"),
        "TYMP_high_vs_low_MyeloidDC": ("Myeloid_DC", "TYMP_log1p"),
        "TYMP_high_vs_low_Treg": ("Treg", "TYMP_log1p"),
    }
    contrast_score_columns = {
        "F1_high_vs_low_Tumor": f1_excluded_score_columns,
        "TYMP_high_vs_low_MyeloidDC": tymp_excluded_score_columns,
        "TYMP_high_vs_low_Treg": tymp_excluded_score_columns,
    }
    effect_rows = []
    correlation_rows = []
    for sample in EXPECTED_SAMPLES:
        sample_frame = merged.obs.loc[merged.obs["sample"] == sample].copy()
        for contrast, (lineage, feature) in contrasts.items():
            active_score_columns = contrast_score_columns[contrast]
            if lineage == "Myeloid_DC":
                subset = sample_frame.loc[
                    sample_frame.cell_type.isin(["Myeloid", "Dendritic"])
                ].copy()
            else:
                subset = sample_frame.loc[sample_frame.cell_type == lineage].copy()
            if len(subset) < 20:
                continue
            pathway_values = subset[active_score_columns].to_numpy(float)
            predictor = subset[feature].to_numpy(float)
            s_score = subset["S_score"].to_numpy(float)
            g2m_score = subset["G2M_score"].to_numpy(float)
            pathway_residual = residualize_matrix(pathway_values, s_score, g2m_score)
            predictor_residual = residualize_matrix(
                predictor.reshape(-1, 1), s_score, g2m_score
            ).ravel()
            raw_masks = quantile_masks(predictor)
            adjusted_masks = quantile_masks(predictor_residual)
            for column_index, score_column in enumerate(active_score_columns):
                info = score_to_pathway[score_column]
                raw_rho, _ = stats.spearmanr(predictor, pathway_values[:, column_index])
                adjusted_rho, _ = stats.spearmanr(
                    predictor_residual, pathway_residual[:, column_index]
                )
                correlation_rows.append({
                    "sample": sample, "contrast": contrast, "lineage": lineage,
                    "feature": feature, "stable_id": info["stable_id"],
                    "pathway": info["pathway"], "category": info["category"],
                    "n_cells": len(subset), "raw_spearman_rho": raw_rho,
                    "cycle_adjusted_spearman_rho": adjusted_rho,
                })
                if raw_masks is None or adjusted_masks is None:
                    continue
                raw_high, raw_low = raw_masks
                adj_high, adj_low = adjusted_masks
                effect_rows.append({
                    "sample": sample, "contrast": contrast, "lineage": lineage,
                    "feature": feature, "stable_id": info["stable_id"],
                    "pathway": info["pathway"], "category": info["category"],
                    "n_cells": len(subset),
                    "raw_high_n": int(raw_high.sum()), "raw_low_n": int(raw_low.sum()),
                    "raw_mean_difference": float(
                        pathway_values[raw_high, column_index].mean()
                        - pathway_values[raw_low, column_index].mean()
                    ),
                    "adjusted_high_n": int(adj_high.sum()), "adjusted_low_n": int(adj_low.sum()),
                    "cycle_adjusted_mean_difference": float(
                        pathway_residual[adj_high, column_index].mean()
                        - pathway_residual[adj_low, column_index].mean()
                    ),
                })

    effects = pd.DataFrame(effect_rows)
    correlations = pd.DataFrame(correlation_rows)
    write_tsv(effects, args.output_dir / "gse207493_metabolic_contrasts_by_sample_v1.tsv.gz")
    write_tsv(correlations, args.output_dir / "gse207493_metabolic_correlations_by_sample_v1.tsv.gz")

    effect_summary_rows = []
    for keys, group in effects.groupby(
        ["contrast", "lineage", "feature", "stable_id", "pathway", "category"],
        observed=True,
    ):
        raw_stat, raw_p = signed_rank(group["raw_mean_difference"])
        adj_stat, adj_p = signed_rank(group["cycle_adjusted_mean_difference"])
        effect_summary_rows.append({
            "contrast": keys[0], "lineage": keys[1], "feature": keys[2],
            "stable_id": keys[3], "pathway": keys[4], "category": keys[5],
            "n_samples": len(group),
            "median_raw_difference": group["raw_mean_difference"].median(),
            "raw_wilcoxon_statistic": raw_stat, "raw_p_value": raw_p,
            "median_cycle_adjusted_difference": group[
                "cycle_adjusted_mean_difference"
            ].median(),
            "adjusted_wilcoxon_statistic": adj_stat, "adjusted_p_value": adj_p,
        })
    effect_summary = pd.DataFrame(effect_summary_rows)
    effect_summary["raw_q_value"] = effect_summary.groupby("contrast", observed=True)[
        "raw_p_value"
    ].transform(bh_adjust)
    effect_summary["adjusted_q_value"] = effect_summary.groupby("contrast", observed=True)[
        "adjusted_p_value"
    ].transform(bh_adjust)
    effect_summary["direction"] = np.where(
        effect_summary["median_cycle_adjusted_difference"] > 0,
        "higher_in_feature_high", "lower_in_feature_high",
    )
    effect_summary["retained"] = (
        (effect_summary["adjusted_q_value"] < 0.05)
        | (effect_summary["adjusted_p_value"] < 0.05)
    )
    write_tsv(effect_summary, args.output_dir / "gse207493_metabolic_contrast_summary_full_v1.tsv")
    write_tsv(
        effect_summary.loc[effect_summary.retained].sort_values(
            ["contrast", "adjusted_q_value", "adjusted_p_value"]
        ),
        args.output_dir / "gse207493_metabolic_contrast_summary_retained_v1.tsv",
    )

    correlation_summary_rows = []
    for keys, group in correlations.groupby(
        ["contrast", "lineage", "feature", "stable_id", "pathway", "category"],
        observed=True,
    ):
        raw_stat, raw_p = signed_rank(group["raw_spearman_rho"])
        adj_stat, adj_p = signed_rank(group["cycle_adjusted_spearman_rho"])
        correlation_summary_rows.append({
            "contrast": keys[0], "lineage": keys[1], "feature": keys[2],
            "stable_id": keys[3], "pathway": keys[4], "category": keys[5],
            "n_samples": len(group),
            "median_raw_spearman_rho": group["raw_spearman_rho"].median(),
            "raw_wilcoxon_statistic": raw_stat, "raw_p_value": raw_p,
            "median_cycle_adjusted_spearman_rho": group[
                "cycle_adjusted_spearman_rho"
            ].median(),
            "adjusted_wilcoxon_statistic": adj_stat, "adjusted_p_value": adj_p,
        })
    correlation_summary = pd.DataFrame(correlation_summary_rows)
    correlation_summary["raw_q_value"] = correlation_summary.groupby(
        "contrast", observed=True
    )["raw_p_value"].transform(bh_adjust)
    correlation_summary["adjusted_q_value"] = correlation_summary.groupby(
        "contrast", observed=True
    )["adjusted_p_value"].transform(bh_adjust)
    correlation_summary["direction"] = np.where(
        correlation_summary["median_cycle_adjusted_spearman_rho"] > 0,
        "positive", "negative",
    )
    correlation_summary["retained"] = (
        (correlation_summary["adjusted_q_value"] < 0.05)
        | (correlation_summary["adjusted_p_value"] < 0.05)
    )
    write_tsv(
        correlation_summary,
        args.output_dir / "gse207493_metabolic_correlation_summary_full_v1.tsv",
    )
    write_tsv(
        correlation_summary.loc[correlation_summary.retained].sort_values(
            ["contrast", "adjusted_q_value", "adjusted_p_value"]
        ),
        args.output_dir / "gse207493_metabolic_correlation_summary_retained_v1.tsv",
    )

    cell_output = merged.obs.reset_index().rename(columns={"index": "cell_barcode"})[
        ["cell_barcode", "sample", "cell_type", "phase", "S_score", "G2M_score",
         "f1_score", "TYMP_log1p"] + score_columns
    ]
    write_tsv(cell_output, args.output_dir / "gse207493_metabolic_scores_by_cell_v1.tsv.gz")

    gates = pd.DataFrame([
        {"gate": "qc_cells", "observed": merged.n_obs, "target": 102723,
         "status": "PASS" if merged.n_obs == 102723 else "FAIL",
         "note": "exact author-reported QC total"},
        {"gate": "cell_metadata_matched", "observed": merged.n_obs - len(missing_metadata),
         "target": merged.n_obs, "status": "PASS" if not len(missing_metadata) else "FAIL",
         "note": "frozen cell type/F1/TYMP/cycle metadata"},
        {"gate": "pathways_requested", "observed": len(PATHWAY_PANEL),
         "target": len(PATHWAY_PANEL), "status": "PASS", "note": "frozen targeted panel"},
        {"gate": "pathways_scored", "observed": len(score_columns),
         "target": len(PATHWAY_PANEL),
         "status": "PASS" if len(score_columns) >= 20 else "FAIL",
         "note": "minimum five detected genes per pathway"},
        {"gate": "anchor_excluded_pathway_scores", "observed": len(
            f1_excluded_score_columns
        ) + len(tymp_excluded_score_columns), "target": 2 * len(score_columns),
         "status": "PASS" if (
             len(f1_excluded_score_columns) == len(score_columns)
             and len(tymp_excluded_score_columns) == len(score_columns)
         ) else "FAIL",
         "note": "F1 genes removed for F1 contrasts; TYMP removed for TYMP contrasts"},
        {"gate": "cell_scores_complete", "observed": int(
            merged.obs[score_columns].notna().all(axis=1).sum()
        ), "target": merged.n_obs,
         "status": "PASS" if merged.obs[score_columns].notna().all().all() else "FAIL",
         "note": "all retained pathway scores finite"},
    ])
    write_tsv(gates, args.output_dir / "gse207493_metabolic_data_gate_v1.tsv")
    print(gates.to_string(index=False), flush=True)
    for contrast in contrasts:
        top = effect_summary.loc[
            effect_summary.contrast == contrast
        ].sort_values(["adjusted_q_value", "adjusted_p_value"]).head(10)
        print(f"\nTOP {contrast}\n{top.to_string(index=False)}", flush=True)


if __name__ == "__main__":
    main()
