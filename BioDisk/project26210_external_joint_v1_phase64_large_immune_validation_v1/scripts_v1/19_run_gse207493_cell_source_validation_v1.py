#!/usr/bin/env python3
"""Reproduce GSE207493 QC and resolve F1/TYMP expression by cell lineage.

The workflow follows the authors' 19-sample QC thresholds and uses Harmony plus
Leiden for unsupervised clusters. Cluster names are assigned from explicit marker
panels because numeric Leiden labels are not portable across software versions.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
from pathlib import Path

for variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(variable, "8")

import anndata as ad
import harmonypy
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

BROAD_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "TRAC", "CD247"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1"],
    "Myeloid": ["LYZ", "LST1", "FCER1G", "TYROBP", "CTSS"],
    "Dendritic": ["CD1C", "FCER1A", "CLEC10A", "HLA-DQA1"],
    "B_cell": ["CD79A", "MS4A1", "CD37", "CD74", "CD22"],
    "Plasma": ["MZB1", "JCHAIN", "SDC1", "IGHG1"],
    "Mast": ["KIT", "TPSAB1", "TPSB2", "CPA3"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "CDH5"],
    "CAF": ["COL1A1", "COL1A2", "DCN", "COL3A1", "PDGFRA"],
    "Tumor": ["CA9", "NDUFA4L2", "EPCAM", "KRT8", "KRT18", "KRT19"],
}
T_MARKERS = {
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"],
    "Exhausted_T": ["PDCD1", "HAVCR2", "LAG3", "TOX", "CXCL13"],
    "Proliferating_T": ["MKI67", "TOP2A", "TYMS", "STMN1"],
    "CD8_T": ["CD8A", "CD8B", "CST7", "CCL5"],
    "CD4_T": ["CD4", "IL7R", "LTB", "MAL"],
}


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def find_matrices(root: Path) -> dict[str, tuple[Path, str]]:
    candidates: dict[str, list[tuple[Path, str]]] = {}
    for matrix in root.rglob("*matrix.mtx*"):
        if not matrix.is_file():
            continue
        prefix = matrix.name.split("matrix.mtx", 1)[0]
        match = re.search(r"RCC\d+", str(matrix), flags=re.IGNORECASE)
        if match is None:
            continue
        sample = match.group(0).upper()
        candidates.setdefault(sample, []).append((matrix.parent, prefix))

    chosen: dict[str, tuple[Path, str]] = {}
    for sample, entries in candidates.items():
        entries.sort(key=lambda item: (
            "filtered_feature_bc_matrix" not in str(item[0]),
            len(str(item[0])),
        ))
        chosen[sample] = entries[0]
    return chosen


def vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=np.float32)
    value = adata[:, gene].X
    if sparse.issparse(value):
        return value.toarray().ravel().astype(np.float32)
    return np.asarray(value).ravel().astype(np.float32)


def mean_score(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=np.float32)
    value = adata[:, present].X.mean(axis=1)
    return np.asarray(value).ravel().astype(np.float32)


def weighted_score(adata: ad.AnnData, genes: list[str], weights: np.ndarray) -> np.ndarray:
    index = [adata.var_names.get_loc(gene) for gene in genes]
    value = adata.X[:, index]
    score = value @ weights
    return np.asarray(score).ravel().astype(np.float32)


def z_columns(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    for column in output:
        sd = output[column].std(ddof=0)
        output[column] = 0.0 if not np.isfinite(sd) or sd == 0 else (
            output[column] - output[column].mean()
        ) / sd
    return output


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


def assign_clusters(adata: ad.AnnData) -> tuple[pd.DataFrame, dict[str, str]]:
    cell_scores: dict[str, np.ndarray] = {}
    for label, genes in {**BROAD_MARKERS, **T_MARKERS}.items():
        cell_scores[label] = mean_score(adata, genes)

    score_frame = pd.DataFrame(cell_scores, index=adata.obs_names)
    score_frame["cluster"] = adata.obs["leiden"].astype(str).values
    cluster_scores = score_frame.groupby("cluster", observed=True).mean()
    broad_z = z_columns(cluster_scores[list(BROAD_MARKERS)])
    mapping = broad_z.idxmax(axis=1).to_dict()

    t_clusters = [cluster for cluster, label in mapping.items() if label == "T_cell"]
    if t_clusters:
        t_z = z_columns(cluster_scores.loc[t_clusters, list(T_MARKERS)])
        for cluster in t_clusters:
            mapping[cluster] = str(t_z.loc[cluster].idxmax())

    cluster_table = cluster_scores.reset_index()
    cluster_table["cell_type"] = cluster_table["cluster"].map(mapping)
    counts = adata.obs["leiden"].astype(str).value_counts()
    cluster_table["n_cells"] = cluster_table["cluster"].map(counts).astype(int)
    cluster_table["marker_method"] = "maximum across-cluster z-scored mean marker expression"
    return cluster_table, mapping


def paired_tests(sample_celltype: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for feature in ["f1_score", "TYMP_log1p"]:
        for reference, members in [
            ("other_T", ["CD4_T", "CD8_T", "Exhausted_T", "Proliferating_T"]),
            ("all_non_Treg", None),
        ]:
            left = sample_celltype[sample_celltype.cell_type == "Treg"]
            if members is None:
                right = sample_celltype[sample_celltype.cell_type != "Treg"]
            else:
                right = sample_celltype[sample_celltype.cell_type.isin(members)]
            right = right.assign(weighted_value=right[feature] * right["n_cells"])
            right = right.groupby("sample", as_index=False, observed=True).agg(
                weighted_value=("weighted_value", "sum"),
                total_cells=("n_cells", "sum"),
            )
            right[feature] = right["weighted_value"] / right["total_cells"]
            merged = left[["sample", feature]].merge(
                right[["sample", feature]], on="sample", suffixes=("_treg", "_reference")
            ).dropna()
            if len(merged) >= 5 and np.any(merged[f"{feature}_treg"] != merged[f"{feature}_reference"]):
                test = stats.wilcoxon(
                    merged[f"{feature}_treg"], merged[f"{feature}_reference"],
                    alternative="two-sided",
                )
                statistic, p_value = float(test.statistic), float(test.pvalue)
            else:
                statistic, p_value = np.nan, np.nan
            rows.append({
                "comparison": f"Treg_vs_{reference}", "feature": feature,
                "n_paired_samples": len(merged),
                "median_paired_difference": (
                    merged[f"{feature}_treg"] - merged[f"{feature}_reference"]
                ).median() if len(merged) else np.nan,
                "wilcoxon_statistic": statistic, "p_value": p_value,
            })
    result = pd.DataFrame(rows)
    result["q_value"] = bh_adjust(result["p_value"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_root", type=Path)
    parser.add_argument("signature_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 8
    np.random.seed(SEED)

    matrices = find_matrices(args.raw_root)
    missing_samples = sorted(set(EXPECTED_SAMPLES) - set(matrices))
    extra_samples = sorted(set(matrices) - set(EXPECTED_SAMPLES))
    if missing_samples:
        raise RuntimeError(f"missing expected sample matrices: {missing_samples}")

    qc_rows: list[dict[str, object]] = []
    samples: list[ad.AnnData] = []
    for sample in EXPECTED_SAMPLES:
        directory, prefix = matrices[sample]
        one = sc.read_10x_mtx(directory, var_names="gene_symbols", make_unique=True, prefix=prefix)
        raw_cells, raw_genes = one.n_obs, one.n_vars
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
        qc_rows.append({
            "sample": sample, "matrix_dir": str(directory), "prefix": prefix,
            "raw_cells": raw_cells, "raw_genes": raw_genes,
            "min_cells_per_gene": MIN_CELLS[sample],
            "max_features": MAX_FEATURES[sample], "qc_cells": one.n_obs,
            "qc_genes": one.n_vars,
        })
        samples.append(one)

    merged = ad.concat(samples, join="outer", merge="same", fill_value=0)
    merged.obs["sample"] = pd.Categorical(merged.obs["sample"], categories=EXPECTED_SAMPLES)
    merged.layers["counts"] = merged.X.copy()
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)

    signature = pd.read_csv(args.signature_tsv, sep="\t")
    signature = signature[signature["signature_id"] == "f1_top100_weighted"].copy()
    signature["gene"] = signature["gene"].astype(str).str.upper()
    signature["raw_weight"] = pd.to_numeric(signature["raw_weight"], errors="coerce")
    present = signature[signature["gene"].isin(merged.var_names)].copy()
    if len(present) < 70:
        raise RuntimeError(f"only {len(present)} F1 genes present")
    weights = present["raw_weight"].to_numpy(float)
    weights = weights / weights.sum()
    merged.obs["f1_score"] = weighted_score(merged, present["gene"].tolist(), weights)
    merged.obs["TYMP_log1p"] = vector(merged, "TYMP")

    sc.pp.highly_variable_genes(merged, flavor="seurat", n_top_genes=2000, subset=False)
    cluster_data = merged[:, merged.var["highly_variable"]].copy()
    sc.pp.scale(cluster_data, max_value=10)
    sc.tl.pca(cluster_data, n_comps=30, svd_solver="arpack", random_state=SEED)
    harmony = harmonypy.run_harmony(
        cluster_data.obsm["X_pca"], cluster_data.obs, "sample",
        max_iter_harmony=20, random_state=SEED,
    )
    harmony_embedding = np.asarray(harmony.Z_corr)
    if harmony_embedding.shape == (cluster_data.n_obs, 30):
        pass
    elif harmony_embedding.shape == (30, cluster_data.n_obs):
        harmony_embedding = harmony_embedding.T
    else:
        raise RuntimeError(
            f"unexpected Harmony matrix shape {harmony_embedding.shape}; "
            f"expected {(cluster_data.n_obs, 30)}"
        )
    cluster_data.obsm["X_pca_harmony"] = harmony_embedding
    sc.pp.neighbors(cluster_data, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.leiden(cluster_data, resolution=0.6, random_state=SEED, key_added="leiden")
    merged.obs["leiden"] = cluster_data.obs["leiden"].astype(str).values

    cluster_table, mapping = assign_clusters(merged)
    merged.obs["cell_type"] = merged.obs["leiden"].astype(str).map(mapping)

    sample_qc = pd.DataFrame(qc_rows)
    write_tsv(sample_qc, args.output_dir / "gse207493_sample_qc_v1.tsv")
    write_tsv(cluster_table, args.output_dir / "gse207493_cluster_annotation_v1.tsv")

    summary = merged.obs.groupby("cell_type", observed=True).agg(
        n_cells=("cell_type", "size"), n_samples=("sample", "nunique"),
        f1_mean=("f1_score", "mean"), f1_median=("f1_score", "median"),
        TYMP_mean=("TYMP_log1p", "mean"), TYMP_median=("TYMP_log1p", "median"),
        TYMP_positive_fraction=("TYMP_log1p", lambda x: float(np.mean(x > 0))),
    ).reset_index()
    write_tsv(summary, args.output_dir / "gse207493_celltype_feature_summary_v1.tsv")

    sample_celltype = merged.obs.groupby(["sample", "cell_type"], observed=True).agg(
        n_cells=("cell_type", "size"), f1_score=("f1_score", "mean"),
        TYMP_log1p=("TYMP_log1p", "mean"),
        TYMP_positive_fraction=("TYMP_log1p", lambda x: float(np.mean(x > 0))),
    ).reset_index()
    sample_celltype = sample_celltype[sample_celltype["n_cells"] >= 20].copy()
    write_tsv(sample_celltype, args.output_dir / "gse207493_sample_celltype_pseudobulk_v1.tsv")
    write_tsv(paired_tests(sample_celltype), args.output_dir / "gse207493_treg_paired_tests_v1.tsv")

    gene_rows: list[dict[str, object]] = []
    genes_to_report = present["gene"].tolist() + (["TYMP"] if "TYMP" in merged.var_names else [])
    for gene in genes_to_report:
        values = vector(merged, gene)
        temp = pd.DataFrame({
            "cell_type": merged.obs["cell_type"].values,
            "value": values,
        })
        for cell_type, group in temp.groupby("cell_type", observed=True):
            gene_rows.append({
                "gene": gene, "cell_type": cell_type, "n_cells": len(group),
                "mean_log1p": group["value"].mean(),
                "positive_fraction": np.mean(group["value"] > 0),
            })
    write_tsv(pd.DataFrame(gene_rows), args.output_dir / "gse207493_f1_tymp_by_celltype_v1.tsv.gz")

    cell_meta = merged.obs.reset_index().rename(columns={"index": "cell_barcode"})[
        ["cell_barcode", "sample", "leiden", "cell_type", "f1_score", "TYMP_log1p"]
    ]
    write_tsv(cell_meta, args.output_dir / "gse207493_cell_metadata_scores_v1.tsv.gz")

    gates = pd.DataFrame([
        {"gate": "expected_samples", "observed": len(matrices), "target": 19,
         "status": "PASS" if len(matrices) == 19 and not extra_samples else "FAIL",
         "note": f"missing={missing_samples}; extra={extra_samples}"},
        {"gate": "qc_cells", "observed": merged.n_obs, "target": 102723,
         "status": "PASS" if abs(merged.n_obs - 102723) / 102723 <= 0.02 else "REVIEW",
         "note": "authors reported 102,723 high-quality cells; 2% reproducibility tolerance"},
        {"gate": "leiden_clusters", "observed": merged.obs["leiden"].nunique(), "target": 22,
         "status": "PASS" if 15 <= merged.obs["leiden"].nunique() <= 30 else "REVIEW",
         "note": "authors reported 22 clusters at resolution 0.6"},
        {"gate": "f1_genes_present", "observed": len(present), "target": 100,
         "status": "PASS" if len(present) >= 90 else "FAIL", "note": "frozen weighted F1"},
        {"gate": "TYMP_present", "observed": int("TYMP" in merged.var_names), "target": 1,
         "status": "PASS" if "TYMP" in merged.var_names else "FAIL", "note": "direct gene"},
    ])
    write_tsv(gates, args.output_dir / "gse207493_data_gate_v1.tsv")
    print(gates.to_string(index=False))
    print(summary.sort_values("n_cells", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
