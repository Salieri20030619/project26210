#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GLMF1、GLMF2 空间生态位独立重算。
旧 obs["niche"] 不参与；两个样本分别构图、选 K、聚类、编号和出图。
"""

from pathlib import Path
import json

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import squidpy as sq


# --------------------------- 参数区 ---------------------------
BASE = Path("/mnt/BioDisk/ludingyu/project26210/results/part4_spatial")
INPUTS = {
    "GLMF1": BASE / "analysis/GLMF1_analyzed.h5ad",
    "GLMF2": BASE / "analysis/GLMF2_analyzed.h5ad",
}
OUT_ROOT = BASE / "niche_separate_v1"

ABUNDANCE_KEYS = [
    "cell_abundance_clean",
    "q05_cell_abundance_w_sf",
    "means_cell_abundance_w_sf",
]
K_RANGE = range(3, 9)
FIXED_K = {"GLMF1": None, "GLMF2": None}  # 如需固定可改为 5 等整数
N_PERM = 500
SEED = 20260723
SAVE_NEW_H5AD = True

COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
]


def read_abundance(adata):
    key = next((k for k in ABUNDANCE_KEYS if k in adata.obsm), None)
    if key is None:
        raise KeyError(f"未找到丰度矩阵；现有 obsm: {list(adata.obsm.keys())}")

    obj = adata.obsm[key]
    if isinstance(obj, pd.DataFrame):
        x = obj.to_numpy(dtype=float)
        celltypes = [str(c) for c in obj.columns]
    else:
        x = np.asarray(obj, dtype=float)
        mod = adata.uns.get("mod", {})
        saved_names = adata.uns.get(f"{key}_columns", [])
        if len(saved_names) == x.shape[1]:
            celltypes = [str(c) for c in saved_names]
        else:
            celltypes = (
                [str(c) for c in mod.get("factor_names", [])]
                if hasattr(mod, "get") else []
            )

    if x.shape != (adata.n_obs, len(celltypes)):
        raise ValueError(
            f"{key} shape={x.shape}，cell type names={len(celltypes)}；"
            "请检查 adata.uns['mod']['factor_names']"
        )
    x[~np.isfinite(x)] = 0
    x = np.clip(x, 0, None)
    if np.any(x.sum(axis=1) == 0):
        raise ValueError("存在总细胞丰度为 0 的 spot")
    return x, celltypes, key


def build_features(adata, abundance):
    # 每个样本单独重建 Visium 一阶六邻域，不复用旧图。
    sq.gr.spatial_neighbors(
        adata, coord_type="grid", n_neighs=6, n_rings=1
    )
    graph = adata.obsp["spatial_connectivities"].copy().tocsr()
    graph.data[:] = 1
    graph = graph.maximum(graph.T)

    degree = np.asarray(graph.sum(axis=1)).ravel()
    if np.median(degree) < 3:
        raise ValueError("空间邻接图异常，请检查 array_row/array_col")

    # spot 本身 + 一阶邻域的细胞丰度，转为组成后做 Hellinger 变换。
    graph_self = graph + sparse.eye(adata.n_obs, format="csr")
    neighborhood_abundance = np.asarray(graph_self @ abundance)
    composition = neighborhood_abundance / neighborhood_abundance.sum(
        axis=1, keepdims=True
    )
    features = np.sqrt(composition)
    return features, composition, degree


def choose_k(features, fixed_k, seed):
    min_required = max(25, int(np.ceil(0.01 * len(features))))
    rows, labels_by_k = [], {}

    for k in K_RANGE:
        model = KMeans(
            n_clusters=k, n_init=50, max_iter=1000, random_state=seed
        )
        labels = model.fit_predict(features)
        labels_by_k[k] = labels

        if len(features) > 2500:
            rng = np.random.default_rng(seed + k)
            idx = np.sort(rng.choice(len(features), 2500, replace=False))
            sil = silhouette_score(features[idx], labels[idx])
        else:
            sil = silhouette_score(features, labels)

        min_size = int(np.bincount(labels, minlength=k).min())
        rows.append({
            "k": k,
            "silhouette": sil,
            "minimum_niche_spots": min_size,
            "minimum_required": min_required,
            "valid": min_size >= min_required,
        })

    table = pd.DataFrame(rows)
    if fixed_k is None:
        valid = table[table["valid"]]
        if valid.empty:
            raise ValueError("K=3~8 均未通过最小生态位大小限制")
        selected = int(valid.sort_values("silhouette").iloc[-1]["k"])
    else:
        selected = int(fixed_k)
        if selected not in labels_by_k:
            raise ValueError(f"FIXED_K 必须位于 {list(K_RANGE)}")

    table["selected"] = table["k"].eq(selected)
    return selected, labels_by_k[selected], table


def relabel(raw_labels, sample):
    # 本样本内按 spot 数从多到少编号，跨样本编号不作对应解释。
    old, counts = np.unique(raw_labels, return_counts=True)
    ordered = old[np.argsort(-counts)]
    mapping = {old_id: new_id for new_id, old_id in enumerate(ordered)}
    labels = np.array([mapping[x] for x in raw_labels], dtype=int)
    names = [f"{sample}-N{i + 1}" for i in range(len(ordered))]
    return labels, names


def draw_heatmap(ax, matrix, rows, cols, title, annotate=False):
    lim = max(1.0, float(np.nanpercentile(np.abs(matrix), 98)))
    im = ax.imshow(
        matrix, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto"
    )
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=8)
    ax.set_title(title, loc="left", fontweight="bold")
    if annotate:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                ax.text(j, i, f"{matrix[i, j]:.1f}",
                        ha="center", va="center", fontsize=7)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def make_figure(adata, sample, labels, names, sizes, enrichment, zscore, path):
    k = len(names)
    colors = COLORS[:k]
    fig, axes = plt.subplots(2, 2, figsize=(16, 13), constrained_layout=True)

    spatial = adata.uns["spatial"]
    lib = sample if sample in spatial else next(iter(spatial))
    info = spatial[lib]
    image = np.asarray(info["images"]["hires"])
    scale = float(info["scalefactors"]["tissue_hires_scalef"])
    diameter = float(info["scalefactors"]["spot_diameter_fullres"]) * scale
    xy = np.asarray(adata.obsm["spatial"]) * scale

    ax = axes[0, 0]
    ax.imshow(image)
    ax.scatter(
        xy[:, 0], xy[:, 1], c=labels,
        cmap=ListedColormap(colors), vmin=-0.5, vmax=k - 0.5,
        s=max(5, diameter ** 2 * 0.18), linewidths=0, alpha=0.82,
    )
    handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="",
            markerfacecolor=colors[i], markeredgecolor="none", label=names[i]
        ) for i in range(k)
    ]
    ax.legend(
        handles=handles, bbox_to_anchor=(1.01, 1),
        loc="upper left", frameon=False
    )
    ax.set_title(f"A  {sample} spatial niches", loc="left", fontweight="bold")
    ax.axis("off")

    ax = axes[0, 1]
    pct = sizes["fraction"].to_numpy() * 100
    ax.barh(np.arange(k), pct[::-1], color=colors[::-1])
    ax.set_yticks(np.arange(k))
    ax.set_yticklabels(names[::-1])
    ax.set_xlabel("Spots (%)")
    ax.set_title("B  Niche size", loc="left", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    for y, value in enumerate(pct[::-1]):
        ax.text(value, y, f" {value:.1f}%", va="center", fontsize=8)

    values = enrichment.to_numpy()
    keep = np.argsort(-np.max(np.abs(values), axis=0))[:18]
    draw_heatmap(
        axes[1, 0], values[:, keep], names,
        enrichment.columns[keep].tolist(),
        "C  Cell-type log2 enrichment"
    )
    draw_heatmap(
        axes[1, 1], zscore, names, names,
        f"D  Niche adjacency ({N_PERM} permutations; z-score)", True
    )

    fig.suptitle(
        f"{sample}: independently computed spatial niches",
        fontsize=17, fontweight="bold"
    )
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def process_sample(sample, path, seed):
    print(f"\n[{sample}] 读取 {path}", flush=True)
    adata = ad.read_h5ad(path)
    out = OUT_ROOT / sample
    out.mkdir(parents=True, exist_ok=True)

    abundance, celltypes, abundance_key = read_abundance(adata)
    features, composition, degree = build_features(adata, abundance)
    selected_k, raw_labels, k_table = choose_k(
        features, FIXED_K[sample], seed
    )
    labels, names = relabel(raw_labels, sample)
    k = len(names)

    counts = np.bincount(labels, minlength=k)
    sizes = pd.DataFrame({
        "niche": names,
        "n_spots": counts,
        "fraction": counts / counts.sum(),
    })
    mean_abundance = np.vstack([
        abundance[labels == i].mean(axis=0) for i in range(k)
    ])
    mean_composition = np.vstack([
        composition[labels == i].mean(axis=0) for i in range(k)
    ])
    global_composition = composition.mean(axis=0)
    log2_enrichment = np.log2(
        (mean_composition + 1e-8) / (global_composition[None, :] + 1e-8)
    )

    abundance_df = pd.DataFrame(
        mean_abundance, index=names, columns=celltypes
    )
    composition_df = pd.DataFrame(
        mean_composition, index=names, columns=celltypes
    )
    enrichment_df = pd.DataFrame(
        log2_enrichment, index=names, columns=celltypes
    )

    suggestions = []
    for i, niche in enumerate(names):
        top = np.argsort(-log2_enrichment[i])[:3]
        suggestions.append({
            "niche": niche,
            "suggested_name": (
                f"{celltypes[top[0]]} / {celltypes[top[1]]}-enriched"
            ),
            "top1": celltypes[top[0]],
            "top1_log2_enrichment": log2_enrichment[i, top[0]],
            "top2": celltypes[top[1]],
            "top2_log2_enrichment": log2_enrichment[i, top[1]],
            "top3": celltypes[top[2]],
            "top3_log2_enrichment": log2_enrichment[i, top[2]],
        })
    suggestions = pd.DataFrame(suggestions)

    adata.obs["niche_independent"] = pd.Categorical(
        [names[i] for i in labels], categories=names, ordered=True
    )
    sq.gr.nhood_enrichment(
        adata,
        cluster_key="niche_independent",
        n_perms=N_PERM,
        seed=seed + 1000,
    )
    nhood = adata.uns["niche_independent_nhood_enrichment"]
    zscore = np.asarray(nhood["zscore"])
    observed = np.asarray(nhood["count"])

    spot_df = adata.obs[["array_row", "array_col"]].copy()
    spot_df.insert(0, "barcode", adata.obs_names.astype(str))
    spot_df["niche_independent"] = adata.obs["niche_independent"].astype(str)

    spot_df.to_csv(out / f"{sample}_spot_niche_assignments.csv", index=False)
    k_table.to_csv(out / f"{sample}_k_selection.csv", index=False)
    sizes.to_csv(out / f"{sample}_niche_sizes.csv", index=False)
    abundance_df.to_csv(out / f"{sample}_mean_cell_abundance.csv")
    composition_df.to_csv(out / f"{sample}_mean_neighborhood_composition.csv")
    enrichment_df.to_csv(out / f"{sample}_celltype_log2_enrichment.csv")
    suggestions.to_csv(out / f"{sample}_niche_name_suggestions.csv", index=False)
    pd.DataFrame(observed, index=names, columns=names).to_csv(
        out / f"{sample}_niche_adjacency_count.csv"
    )
    pd.DataFrame(zscore, index=names, columns=names).to_csv(
        out / f"{sample}_niche_adjacency_zscore.csv"
    )

    make_figure(
        adata, sample, labels, names, sizes, enrichment_df, zscore,
        out / f"{sample}_spatial_niche_independent.png"
    )

    method = {
        "sample": sample,
        "input": str(path),
        "abundance_key": abundance_key,
        "n_spots": int(adata.n_obs),
        "n_celltypes": int(len(celltypes)),
        "median_first_ring_degree": float(np.median(degree)),
        "selected_k": int(selected_k),
        "k_selection_mode": (
            "auto_silhouette" if FIXED_K[sample] is None else "fixed"
        ),
        "feature": (
            "Hellinger-transformed cell-type composition of each spot "
            "plus first-ring Visium neighbors"
        ),
        "clustering": "KMeans independently fitted within sample",
        "old_obs_niche_used": False,
        "seed": int(seed),
    }
    with open(out / f"{sample}_method.json", "w", encoding="utf-8") as f:
        json.dump(method, f, ensure_ascii=False, indent=2)

    if SAVE_NEW_H5AD:
        adata.uns["niche_independent_method"] = method
        adata.write_h5ad(
            out / f"{sample}_analyzed_with_independent_niche.h5ad",
            compression="gzip"
        )

    print(
        f"[{sample}] 完成：K={selected_k}，丰度={abundance_key}，输出={out}",
        flush=True
    )
    return method


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = []
    for i, (sample, path) in enumerate(INPUTS.items()):
        if not path.exists():
            raise FileNotFoundError(path)
        summary.append(process_sample(sample, path, SEED + i * 10000))

    pd.DataFrame(summary).to_csv(
        OUT_ROOT / "independent_niche_run_summary.csv", index=False
    )
    print(f"\n全部完成：{OUT_ROOT}", flush=True)


if __name__ == "__main__":
    main()
