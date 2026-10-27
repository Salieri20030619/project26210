#!/usr/bin/env python
# Part4-02: 构建带 TYMP 亚型 + TAM 亚型 + 内皮的 cell2location reference 并训练 NB 回归
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
import scipy.sparse as sp, torch
import cell2location
from cell2location.models import RegressionModel

torch.set_num_threads(32)            # 256核，intra-op 32线程，留余量
OBJ="/mnt/BioDisk/ludingyu/project26210/object"
OUT="/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/ref"
os.makedirs(OUT, exist_ok=True)

# ---- 1. 读 annotation，用 counts 作 X ----
a=sc.read_h5ad(f"{OBJ}/01_annotation_adata_annotated.h5ad")
assert "counts" in a.layers, "缺 counts layer"
a.X=a.layers["counts"].copy()        # cell2location 必须 raw counts
print("annotation:", a.shape)

# ---- 2. merge 三套标签 ----
cg=pd.read_csv("/tmp/label_cellchat_group.csv", index_col=0)["cellchat_group"]
ms=pd.read_csv("/tmp/label_macro_subtype.csv", index_col=0)["macro_subtype"]
a.obs["cellchat_group"]=a.obs_names.map(cg)
a.obs["macro_subtype"]=a.obs_names.map(ms)

# 只保留有 cellchat_group 的细胞（明确归类）
a=a[~a.obs["cellchat_group"].isna()].copy()
print("保留有明确骨架标签的细胞:", a.shape)

# 建统一 c2l_celltype
def merge_label(row):
    g=row["cellchat_group"]
    if g=="Macrophage":
        sub=row["macro_subtype"]
        if pd.isna(sub): return "Other Macrophage"
        if sub=="SPP1+Macro": return "SPP1+ TAM"
        if sub=="APOE+C1QC+Macro": return "APOE+C1QC+ TAM"
        return "Other Macrophage"          # 其余5个亚型并入
    return g                                # 非巨噬原样（含TYMP拆分、内皮、正常上皮等）
a.obs["c2l_celltype"]=a.obs.apply(merge_label, axis=1).astype("category")
print("\n=== c2l_celltype 分布 ===")
print(a.obs["c2l_celltype"].value_counts())

# batch：用 study 校正数据集间技术差异
a.obs["c2l_batch"]=a.obs["study"].astype(str)
print("\n=== batch(study) ===\n", a.obs["c2l_batch"].value_counts())

# ---- 3. 基因过滤（cell2location 官方推荐的表达筛选）----
from cell2location.utils.filtering import filter_genes
sel=filter_genes(a, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
a=a[:, sel].copy()
print("\n基因过滤后:", a.shape)

# ---- 4. 训 RegressionModel（全量，CPU）----
RegressionModel.setup_anndata(adata=a, batch_key="c2l_batch", labels_key="c2l_celltype")
mod=RegressionModel(a)
mod.train(max_epochs=250, batch_size=2500,
          train_size=1.0, accelerator="cpu",
          datasplitter_kwargs={"num_workers":16, "persistent_workers":True})

# ---- 5. 导出签名 ----
a=mod.export_posterior(a, sample_kwargs={"num_samples":1000, "batch_size":2500, "accelerator":"cpu"})
mod.save(f"{OUT}/regression_model", overwrite=True)
a.write(f"{OUT}/ref_signatures.h5ad")

# 各细胞类型每基因表达签名（供映射用）
if "means_per_cluster_mu_fg" in a.varm:
    sig=a.varm["means_per_cluster_mu_fg"]
else:
    keys=[k for k in a.varm if "per_cluster_mu_fg" in k]
    sig=a.varm[keys[0]]
sig=pd.DataFrame(sig, index=a.var_names)
sig.to_csv(f"{OUT}/ref_signatures_matrix.csv")

# 训练诊断图
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
mod.plot_history(20); plt.savefig(f"{OUT}/reg_train_history.png", dpi=120, bbox_inches="tight")
print("\n完成，输出:", OUT)
print("签名矩阵:", sig.shape)
