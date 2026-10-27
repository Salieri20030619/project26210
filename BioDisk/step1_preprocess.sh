#!/bin/bash
# ============================================================
# Step 1: 空间数据预处理（单样本版，手动构建 Visium 结构）
# 用法：bash step1_preprocess.sh GLMF1
#       bash step1_preprocess.sh GLMF2
# 功能：解压数据 -> 手动构建 AnnData -> 质控过滤 -> 保存 .h5ad
# 输出：results/part4_spatial/preprocessed/${SAMPLE}_preprocessed.h5ad
# ============================================================

set -euo pipefail

SAMPLE=${1:-}
if [ -z "$SAMPLE" ]; then
    echo "错误：请指定样本名，例如：bash step1_preprocess.sh GLMF1"
    exit 1
fi
if [[ "$SAMPLE" != "GLMF1" && "$SAMPLE" != "GLMF2" ]]; then
    echo "错误：样本名必须为 GLMF1 或 GLMF2"
    exit 1
fi

echo "=== 开始处理样本：${SAMPLE} ==="

BASE_DATA="/mnt/BioDisk/ludingyu/project26210/external/GSE210041/unpacked"
OUTDIR="/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/preprocessed"
mkdir -p "$OUTDIR"

SAMPLE_DIR="${BASE_DATA}/${SAMPLE}"
mkdir -p "$SAMPLE_DIR"

# -------- 1. 准备数据文件（链接 H5，解压 spatial） --------
echo ">>> 检查并准备数据文件..."

H5_FILE=$(find "$BASE_DATA" -maxdepth 1 -name "*${SAMPLE}_filtered_feature_bc_matrix.h5" | head -n1)
if [ -z "$H5_FILE" ]; then
    echo "错误：未找到 ${SAMPLE} 的 filtered_feature_bc_matrix.h5 文件"
    exit 1
fi

if [ ! -L "${SAMPLE_DIR}/filtered_feature_bc_matrix.h5" ] && [ ! -f "${SAMPLE_DIR}/filtered_feature_bc_matrix.h5" ]; then
    echo "  链接 H5 文件：$H5_FILE -> ${SAMPLE_DIR}/filtered_feature_bc_matrix.h5"
    ln -s "$H5_FILE" "${SAMPLE_DIR}/filtered_feature_bc_matrix.h5"
fi

SPATIAL_TAR=$(find "$BASE_DATA" -maxdepth 1 -name "*${SAMPLE}_spatial.tar.gz" | head -n1)
if [ -z "$SPATIAL_TAR" ]; then
    echo "错误：未找到 ${SAMPLE} 的 spatial.tar.gz 文件"
    exit 1
fi

if [ ! -d "${SAMPLE_DIR}/spatial" ]; then
    echo "  解压空间图像/坐标数据：$SPATIAL_TAR"
    TMP_EXTRACT="${BASE_DATA}/tmp_extract_${SAMPLE}"
    mkdir -p "$TMP_EXTRACT"
    tar -xzf "$SPATIAL_TAR" -C "$TMP_EXTRACT"
    if [ -d "${TMP_EXTRACT}/spatial" ]; then
        mv "${TMP_EXTRACT}/spatial" "${SAMPLE_DIR}/"
    else
        mkdir -p "${SAMPLE_DIR}/spatial"
        mv "${TMP_EXTRACT}"/* "${SAMPLE_DIR}/spatial/" 2>/dev/null || true
    fi
    rm -rf "$TMP_EXTRACT"
else
    echo "  spatial 目录已存在，跳过解压"
fi

echo ">>> 数据准备完成，目录结构："
ls -lh "${SAMPLE_DIR}/"

# -------- 2. 运行 Python 预处理（手动构建，避免 squidpy 版本问题） --------
echo ">>> 开始运行 Python 预处理脚本..."

python3 << EOF
import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

sample = "$SAMPLE"
sample_dir = Path("$SAMPLE_DIR")
out_dir = Path("$OUTDIR")

print(f"读取样本：{sample}")
print(f"数据目录：{sample_dir}")

# ---------- 2.1 读取矩阵 ----------
h5_path = sample_dir / "filtered_feature_bc_matrix.h5"
adata = sc.read_10x_h5(h5_path)
if 'gene_symbols' in adata.var.columns:
    adata.var_names = adata.var['gene_symbols']
    adata.var.drop(columns='gene_symbols', inplace=True)

# ---------- 2.2 读取空间坐标 ----------
spatial_dir = sample_dir / "spatial"
pos_file = list(spatial_dir.glob("tissue_positions*.csv")) + list(spatial_dir.glob("tissue_positions*.txt"))
if not pos_file:
    raise FileNotFoundError(f"未找到 tissue_positions 文件，在 {spatial_dir} 中")
pos_df = pd.read_csv(pos_file[0], header=None, index_col=0)
pos_df.columns = ["in_tissue", "array_row", "array_col", "px", "py"]
pos_df = pos_df[pos_df["in_tissue"] == 1]
common_barcodes = adata.obs_names.intersection(pos_df.index)
adata = adata[common_barcodes, :]
spatial_coords = pos_df.loc[common_barcodes, ["px", "py"]].values
adata.obsm['spatial'] = spatial_coords

# ---------- 2.3 读取图像和缩放因子（转换为字符串，避免序列化问题） ----------
sf_file = spatial_dir / "scalefactors_json.json"
if sf_file.exists():
    with open(sf_file, 'r') as f:
        scalefactors = json.load(f)
else:
    scalefactors = {}

# 查找图像文件，转换为字符串
img_hires = list(spatial_dir.glob("tissue_hires_image.png"))
img_lowres = list(spatial_dir.glob("tissue_lowres_image.png"))
hires_path = str(img_hires[0]) if img_hires else ""
lowres_path = str(img_lowres[0]) if img_lowres else ""

# 构建 adata.uns['spatial']
adata.uns['spatial'] = {
    sample: {
        'images': {
            'hires': hires_path,
            'lowres': lowres_path,
        },
        'scalefactors': scalefactors,
    }
}

# ---------- 2.4 质控过滤 ----------
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

print(f"原始组织 spots: {adata.n_obs}, 基因: {adata.n_vars}")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=6000)
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
sc.pp.filter_genes(adata, min_cells=5)

print(f"过滤后 spots: {adata.n_obs}, 基因: {adata.n_vars}")

# ---------- 2.5 保存 ----------
adata.obs['sample'] = sample
out_path = out_dir / f"{sample}_preprocessed.h5ad"
adata.write(out_path)
print(f"预处理完成，保存至：{out_path}")
print(f"文件大小：{out_path.stat().st_size / (1024**3):.2f} GB")
print("空间坐标预览 (前5行):")
print(adata.obsm['spatial'][:5])
EOF

# -------- 3. 检查产出 --------
echo ">>> 预处理完成，产出文件列表："
ls -lh "$OUTDIR/${SAMPLE}_preprocessed.h5ad"

echo "=== 样本 ${SAMPLE} Step 1 全部完成 ==="
