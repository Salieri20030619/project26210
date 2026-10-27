#!/usr/bin/env python

import os
import pandas as pd
import scipy.io
import anndata as ad


out_dir = os.environ.get("OUT_DIR", "results/scDRS_real_GWAS_20260713")
input_dir = os.path.join(out_dir, "inputs")

mtx_file = os.path.join(input_dir, "malignant_epithelial_counts_cells_by_genes.mtx")
obs_file = os.path.join(input_dir, "malignant_epithelial_obs.tsv")
genes_file = os.path.join(input_dir, "malignant_epithelial_genes.tsv")
h5ad_file = os.path.join(input_dir, "malignant_epithelial_counts.h5ad")

print(f"[h5ad] Reading {mtx_file}")
X = scipy.io.mmread(mtx_file).tocsr()
obs = pd.read_csv(obs_file, sep="\t")
var = pd.read_csv(genes_file, sep="\t")

obs.index = obs["cell"].astype(str)
var.index = var["gene"].astype(str)
obs.index.name = None
var.index.name = None

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.var_names_make_unique()
print(f"[h5ad] Writing {h5ad_file}: cells={adata.n_obs}, genes={adata.n_vars}")
adata.write_h5ad(h5ad_file, compression="gzip")
print("[h5ad] Done.")
