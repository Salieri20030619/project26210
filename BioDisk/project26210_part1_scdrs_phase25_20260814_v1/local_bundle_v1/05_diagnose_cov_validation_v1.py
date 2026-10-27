#!/usr/bin/env python3
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

parser = argparse.ArgumentParser()
parser.add_argument("--input-dir", required=True, type=Path)
args = parser.parse_args()

data = ad.read_h5ad(args.input_dir / "malignant_eligible_raw_counts_phase25_v1.h5ad", backed="r")
ids = pd.Index(data.obs_names.astype(str))
library = []
detected = []
for start in range(0, data.n_obs, 512):
    block = sp.csr_matrix(data.X[start:start + 512, :])
    library.append(np.asarray(block.sum(axis=1)).ravel())
    detected.append(np.diff(block.indptr))
library = np.concatenate(library)
detected = np.concatenate(detected)
cov = pd.read_csv(args.input_dir / "malignant_eligible_phase25_v1.cov", sep="\t", index_col=0)
cov.index = cov.index.astype(str)
delta = cov["log_n_counts"].to_numpy(float) - np.log1p(library)
print("index_equal", cov.index.equals(ids))
print("index_set_equal", set(cov.index) == set(ids))
print("n_genes_equal", np.array_equal(cov["n_genes"].to_numpy(int), detected))
print("log_allclose_1e12", np.allclose(delta, 0, rtol=0, atol=1e-12))
print("log_max_abs_diff", float(np.max(np.abs(delta))))
print("log_nonzero_count", int((delta != 0).sum()))
print("const_unique", sorted(cov["const"].unique().tolist()))
data.file.close()
