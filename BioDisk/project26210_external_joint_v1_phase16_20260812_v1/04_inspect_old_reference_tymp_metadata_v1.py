#!/usr/bin/env python3
import anndata as ad


path = "/mnt/BioDisk/ludingyu/project26210/results/part4_spatial/ref/ref_signatures.h5ad"
obj = ad.read_h5ad(path, backed="r")
print("shape", obj.shape)
columns = [
    column for column in obj.obs.columns
    if any(key in column.lower() for key in ["tymp", "factor", "celltype", "c2l", "malig"])
]
print("columns", columns)
for column in columns:
    counts = obj.obs[column].astype(str).value_counts().head(30).to_dict()
    print("COLUMN", column, "DTYPE", obj.obs[column].dtype, "VALUES", counts)
obj.file.close()
