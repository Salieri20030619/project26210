#!/usr/bin/env python3
import json
import sys
import anndata as ad

p = sys.argv[1]
a = ad.read_h5ad(p, backed="r")
out = {
    "shape": list(a.shape),
    "obs_columns": list(a.obs.columns),
    "layers": list(a.layers.keys()),
    "raw_present": a.raw is not None,
    "x_dtype": str(a.X.dtype),
    "obs_preview": {c: a.obs[c].astype(str).value_counts().head(20).to_dict() for c in a.obs.columns if any(k in c.lower() for k in ["cell", "annot", "type", "malig", "dataset", "sample", "patient"])},
}
print(json.dumps(out, ensure_ascii=False, indent=2))
a.file.close()
