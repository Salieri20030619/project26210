#!/usr/bin/env python3
"""Build the Phase37 full-gene raw-count manifest from frozen Phase17 preflights."""
import argparse
from pathlib import Path
import pandas as pd

ap=argparse.ArgumentParser(); ap.add_argument("--glmf",type=Path,required=True); ap.add_argument("--tjrcc",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
if a.output.exists(): raise FileExistsError(a.output)
a.output.mkdir(parents=True)
g=pd.read_csv(a.glmf,sep="\t"); g.insert(0,"cohort","GLMF")
t=pd.read_csv(a.tjrcc,sep="\t"); t.insert(0,"cohort","TJRCC")
x=pd.concat([g,t],ignore_index=True,sort=False).rename(columns={"input_sha256":"input_h5ad_sha256"})
keep=["cohort","sample","input_h5ad","input_h5ad_sha256","n_spots","n_input_genes","sample_counts_nonnegative","sample_counts_integer_like","spatial_shape"]
x=x[keep]
part=a.output/"full_expression_input_manifest_v1.tsv.part"; x.to_csv(part,sep="\t",index=False); part.replace(a.output/"full_expression_input_manifest_v1.tsv")
print(x[["cohort","sample","n_spots","n_input_genes"]].to_string(index=False))
