#!/usr/bin/env python3
import argparse
from pathlib import Path
import anndata as ad
import pandas as pd

PAIRS={"MST1_MST1R_PRIMARY":["MST1","MST1R"],"ICAM5_ITGAL_ITGB2_PRIMARY":["ICAM5","ITGAL","ITGB2"],"ICAM5_CD209_PRIMARY":["ICAM5","CD209"],"ICAM5_PRTG_PRIMARY":["ICAM5","PRTG"],"ICAM5_IGDCC4_PRIMARY":["ICAM5","IGDCC4"],"ICAM5_ITGAM_ITGB2_SENSITIVITY":["ICAM5","ITGAM","ITGB2"]}
ap=argparse.ArgumentParser(); ap.add_argument("--manifest",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
if a.output.exists(): raise FileExistsError(a.output)
a.output.mkdir(parents=True)
m=pd.read_csv(a.manifest,sep="\t"); rows=[]
for r in m.itertuples(index=False):
    h=ad.read_h5ad(r.input_h5ad,backed="r"); genes=set(map(str,h.var_names)); h.file.close()
    for pair,required in PAIRS.items():
        present=[g for g in required if g in genes]; missing=[g for g in required if g not in genes]
        rows.append({"cohort":r.cohort,"sample":r.sample,"pair":pair,"required_genes":";".join(required),"present_genes":";".join(present),"missing_genes":";".join(missing),"eligible":not missing})
x=pd.DataFrame(rows); p=a.output/"lr_pair_gene_coverage_v1.tsv.part"; x.to_csv(p,sep="\t",index=False); p.replace(a.output/"lr_pair_gene_coverage_v1.tsv")
s=x.groupby(["cohort","pair"],as_index=False).agg(n_sections=("sample","size"),n_eligible=("eligible","sum")); p=a.output/"lr_pair_scope_summary_v1.tsv.part"; s.to_csv(p,sep="\t",index=False); p.replace(a.output/"lr_pair_scope_summary_v1.tsv"); print(s.to_string(index=False))
