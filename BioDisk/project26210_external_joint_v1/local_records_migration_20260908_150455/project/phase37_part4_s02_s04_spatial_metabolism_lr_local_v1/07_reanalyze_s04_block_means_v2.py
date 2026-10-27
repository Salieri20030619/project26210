#!/usr/bin/env python3
from __future__ import annotations
import argparse, importlib.util, json, math, sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("p37",HERE/"06_analyze_s02_s04_spatial_v1.py"); p37=importlib.util.module_from_spec(spec); spec.loader.exec_module(p37)

def mean_block_rho(x,y,b,mask):
 f=pd.DataFrame({"b":b[mask],"x":np.asarray(x)[mask],"y":np.asarray(y)[mask]}).groupby("b",as_index=False).agg(x=("x","mean"),y=("y","mean"),n=("x","size"))
 if len(f)<5 or np.ptp(f.x)==0 or np.ptp(f.y)==0: return math.nan,len(f),int(mask.sum())
 return float(spearmanr(f.x,f.y).statistic),len(f),int(mask.sum())

def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--phase18-spots",type=Path,required=True); ap.add_argument("--phase19-spots",type=Path,required=True); ap.add_argument("--s04-spot-scores",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
 if a.output.exists(): raise FileExistsError(a.output)
 a.output.mkdir(parents=True); tables=a.output/"tables"; tables.mkdir()
 p18=pd.read_csv(a.phase18_spots,sep="\t"); p19=pd.read_csv(a.phase19_spots,sep="\t"); scores=pd.read_csv(a.s04_spot_scores,sep="\t")
 rows=[]; sparsity=[]
 for (cohort,sample),base0 in p18.groupby(["cohort","sample"],sort=True):
  base=base0.merge(p19[p19["sample"].eq(sample)][["spot_id","tymp_exposure"]],on="spot_id",validate="one_to_one").set_index("spot_id"); coords=base[["x","y"]].to_numpy(float); b=p37.blocks(coords); mal=base["fraction__Malignant"].to_numpy(float); macro=base[["fraction__C1QC-state Macrophage","fraction__SPP1-state Macrophage","fraction__Other Macrophage"]].sum(axis=1).to_numpy(float); endo=base["fraction__Endothelial"].to_numpy(float)
  sf=scores[(scores.cohort==cohort)&(scores["sample"]==sample)]
  for pair,(_,_,receiver,role) in p37.PAIRS.items():
   z=sf[sf.pair==pair]
   if z.empty:
    for gate in ["q25","q50"]: rows.append({"cohort":cohort,"sample":sample,"pair":pair,"role":role,"gate":gate,"eligible":False,"n_spots":0,"n_blocks":0,"n_nonzero_lr_spots":0,"fraction_nonzero_lr_spots":0.,"spearman_rho":np.nan})
    continue
   score=z.set_index("spot_id").loc[base.index,"spatial_lr_score"].to_numpy(float); rv=macro if receiver=="macrophage" else endo
   sparsity.append({"cohort":cohort,"sample":sample,"pair":pair,"n_spots":len(score),"n_nonzero_lr_spots":int((score>0).sum()),"fraction_nonzero_lr_spots":float((score>0).mean()),"maximum_lr_score":float(score.max())})
   for gate,q in [("q25",.25),("q50",.5)]:
    mask=(mal>=np.quantile(mal,q))&(rv>=np.quantile(rv,q)); rho,nb,ns=mean_block_rho(base.tymp_exposure,score,b,mask); rows.append({"cohort":cohort,"sample":sample,"pair":pair,"role":role,"gate":gate,"eligible":True,"n_spots":ns,"n_blocks":nb,"n_nonzero_lr_spots":int((score[mask]>0).sum()),"fraction_nonzero_lr_spots":float((score[mask]>0).mean()),"spearman_rho":rho})
 section=pd.DataFrame(rows); primary=section[(section.cohort=="TJRCC")&(section.gate=="q25")]; summary,loo=p37.summarize(primary,"S04"); prim=summary.role.eq("PRIMARY"); summary.loc[prim,"BH_FDR_prespecified5primary"]=p37.bh(summary.loc[prim,"wilcoxon_p_two_sided"]); summary.loc[~prim,"BH_FDR_prespecified5primary"]=np.nan; summary["statistical_support"]=np.where(summary.n_sections>=10,np.where(prim,summary.BH_FDR_prespecified5primary<.05,summary.wilcoxon_p_two_sided<.05),np.where(summary.n_sections>=3,summary.wilcoxon_p_two_sided<.05,False)); summary["direction"]=np.where(summary.median_rho>0,"positive",np.where(summary.median_rho<0,"negative","unevaluable"))
 p37.atomic_tsv(section,tables/"s04_section_statistics_block_mean_v2.tsv"); p37.atomic_tsv(summary,tables/"s04_tjrcc12_inference_block_mean_v2.tsv"); p37.atomic_tsv(loo,tables/"s04_leave_one_section_out_block_mean_v2.tsv"); p37.atomic_tsv(pd.DataFrame(sparsity),tables/"s04_lr_sparsity_audit_v2.tsv")
 formal=primary[(primary.role=="PRIMARY")&primary.eligible]; sens=primary[(primary.role=="SENSITIVITY")&primary.eligible]; ok=(formal.pair.nunique()==3 and len(formal)==36 and formal.spearman_rho.notna().all() and len(sens)==12 and sens.spearman_rho.notna().all())
 payload={"status":"PASS_S04_BLOCK_MEAN_REANALYSIS_COMPLETE_PRESENTATION_PENDING" if ok else "FAIL_S04_BLOCK_MEAN_REANALYSIS","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"supersedes_s04_in":"analysis_v1","reason":"block medians were zero-inflated and produced no estimable TJRCC section statistics","n_tjrcc_primary_evaluable_pairs":3,"n_tjrcc_primary_section_effects":int(len(formal)),"n_tjrcc_sensitivity_section_effects":int(len(sens)),"two_sided_testing":True,"direction_prespecified":False,"figures_created":False}; p37.atomic_json(payload,a.output/"status.json"); print(json.dumps(payload,indent=2)); return 0 if ok else 1
if __name__=="__main__": sys.exit(main())
