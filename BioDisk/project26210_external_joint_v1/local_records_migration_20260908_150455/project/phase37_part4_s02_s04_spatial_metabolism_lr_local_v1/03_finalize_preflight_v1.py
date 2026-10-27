#!/usr/bin/env python3
import argparse, json, sys
from datetime import datetime
from pathlib import Path
import pandas as pd

ap=argparse.ArgumentParser(); ap.add_argument("--preflight",type=Path,required=True); ap.add_argument("--lr-audit",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
if a.output.exists(): raise FileExistsError(a.output)
a.output.mkdir(parents=True)
h=pd.read_csv(a.preflight/"spatial_h5ad_preflight_v1.tsv",sep="\t"); c=pd.read_csv(a.preflight/"pathway_gene_coverage_v1.tsv",sep="\t"); b=pd.read_csv(a.preflight/"boundary_annotation_preflight_v1.tsv",sep="\t"); lr=pd.read_csv(a.lr_audit/"lr_pair_scope_summary_v1.tsv",sep="\t")
primary=lr[(lr.cohort=="TJRCC") & lr.pair.str.endswith("PRIMARY")]; sensitivity=lr[(lr.cohort=="TJRCC") & lr.pair.str.endswith("SENSITIVITY")]
checks=[
 ("14 full-gene raw-count sections",len(h)==14 and h.n_spots.sum()==41688),
 ("source hashes match frozen Phase17 preflights",h.hash_match.all()),
 ("Phase18 and Phase19 spot identities match",h.phase18_spot_match.all() and h.phase19_spot_match.all()),
 ("candidate counts finite nonnegative integer",h.candidate_values_finite.all() and h.candidate_values_nonnegative.all() and h.candidate_values_integer.all()),
 ("four pathways in two sources",len(c)==112 and c.pathway.nunique()==4 and c.source.nunique()==2),
 ("all pathway sets have >=10 covered genes",(c.n_present>=10).all()),
 ("at least three primary LR pairs cover all TJRCC12",(primary.n_eligible==12).sum()>=3),
 ("sensitivity LR pair covers all TJRCC12",len(sensitivity)==1 and sensitivity.n_eligible.iloc[0]==12),
 ("GLMF1/2 cover all six LR pairs",(lr[lr.cohort=="GLMF"].n_eligible==2).all()),
 ("formal boundary scope explicitly zero",b.formal_pathology_ready.sum()==0),
 ("two AI-draft boundary sensitivity sections retained",(b.has_tumor & b.has_nat).sum()==2),
]
df=pd.DataFrame(checks,columns=["check","passed"]); df["detail"]=["" for _ in range(len(df))]; p=a.output/"preflight_final_checks_v1.tsv.part"; df.to_csv(p,sep="\t",index=False); p.replace(a.output/"preflight_final_checks_v1.tsv")
n=int(df.passed.sum()); payload={"status":"PASS_S02_S04_SPATIAL_PREFLIGHT_WITH_PARTIAL_LR_COVERAGE" if n==len(df) else "FAIL_S02_S04_SPATIAL_PREFLIGHT_FINAL","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(df),"n_sections":14,"n_spots":41688,"primary_lr_pairs_formally_evaluable":primary.loc[primary.n_eligible.eq(12),"pair"].tolist(),"primary_lr_pairs_not_evaluable_in_tjrcc":primary.loc[primary.n_eligible.ne(12),"pair"].tolist(),"sensitivity_lr_pairs_evaluable":sensitivity.loc[sensitivity.n_eligible.eq(12),"pair"].tolist(),"formal_pathology_boundary_sections":0,"draft_boundary_sensitivity_sections":2,"figures_created":False}
q=a.output/"status.json.part"; q.write_text(json.dumps(payload,indent=2),encoding="utf-8"); q.replace(a.output/"status.json"); print(json.dumps(payload,indent=2)); sys.exit(0 if n==len(df) else 1)
