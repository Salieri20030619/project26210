#!/usr/bin/env python3
import argparse, json, sys
from datetime import datetime
from pathlib import Path
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser(); ap.add_argument("--reanalyzed",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
if a.output.exists(): raise FileExistsError(a.output)
a.output.mkdir(parents=True)
t=a.reanalyzed/"tables"; sec=pd.read_csv(t/"s04_section_statistics_block_mean_v2.tsv",sep="\t"); sm=pd.read_csv(t/"s04_tjrcc12_inference_block_mean_v2.tsv",sep="\t")
expected={"ICAM5_ITGAL_ITGB2_PRIMARY":10,"ICAM5_CD209_PRIMARY":7,"ICAM5_PRTG_PRIMARY":2,"ICAM5_IGDCC4_PRIMARY":0,"MST1_MST1R_PRIMARY":0,"ICAM5_ITGAM_ITGB2_SENSITIVITY":10}
actual=dict(zip(sm.pair,sm.n_sections))
checks=[
 ("all six prespecified pairs retained",set(actual)==set(expected)),
 ("effective section counts reproduced",all(int(actual[k])==v for k,v in expected.items())),
 ("primary family retains five rows",(sm.role=="PRIMARY").sum()==5),
 ("unevaluable primary pairs assigned P=1 and FDR=1",((sm.n_sections==0)&(sm.role=="PRIMARY")&(sm.wilcoxon_p_two_sided==1)&(sm.BH_FDR_prespecified5primary==1)).sum()==2),
 ("n>=10 primary inference uses five-pair FDR",bool((sm.loc[sm.pair=="ICAM5_ITGAL_ITGB2_PRIMARY","statistical_support"]==False).all())),
 ("n=7 primary inference uses nominal P",float(sm.loc[sm.pair=="ICAM5_CD209_PRIMARY","wilcoxon_p_two_sided"].iloc[0])>=.05),
 ("n=2 pair has no statistical support",bool((sm.loc[sm.pair=="ICAM5_PRTG_PRIMARY","statistical_support"]==False).all())),
 ("sensitivity pair explicitly labeled",sm.loc[sm.pair=="ICAM5_ITGAM_ITGB2_SENSITIVITY","role"].iloc[0]=="SENSITIVITY"),
 ("no S04 pair supported",not sm.statistical_support.astype(bool).any()),
 ("section grid retained",len(sec)==168 and sec["sample"].nunique()==14),
]
df=pd.DataFrame(checks,columns=["check","passed"]); df["detail"]=""; p=a.output/"s04_expression_limited_final_checks_v3.tsv.part"; df.to_csv(p,sep="\t",index=False); p.replace(a.output/"s04_expression_limited_final_checks_v3.tsv")
n=int(df.passed.sum()); payload={"status":"PASS_S04_ANALYSIS_COMPLETE_EXPRESSION_LIMITED_PRESENTATION_PENDING" if n==len(df) else "FAIL_S04_FINALIZATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(df),"authoritative_analysis":str(a.reanalyzed),"effective_tjrcc_sections":expected,"statistically_supported_pairs":[],"two_sided_testing":True,"direction_prespecified":False,"figures_created":False}
q=a.output/"status.json.part"; q.write_text(json.dumps(payload,indent=2),encoding="utf-8"); q.replace(a.output/"status.json"); print(json.dumps(payload,indent=2)); sys.exit(0 if n==len(df) else 1)
