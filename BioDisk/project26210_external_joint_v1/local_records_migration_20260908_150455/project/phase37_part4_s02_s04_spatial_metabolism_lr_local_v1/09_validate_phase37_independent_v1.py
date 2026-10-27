#!/usr/bin/env python3
import hashlib,json,sys
from datetime import datetime
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
ROOT=Path(__file__).resolve().parent; SYNC=ROOT/"server_results_sync_v2"; OUT=ROOT/"independent_validation_local_v1"; checks=[]
def ck(n,o,d): checks.append({"check":n,"passed":bool(o),"detail":str(d)})
def bh(v):
 p=np.asarray(v,float); o=np.argsort(p); q=p[o]*len(p)/np.arange(1,len(p)+1); q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(p)); z[o]=np.minimum(q,1); return z
def sr(v):
 v=np.asarray(v,float); v=v[np.isfinite(v)]
 return 1. if len(v)<3 or np.all(v==0) else float(wilcoxon(v,alternative="two-sided",zero_method="wilcox",correction=False,method="auto").pvalue)
def close(a,b): return np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=1e-8,atol=1e-10,equal_nan=True)
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for x in iter(lambda:f.read(1048576),b""):h.update(x)
 return h.hexdigest()
def main():
 pre=json.loads((SYNC/"preflight_final_v3"/"status.json").read_text()); ana=json.loads((SYNC/"analysis_v1"/"status.json").read_text()); final=json.loads((SYNC/"s04_final_v3"/"status.json").read_text())
 ck("final preflight PASS",pre["status"].startswith("PASS"),pre["status"]); ck("preflight checks 11/11",pre["checks_passed"]==pre["checks_total"]==11,"11/11"); ck("14 sections 41688 spots",pre["n_sections"]==14 and pre["n_spots"]==41688,"14/41688"); ck("S02/S04 analysis marker PASS",ana["status"].startswith("PASS"),ana["status"]); ck("two-sided no-direction no-figure",ana["two_sided_testing"] and not ana["direction_prespecified"] and not ana["figures_created"],"contract")
 t=SYNC/"analysis_v1"/"tables"; spot=pd.read_csv(t/"s02_spot_gsva_spatial_flux_scores_v1.tsv.gz",sep="\t"); sec=pd.read_csv(t/"s02_section_statistics_v1.tsv",sep="\t"); sm=pd.read_csv(t/"s02_tjrcc12_inference_v1.tsv",sep="\t")
 ck("S02 spot grid conserved",len(spot)==41688*8 and spot["sample"].nunique()==14 and not spot.duplicated(["sample","spot_id","source","pathway"]).any(),spot.shape); ck("S02 scores finite",np.isfinite(spot[["gsva_ssgsea_score","spatial_flux_score"]].to_numpy()).all(),"all finite"); ck("S02 section grid conserved",len(sec)==14*8*2 and sec["sample"].nunique()==14,sec.shape); ck("S02 inference has two sources four pathways",len(sm)==8 and sm.source.nunique()==2 and sm.pathway.nunique()==4,sm.shape)
 rb=[]
 for (source,pathway),f in sec[(sec.cohort=="TJRCC")&(sec.gate=="q25")].groupby(["source","pathway"],sort=True):
  v=f.spearman_rho.dropna().to_numpy(); rb.append([source,pathway,len(v),(v>0).sum(),(v<0).sum(),np.median(v),sr(v)])
 rb=pd.DataFrame(rb,columns=["source","pathway","n","pos","neg","med","p"]); rb["q"]=rb.groupby("source",sort=False).p.transform(bh); s=sm.sort_values(["source","pathway"]).reset_index(drop=True); rb=rb.sort_values(["source","pathway"]).reset_index(drop=True)
 ck("S02 counts and medians reproduced",np.array_equal(s.n_sections,rb.n) and np.array_equal(s.n_positive,rb.pos) and np.array_equal(s.n_negative,rb.neg) and close(s.median_rho,rb.med),"8/8"); ck("S02 Wilcoxon and BH reproduced",close(s.wilcoxon_p_two_sided,rb.p) and close(s.BH_FDR_within_source4pathways,rb.q),"8/8"); ck("S02 support flags follow FDR",np.array_equal(s.statistical_support.astype(bool),rb.q<.05),"8/8")
 rt=SYNC/"s04_reanalysis_v2"/"tables"; s4=pd.read_csv(rt/"s04_section_statistics_block_mean_v2.tsv",sep="\t"); m4=pd.read_csv(rt/"s04_tjrcc12_inference_block_mean_v2.tsv",sep="\t"); sp4=pd.read_csv(rt/"s04_lr_sparsity_audit_v2.tsv",sep="\t")
 ck("S04 section grid conserved",len(s4)==14*6*2 and s4["sample"].nunique()==14,s4.shape); ck("S04 all six pairs retained",m4.pair.nunique()==6 and len(m4)==6,list(m4.pair)); ck("S04 sparsity audit retained",len(sp4)>0 and (sp4.fraction_nonzero_lr_spots.between(0,1)).all(),sp4.shape)
 rb4=[]
 for (pair,role),f in s4[(s4.cohort=="TJRCC")&(s4.gate=="q25")].groupby(["pair","role"],sort=True):
  v=f.spearman_rho.dropna().to_numpy(); rb4.append([pair,role,len(v),(v>0).sum(),(v<0).sum(),np.median(v) if len(v) else np.nan,sr(v) if len(v) else 1.])
 rb4=pd.DataFrame(rb4,columns=["pair","role","n","pos","neg","med","p"]); q=bh(rb4.loc[rb4.role=="PRIMARY","p"]); rb4.loc[rb4.role=="PRIMARY","q"]=q; x=m4.sort_values("pair").reset_index(drop=True); y=rb4.sort_values("pair").reset_index(drop=True)
 ck("S04 effective section counts reproduced",np.array_equal(x.n_sections,y.n) and np.array_equal(x.n_positive,y.pos) and np.array_equal(x.n_negative,y.neg),"6/6"); ck("S04 medians and Wilcoxon reproduced",close(x.median_rho,y.med) and close(x.wilcoxon_p_two_sided,y.p),"6/6"); ck("S04 five-pair BH reproduced",close(x.loc[x.role=="PRIMARY","BH_FDR_prespecified5primary"],y.loc[y.role=="PRIMARY","q"]),"5/5"); expected=np.where(x.n_sections>=10,np.where(x.role=="PRIMARY",x.BH_FDR_prespecified5primary<.05,x.wilcoxon_p_two_sided<.05),np.where(x.n_sections>=3,x.wilcoxon_p_two_sided<.05,False)); ck("S04 reporting rule reproduced",np.array_equal(x.statistical_support.astype(bool),expected),"6/6"); ck("S04 no supported pair",not x.statistical_support.astype(bool).any(),"0/6"); ck("S04 finalizer PASS 10/10",final["status"].startswith("PASS") and final["checks_passed"]==final["checks_total"]==10,final["status"]); ck("invalid median and strict-v2 attempts preserved",(SYNC/"analysis_v1").exists() and (SYNC/"s04_reanalysis_v2"/"status.json").exists(),"preserved")
 figs=[p for p in ROOT.rglob("*") if p.suffix.lower() in {".png",".jpg",".jpeg",".pdf",".svg",".tif",".tiff"}]; ck("no figures generated",not figs,len(figs))
 OUT.mkdir(parents=False,exist_ok=False); df=pd.DataFrame(checks); df.to_csv(OUT/"phase37_independent_validation_checks_v1.tsv",sep="\t",index=False); n=int(df.passed.sum()); files={str(p.relative_to(ROOT)):sha(p) for p in SYNC.rglob("*") if p.is_file()}; payload={"status":"PASS_PHASE37_S02_S04_INDEPENDENT_VALIDATION" if n==len(df) else "FAIL_PHASE37_S02_S04_INDEPENDENT_VALIDATION","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"checks_passed":n,"checks_total":len(df),"validated_files":files,"figures_created":False}; (OUT/"status.json").write_text(json.dumps(payload,indent=2),encoding="utf-8"); print(json.dumps(payload,indent=2)); return 0 if n==len(df) else 1
if __name__=="__main__":sys.exit(main())
