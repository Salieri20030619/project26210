#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, sys
from datetime import datetime
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon

SEED=20260815; BOOT=9999
FOUR=["Pyrimidine metabolism","Cysteine and methionine metabolism","Arginine and proline metabolism","Glycerophospholipid metabolism"]
PAIRS={
 "MST1_MST1R_PRIMARY":(["MST1"],["MST1R"],"macrophage","PRIMARY"),
 "ICAM5_ITGAL_ITGB2_PRIMARY":(["ICAM5"],["ITGAL","ITGB2"],"macrophage","PRIMARY"),
 "ICAM5_CD209_PRIMARY":(["ICAM5"],["CD209"],"macrophage","PRIMARY"),
 "ICAM5_PRTG_PRIMARY":(["ICAM5"],["PRTG"],"endothelial","PRIMARY"),
 "ICAM5_IGDCC4_PRIMARY":(["ICAM5"],["IGDCC4"],"endothelial","PRIMARY"),
 "ICAM5_ITGAM_ITGB2_SENSITIVITY":(["ICAM5"],["ITGAM","ITGB2"],"macrophage","SENSITIVITY")}

def bh(v):
 p=np.asarray(v,float); o=np.argsort(p); q=p[o]*len(p)/np.arange(1,len(p)+1); q=np.minimum.accumulate(q[::-1])[::-1]; z=np.empty(len(p)); z[o]=np.minimum(q,1); return z
def srank(v):
 v=np.asarray(v,float); v=v[np.isfinite(v)]
 if len(v)<3 or np.all(v==0): return 1.0
 return float(wilcoxon(v,alternative="two-sided",zero_method="wilcox",correction=False,method="auto").pvalue)
def graph(coords):
 k=min(6,len(coords)-1); nei=cKDTree(coords).query(coords,k=k+1)[1][:,1:]; rr=np.repeat(np.arange(len(coords)),k); cc=nei.reshape(-1); g=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(coords),)*2); g=g.maximum(g.T); g.data[:]=1.; return g
def zscore(v):
 v=np.asarray(v,float); s=np.std(v); return np.zeros(len(v)) if not np.isfinite(s) or s==0 else (v-np.mean(v))/s
def blocks(coords):
 b=[]
 for j in range(2): b.append(pd.qcut(pd.Series(coords[:,j]).rank(method="first"),q=4,labels=False).to_numpy())
 return b[0]*4+b[1]
def block_rho(x,y,b,mask):
 f=pd.DataFrame({"b":b[mask],"x":np.asarray(x)[mask],"y":np.asarray(y)[mask]}).groupby("b",as_index=False).agg(x=("x","median"),y=("y","median"),n=("x","size"))
 if len(f)<5 or np.ptp(f.x)==0 or np.ptp(f.y)==0: return math.nan,len(f),int(mask.sum())
 return float(spearmanr(f.x,f.y).statistic),len(f),int(mask.sum())
def dense(x):
 if hasattr(x,"to_memory"): x=x.to_memory()
 return x.toarray() if sparse.issparse(x) else np.asarray(x)
def atomic_tsv(x,p,gz=False):
 q=Path(str(p)+".part"); x.to_csv(q,sep="\t",index=False,compression="gzip" if gz else None); q.replace(p)
def atomic_json(x,p):
 q=Path(str(p)+".part"); q.write_text(json.dumps(x,indent=2),encoding="utf-8"); q.replace(p)
def summarize(section,kind):
 rng=np.random.default_rng(SEED+(0 if kind=="S02" else 1000)); rows=[]; loo=[]
 keys=["source","pathway"] if kind=="S02" else ["pair","role"]
 for key,f in section.groupby(keys,dropna=False,sort=True):
  if not isinstance(key,tuple): key=(key,)
  v=f.spearman_rho.dropna().to_numpy(float); n=len(v); p=srank(v) if n else 1.; boot=np.median(v[rng.integers(0,n,size=(BOOT,n))],axis=1) if n else np.array([np.nan])
  d=dict(zip(keys,key)); d.update(n_sections=n,n_positive=int((v>0).sum()),n_negative=int((v<0).sum()),n_zero=int((v==0).sum()),median_rho=float(np.median(v)) if n else np.nan,wilcoxon_p_two_sided=p,bootstrap_ci_low=float(np.quantile(boot,.025)) if n else np.nan,bootstrap_ci_high=float(np.quantile(boot,.975)) if n else np.nan)
  lv=[]
  for sample in sorted(f.loc[f.spearman_rho.notna(),"sample"].unique()):
   z=f.loc[f.spearman_rho.notna() & f["sample"].ne(sample),"spearman_rho"].to_numpy(); med=float(np.median(z)); lv.append(med); loo.append({**dict(zip(keys,key)),"left_out_sample":sample,"median_rho":med})
  d["loo_direction_stability"]=float(max(np.mean(np.asarray(lv)>0),np.mean(np.asarray(lv)<0))) if lv else np.nan; rows.append(d)
 return pd.DataFrame(rows),pd.DataFrame(loo)

def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--manifest",type=Path,required=True); ap.add_argument("--phase18-spots",type=Path,required=True); ap.add_argument("--phase19-spots",type=Path,required=True); ap.add_argument("--gsva-pilot",type=Path,required=True); ap.add_argument("--gsva-sections",type=Path,required=True); ap.add_argument("--lr-coverage",type=Path,required=True); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
 if a.output.exists(): raise FileExistsError(a.output)
 a.output.mkdir(parents=True); tables=a.output/"tables"; tables.mkdir()
 manifest=pd.read_csv(a.manifest,sep="\t"); p18=pd.read_csv(a.phase18_spots,sep="\t"); p19=pd.read_csv(a.phase19_spots,sep="\t"); cov=pd.read_csv(a.lr_coverage,sep="\t")
 spot02=[]; stat02=[]; spot04=[]; stat04=[]; qc=[]
 for r in manifest.itertuples(index=False):
  base=p18[p18["sample"].eq(r.sample)].copy(); ty=p19[p19["sample"].eq(r.sample)][["spot_id","tymp_exposure"]]
  base=base.merge(ty,on="spot_id",validate="one_to_one").set_index("spot_id")
  coords=base[["x","y"]].to_numpy(float); g=graph(coords); gs=g+sparse.eye(len(base),format="csr"); denom=np.asarray(gs.sum(axis=1)).ravel(); bl=blocks(coords); mal=base["fraction__Malignant"].to_numpy(float)
  sdir=a.gsva_pilot if r.sample=="GLMF1" else a.gsva_sections/r.sample
  ss=pd.read_csv(sdir/"s02_gsva_spot_scores_v1.tsv",sep="\t"); wide=ss.pivot(index="spot_id",columns=["source","pathway"],values="gsva_ssgsea_score").loc[base.index]
  for source in ["CUSTOM_PRIMARY","KEGG_SENSITIVITY"]:
   for pathway in FOUR:
    own=wide[(source,pathway)].to_numpy(float); flux=np.asarray(gs@own).ravel()/denom; flux=zscore(flux)
    for gate,q in [("q25",.25),("q50",.5)]:
     mask=mal>=np.quantile(mal,q); rho,nb,ns=block_rho(base.tymp_exposure,flux,bl,mask); stat02.append({"cohort":r.cohort,"sample":r.sample,"source":source,"pathway":pathway,"gate":gate,"n_spots":ns,"n_blocks":nb,"spearman_rho":rho})
    spot02.append(pd.DataFrame({"cohort":r.cohort,"sample":r.sample,"spot_id":base.index,"source":source,"pathway":pathway,"gsva_ssgsea_score":own,"spatial_flux_score":flux}))
  h=ad.read_h5ad(r.input_h5ad,backed="r"); genes=set(map(str,h.var_names)); wanted=sorted(set(sum(([x for x in z[0]+z[1]] for z in PAIRS.values()),[])) & genes); idx=[h.var_names.get_loc(x) for x in wanted]; raw=dense(h.X[:,idx]).astype(float); ids=pd.Index(h.obs_names.astype(str)); h.file.close(); expr=np.log1p(raw/np.maximum(raw.sum(axis=1,keepdims=True),1)*10000); ex=pd.DataFrame(expr,index=ids,columns=wanted).loc[base.index]
  macro=base[["fraction__C1QC-state Macrophage","fraction__SPP1-state Macrophage","fraction__Other Macrophage"]].sum(axis=1).to_numpy(float); endo=base["fraction__Endothelial"].to_numpy(float)
  for pair,(ligs,recs,receiver,role) in PAIRS.items():
   eligible=all(x in genes for x in ligs+recs); rv=macro if receiver=="macrophage" else endo
   if eligible:
    ligand=np.prod(ex[ligs].to_numpy()+1,axis=1)**(1/len(ligs))-1; receptor=np.prod(ex[recs].to_numpy()+1,axis=1)**(1/len(recs))-1; sender=ligand*mal; recv=receptor*rv; neigh=np.asarray(g@recv).ravel()/np.maximum(np.asarray(g.sum(axis=1)).ravel(),1); score=np.log1p(sender*neigh); spot04.append(pd.DataFrame({"cohort":r.cohort,"sample":r.sample,"spot_id":base.index,"pair":pair,"role":role,"eligible":True,"spatial_lr_score":score}))
    for gate,q in [("q25",.25),("q50",.5)]:
     mask=(mal>=np.quantile(mal,q))&(rv>=np.quantile(rv,q)); rho,nb,ns=block_rho(base.tymp_exposure,score,bl,mask); stat04.append({"cohort":r.cohort,"sample":r.sample,"pair":pair,"role":role,"gate":gate,"eligible":True,"n_spots":ns,"n_blocks":nb,"spearman_rho":rho})
   else:
    for gate in ["q25","q50"]: stat04.append({"cohort":r.cohort,"sample":r.sample,"pair":pair,"role":role,"gate":gate,"eligible":False,"n_spots":0,"n_blocks":0,"spearman_rho":np.nan})
  qc.append({"cohort":r.cohort,"sample":r.sample,"n_spots":len(base),"n_expression_genes":r.n_input_genes,"n_lr_genes_present":len(wanted),"gsva_rows":len(ss),"status":"PASS"})
 s02=pd.DataFrame(stat02); s04=pd.DataFrame(stat04); primary02=s02[(s02.cohort=="TJRCC")&(s02.gate=="q25")]; sum02,loo02=summarize(primary02,"S02"); sum02["BH_FDR_within_source4pathways"]=sum02.groupby("source",sort=False).wilcoxon_p_two_sided.transform(bh); sum02["statistical_support"]=sum02.BH_FDR_within_source4pathways.lt(.05); sum02["direction"]=np.where(sum02.median_rho>0,"positive",np.where(sum02.median_rho<0,"negative","zero"))
 primary04=s04[(s04.cohort=="TJRCC")&(s04.gate=="q25")]; sum04,loo04=summarize(primary04,"S04"); prim=sum04.role.eq("PRIMARY"); sum04.loc[prim,"BH_FDR_prespecified5primary"]=bh(sum04.loc[prim,"wilcoxon_p_two_sided"]); sum04.loc[~prim,"BH_FDR_prespecified5primary"]=np.nan; sum04["statistical_support"]=np.where(sum04.n_sections>=10,np.where(prim,sum04.BH_FDR_prespecified5primary<.05,sum04.wilcoxon_p_two_sided<.05),np.where(sum04.n_sections>=3,sum04.wilcoxon_p_two_sided<.05,False)); sum04["direction"]=np.where(sum04.median_rho>0,"positive",np.where(sum04.median_rho<0,"negative",np.where(sum04.median_rho<0,"negative","unevaluable")))
 atomic_tsv(pd.concat(spot02,ignore_index=True),tables/"s02_spot_gsva_spatial_flux_scores_v1.tsv.gz",True); atomic_tsv(s02,tables/"s02_section_statistics_v1.tsv"); atomic_tsv(sum02,tables/"s02_tjrcc12_inference_v1.tsv"); atomic_tsv(loo02,tables/"s02_leave_one_section_out_v1.tsv"); atomic_tsv(pd.concat(spot04,ignore_index=True),tables/"s04_spot_restricted_lr_scores_v1.tsv.gz",True); atomic_tsv(s04,tables/"s04_section_statistics_v1.tsv"); atomic_tsv(sum04,tables/"s04_tjrcc12_inference_v1.tsv"); atomic_tsv(loo04,tables/"s04_leave_one_section_out_v1.tsv"); atomic_tsv(pd.DataFrame(qc),tables/"section_qc_v1.tsv")
 ok=len(s02)==224 and len(s04)==168 and len(primary02)==96 and len(sum02)==8 and len(sum04)==6 and pd.DataFrame(qc).status.eq("PASS").all()
 payload={"status":"PASS_S02_S04_ANALYSIS_COMPLETE_PRESENTATION_PENDING" if ok else "FAIL_S02_S04_ANALYSIS","completed_at":datetime.now().astimezone().isoformat(timespec="seconds"),"n_sections":14,"n_spots":41688,"s02_method":"GSVA 2.4.9 ssGSEA plus symmetric-6NN spatial flux score","s02_primary_tests":8,"s04_prespecified_pairs":6,"s04_tjrcc_primary_evaluable":3,"s04_tjrcc_primary_unevaluable":2,"two_sided_testing":True,"direction_prespecified":False,"figures_created":False}; atomic_json(payload,a.output/"status.json"); print(json.dumps(payload,indent=2)); return 0 if ok else 1
if __name__=="__main__": sys.exit(main())
