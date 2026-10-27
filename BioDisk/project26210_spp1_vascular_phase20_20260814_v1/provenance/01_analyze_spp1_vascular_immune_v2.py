#!/usr/bin/env python3
from __future__ import annotations
import argparse,gzip,hashlib,json,math,os,platform
from datetime import datetime
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,scipy
from scipy.sparse import issparse
from scipy.stats import binomtest,spearmanr,wilcoxon

SEED=20260814; BOOTSTRAPS=9999
SIGS=["cytotoxicity","hypoxia_angiogenesis","invasion_emt"]
PRIMARY=[("SPP1_polarization_to_vascular_compartment",1),("SPP1_vascular_niche_to_T_cell_balance",-1),("SPP1_vascular_niche_to_cytotoxicity",-1),("SPP1_vascular_niche_to_hypoxia_angiogenesis",1),("SPP1_vascular_niche_to_invasion_EMT",1)]
def now():return datetime.now().astimezone().isoformat(timespec="seconds")
def sha256(p):
 d=hashlib.sha256()
 with Path(p).open("rb") as h:
  for b in iter(lambda:h.read(8*1024*1024),b""):d.update(b)
 return d.hexdigest()
def z(x):
 x=np.asarray(x,float);s=np.nanstd(x);return np.zeros(len(x)) if not np.isfinite(s) or s<=0 else (x-np.nanmean(x))/s
def resid(y,c):
 x=np.column_stack([np.ones(len(y)),*[z(v) for v in c]]);t=z(y);return z(t-x@np.linalg.lstsq(x,t,rcond=None)[0])
def rho(x,y):
 x,y=np.asarray(x,float),np.asarray(y,float);o=np.isfinite(x)&np.isfinite(y);return math.nan if o.sum()<5 or np.ptp(x[o])==0 or np.ptp(y[o])==0 else float(spearmanr(x[o],y[o]).statistic)
def percentile(x):return pd.Series(x).rank(method="average",pct=True).to_numpy(float)
def bh(p):
 p=np.asarray(p,float);o=np.argsort(p);q=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];r=np.empty(len(p));r[o]=np.minimum(q,1);return r
def atomic_tsv(f,p):
 part=Path(str(p)+".part")
 if str(p).endswith(".gz"):
  with gzip.open(part,"wt",encoding="utf-8",newline="") as h:f.to_csv(h,sep="\t",index=False,lineterminator="\n")
 else:f.to_csv(part,sep="\t",index=False,lineterminator="\n")
 os.replace(part,p)
def atomic_json(p,x):
 q=Path(str(p)+".part");q.write_text(json.dumps(x,indent=2)+"\n",encoding="utf-8");os.replace(q,p)
def blocks(coords):
 b=[]
 for j in range(2):b.append(pd.qcut(pd.Series(coords[:,j]).rank(method="first"),q=4,labels=False).to_numpy())
 return b[0]*4+b[1]
def blockrho(x,y,b,m):
 f=pd.DataFrame({"b":b[m],"x":np.asarray(x)[m],"y":np.asarray(y)[m]}).groupby("b",as_index=False).median(numeric_only=True);return rho(f.x,f.y),len(f),int(m.sum())
def discover(root):
 r=[]
 for cohort,sub in [("GLMF","glmf12_cell2location_tymp_continuous_v1"),("TJRCC","tjrcc12_cell2location_tymp_continuous_v1")]:
  for sp in sorted((root/sub).glob("*/status.json")):
   s=json.loads(sp.read_text());h=Path(s["output_h5ad"])
   if s["status"]!="completed" or int(s["epoch"])!=10000 or h.stat().st_size!=int(s["output_size_bytes"]) or sha256(h)!=s["output_sha256"]:raise RuntimeError(f"Phase17 gate {sp}")
   r.append((cohort,s["sample"],h,s["output_sha256"]))
 if len(r)!=14:return (_ for _ in ()).throw(RuntimeError("14 inputs required"))
 return r
def score_programs(h5ad,genes):
 a=ad.read_h5ad(h5ad)
 try:
  wanted=genes[genes.signature_id.isin(SIGS)].gene.drop_duplicates().tolist();present=[g for g in wanted if g in a.var_names];idx=[a.var_names.get_loc(g) for g in present];x=a.X[:,idx]
  if issparse(x):x=x.toarray()
  x=np.asarray(x,float);tot=np.asarray(a.X.sum(axis=1)).ravel().astype(float);ln=np.log1p(x*(10000/tot[:,None]));pos={g:i for i,g in enumerate(present)};spots=a.obs_names.astype(str).tolist()
 finally:del a
 out=pd.DataFrame({"spot_id":spots,"total_counts":tot});qc=[]
 for sid in SIGS:
  d=genes[genes.signature_id==sid].sort_values("gene_order");use=[g for g in d.gene if g in pos];minimum=int(d.minimum_present_genes.iloc[0])
  if len(use)<minimum:raise RuntimeError(f"coverage {sid}")
  out[sid]=np.column_stack([z(ln[:,pos[g]]) for g in use]).mean(axis=1);qc.append({"signature_id":sid,"defined":len(d),"minimum_present":minimum,"present":len(use),"missing":";".join([g for g in d.gene if g not in pos])})
 return out,pd.DataFrame(qc)
def main():
 p=argparse.ArgumentParser();p.add_argument("--phase17",type=Path,required=True);p.add_argument("--phase18",type=Path,required=True);p.add_argument("--phase18-validation",type=Path,required=True);p.add_argument("--phase19",type=Path,required=True);p.add_argument("--phase19-validation",type=Path,required=True);p.add_argument("--signature-genes",type=Path,required=True);p.add_argument("--signature-registry",type=Path,required=True);p.add_argument("--analysis-registry",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
 if a.output.exists():raise FileExistsError(a.output)
 a.output.mkdir(parents=True);t=a.output/"tables";t.mkdir();atomic_json(a.output/"RUNNING.json",{"status":"RUNNING","started_at":now()})
 if json.loads((a.phase18_validation/"PHASE18_INDEPENDENT_VALIDATION_SUMMARY_v1.json").read_text())["status"]!="PASS" or json.loads((a.phase19_validation/"PHASE19_INDEPENDENT_VALIDATION_SUMMARY_v1.json").read_text())["status"]!="PASS":raise RuntimeError("upstream validation gate")
 wide=pd.read_csv(a.phase18/"tables/spot_composition_wide.tsv.gz",sep="\t");p19=pd.read_csv(a.phase19/"tables/spot_scores_v1.tsv.gz",sep="\t",usecols=["sample","spot_id","neighbor_tymp_exposure","neighbor_mta_context"]);reg=pd.read_csv(a.signature_registry,sep="\t");genes=pd.read_csv(a.signature_genes,sep="\t").merge(reg[["signature_id","minimum_present_genes"]],on="signature_id")
 stats=[];spots=[];qcs=[];niches=[];prov=[]
 for cohort,sample,h,hsh in discover(a.phase17):
  print(f"[{now()}] {cohort} {sample}",flush=True);f=wide[wide["sample"]==sample].copy();sc,qc=score_programs(h,genes);qc.insert(0,"sample",sample);qc.insert(0,"cohort",cohort);qcs.append(qc);f=f.merge(sc,on="spot_id",validate="one_to_one").merge(p19[p19["sample"]==sample].drop(columns="sample"),on="spot_id",validate="one_to_one")
  macro=f[["fraction__C1QC-state Macrophage","fraction__SPP1-state Macrophage","fraction__Other Macrophage"]].sum(axis=1).to_numpy();vascular_total=f[["fraction__Endothelial","fraction__SMC/Pericyte"]].sum(axis=1).to_numpy();spp1=f["fraction__SPP1-state Macrophage"].to_numpy()/np.maximum(macro,1e-12);endo=f["fraction__Endothelial"].to_numpy();smc=f["fraction__SMC/Pericyte"].to_numpy();vascular=np.sqrt(percentile(endo)*percentile(smc));niche=np.cbrt(percentile(spp1)*percentile(endo)*percentile(smc));tb=f["fraction__T cell"].to_numpy()/np.maximum(f[["fraction__T cell","fraction__NK cell","fraction__B cell"]].sum(axis=1).to_numpy(),1e-12);depth=np.log1p(f.total_counts.to_numpy());mal=f["fraction__Malignant"].to_numpy();tfrac=f["fraction__T cell"].to_numpy()
  cy=resid(f.cytotoxicity,[depth,tfrac]);hy=resid(f.hypoxia_angiogenesis,[depth,mal]);iv=resid(f.invasion_emt,[depth,mal]);b=blocks(f[["x","y"]].to_numpy());mask=(macro>=np.quantile(macro,.25))&(vascular_total>=np.quantile(vascular_total,.25))
  rel=[("SPP1_polarization_to_vascular_compartment",spp1,vascular,1,"primary"),("SPP1_vascular_niche_to_T_cell_balance",niche,tb,-1,"primary"),("SPP1_vascular_niche_to_cytotoxicity",niche,cy,-1,"primary"),("SPP1_vascular_niche_to_hypoxia_angiogenesis",niche,hy,1,"primary"),("SPP1_vascular_niche_to_invasion_EMT",niche,iv,1,"primary"),("TYMP_context_to_SPP1_vascular_niche",f.neighbor_tymp_exposure,niche,1,"secondary_context"),("MTA_context_to_SPP1_vascular_niche",f.neighbor_mta_context,niche,1,"secondary_context")]
  for name,x,y,direction,role in rel:
   r,nb,ns=blockrho(x,y,b,mask);stats.append({"cohort":cohort,"sample":sample,"relation":name,"expected_direction":direction,"role":role,"gate":"macro_and_vascular_q25","grid":"4x4","n_spots":ns,"n_blocks":nb,"spearman_rho":r})
  f["spp1_polarization"]=spp1;f["vascular_compartment"]=vascular;f["spp1_vascular_niche"]=niche;f["t_cell_balance"]=tb;f["cytotoxicity_residual"]=cy;f["hypoxia_residual"]=hy;f["invasion_residual"]=iv
  for n,g in f.groupby("niche"):
   niches.append({"cohort":cohort,"sample":sample,"niche":n,"n_spots":len(g),"median_spp1_vascular_niche":float(g.spp1_vascular_niche.median()),"median_t_cell_balance":float(g.t_cell_balance.median()),"median_cytotoxicity_residual":float(g.cytotoxicity_residual.median()),"median_hypoxia_residual":float(g.hypoxia_residual.median()),"median_invasion_residual":float(g.invasion_residual.median())})
  spots.append(f[["cohort","sample","spot_id","niche","spp1_polarization","vascular_compartment","spp1_vascular_niche","t_cell_balance","cytotoxicity_residual","hypoxia_residual","invasion_residual","neighbor_tymp_exposure","neighbor_mta_context"]]);prov.append({"cohort":cohort,"sample":sample,"h5ad":str(h),"sha256":hsh,"n_spots":len(f)})
 section=pd.DataFrame(stats);primary=section[(section.cohort=="TJRCC")&(section.role=="primary")];rng=np.random.default_rng(SEED);agg=[];loo=[]
 if len(primary)!=60:raise RuntimeError("primary rows")
 for relation,direction in PRIMARY:
  v=primary[primary.relation==relation].spearman_rho.to_numpy();oriented=v*direction;boot=np.median(v[rng.integers(0,12,size=(BOOTSTRAPS,12))],axis=1)
  try:wp=float(wilcoxon(v).pvalue)
  except ValueError:wp=1.;sp=float(binomtest(int((oriented>0).sum()),12,.5).pvalue);lv=[]
  sp=float(binomtest(int((oriented>0).sum()),12,.5).pvalue);lv=[]
  for sample in sorted(primary["sample"].unique()):
   med=float(np.median(primary[(primary.relation==relation)&(primary["sample"]!=sample)].spearman_rho));lv.append(med*direction);loo.append({"relation":relation,"left_out_sample":sample,"median_rho":med,"expected_direction":direction})
  agg.append({"relation":relation,"expected_direction":direction,"sections":12,"median_rho":float(np.median(v)),"expected_direction_sections":int((oriented>0).sum()),"opposite_direction_sections":int((oriented<0).sum()),"wilcoxon_p":wp,"sign_test_p":sp,"bootstrap_ci_low":float(np.quantile(boot,.025)),"bootstrap_ci_high":float(np.quantile(boot,.975)),"loo_expected_direction_fraction":float((np.array(lv)>0).mean())})
 aggregate=pd.DataFrame(agg);aggregate["wilcoxon_fdr_bh"]=bh(aggregate.wilcoxon_p);support=(aggregate.median_rho*aggregate.expected_direction>=.1)&(aggregate.expected_direction_sections>=8)&(aggregate.wilcoxon_fdr_bh<.05)&(aggregate.bootstrap_ci_low*aggregate.expected_direction>0)&(aggregate.bootstrap_ci_high*aggregate.expected_direction>0)&(aggregate.loo_expected_direction_fraction==1);aggregate["decision"]=np.where(support,"PRIMARY_LAYER_SUPPORT",np.where(aggregate.expected_direction_sections>=8,"DIRECTIONAL_NOT_FDR",np.where((aggregate.expected_direction_sections>=4)&(aggregate.opposite_direction_sections>=4),"HETEROGENEOUS","NOT_SUPPORTED")))
 outs={"spot_vascular_immune_scores_v1.tsv.gz":pd.concat(spots,ignore_index=True),"section_statistics_v1.tsv":section,"primary_tjrcc_section_statistics_v1.tsv":primary,"primary_across_section_inference_v1.tsv":aggregate,"primary_leave_one_section_out_v1.tsv":pd.DataFrame(loo),"niche_vascular_immune_summary_v1.tsv":pd.DataFrame(niches),"signature_gene_coverage_v1.tsv":pd.concat(qcs,ignore_index=True),"input_provenance_v1.tsv":pd.DataFrame(prov)}
 for n,f in outs.items():atomic_tsv(f,t/n)
 atomic_json(a.output/"software_versions_v1.json",{"python":platform.python_version(),"anndata":ad.__version__,"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"script_sha256":sha256(Path(__file__))});m=[]
 for f in sorted(x for x in a.output.rglob("*") if x.is_file() and x.name!="RUNNING.json"):m.append({"file":f.relative_to(a.output).as_posix(),"size_bytes":f.stat().st_size,"sha256":sha256(f)})
 atomic_tsv(pd.DataFrame(m),a.output/"OUTPUT_MANIFEST_v1.tsv");status={"schema_version":"phase20_spp1_vascular_immune_v1","status":"PASS_ANALYSIS_COMPLETE","completed_at":now(),"sections":14,"spots":len(outs["spot_vascular_immune_scores_v1.tsv.gz"]),"primary_rows":len(primary),"primary_relations":[x[0] for x in PRIMARY],"route_independent_primary":True,"tymp_mta_secondary_context_only":True,"section_as_biological_unit":True,"spot_p_values":False,"registry_sha256":sha256(a.analysis_registry)};atomic_json(a.output/"PHASE20_SPP1_VASCULAR_STATUS_v1.json",status);(a.output/"RUNNING.json").unlink();print(json.dumps(status,indent=2))
if __name__=="__main__":main()
