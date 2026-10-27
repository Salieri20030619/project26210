#!/usr/bin/env python3
from __future__ import annotations
import argparse,gzip,hashlib,json,math,os,platform
from datetime import datetime
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd,scanpy as sc,scipy
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.stats import binomtest,spearmanr,wilcoxon

SEED=20260814; BOOTSTRAPS=9999
STATES=["Monocyte","Other Macrophage","C1QC-state Macrophage","SPP1-state Macrophage"]
PRIMARY=[("monocyte_program",-1),("spp1_state_no_spp1",1),("spp1_vs_c1qc",1)]
SUPPORTING=["spp1_state_full","vascular_remodeling","myeloid_immunosuppression","polyamine_mta_context"]
INPUT_H5AD_SHA="e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"
ASSIGNMENTS_SHA="f8d0d01dcb76b3a19055accae42be3a2e14c6ba3a9da85db30c9592f2e635e6d"
def now():return datetime.now().astimezone().isoformat(timespec="seconds")
def sha256(p):
 d=hashlib.sha256()
 with Path(p).open("rb") as h:
  for b in iter(lambda:h.read(8*1024*1024),b""):d.update(b)
 return d.hexdigest()
def z(x):
 x=np.asarray(x,float);s=np.nanstd(x);return np.zeros(len(x)) if not np.isfinite(s) or s<=0 else (x-np.nanmean(x))/s
def rho(x,y):
 x,y=np.asarray(x,float),np.asarray(y,float);m=np.isfinite(x)&np.isfinite(y)
 return math.nan if m.sum()<5 or np.ptp(x[m])==0 or np.ptp(y[m])==0 else float(spearmanr(x[m],y[m]).statistic)
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);m=np.isfinite(p)
 if not m.any():return out
 v=p[m];o=np.argsort(v);q=np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];r=np.empty(len(v));r[o]=np.minimum(q,1);out[m]=r;return out
def atomic_tsv(f,p):
 q=Path(str(p)+".part")
 if str(p).endswith(".gz"):
  with gzip.open(q,"wt",encoding="utf-8",newline="") as h:f.to_csv(h,sep="\t",index=False,lineterminator="\n")
 else:f.to_csv(q,sep="\t",index=False,lineterminator="\n")
 os.replace(q,p)
def atomic_json(p,x):
 q=Path(str(p)+".part");q.write_text(json.dumps(x,indent=2)+"\n",encoding="utf-8");os.replace(q,p)
def score(obj,genes,registry):
 ids=[x[0] for x in PRIMARY]+SUPPORTING+["c1qc_state"]
 wanted=genes[genes.signature_id.isin(ids)].gene.drop_duplicates().tolist();present=[g for g in wanted if g in obj.var_names];idx=[obj.var_names.get_loc(g) for g in present];x=obj.layers["counts"][:,idx]
 if sparse.issparse(x):x=x.toarray()
 x=np.asarray(x,float);tot=np.asarray(obj.layers["counts"].sum(axis=1)).ravel().astype(float)
 if not np.isfinite(tot).all() or (tot<=0).any():raise RuntimeError("invalid library sizes")
 ln=np.log1p(x*(10000/tot[:,None]));pos={g:i for i,g in enumerate(present)};out={};qc=[]
 limits=registry.set_index("signature_id").minimum_present_genes.astype(int)
 for sid in ids:
  d=genes[genes.signature_id==sid].sort_values("gene_order");use=[g for g in d.gene if g in pos];minimum=int(limits[sid])
  if len(use)<minimum:raise RuntimeError(f"coverage {sid} {len(use)}/{minimum}")
  out[sid]=np.column_stack([z(ln[:,pos[g]]) for g in use]).mean(axis=1);qc.append({"signature_id":sid,"defined":len(d),"minimum_present":minimum,"present":len(use),"missing":";".join(g for g in d.gene if g not in pos)})
 out["spp1_vs_c1qc"]=z(out["spp1_state_no_spp1"])-z(out["c1qc_state"])
 return pd.DataFrame(out,index=obj.obs_names),pd.DataFrame(qc)
def main():
 p=argparse.ArgumentParser();p.add_argument("--h5ad",type=Path,required=True);p.add_argument("--assignments",type=Path,required=True);p.add_argument("--signature-registry",type=Path,required=True);p.add_argument("--signature-genes",type=Path,required=True);p.add_argument("--analysis-registry",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
 if a.output.exists():raise FileExistsError(a.output)
 a.output.mkdir(parents=True);t=a.output/"tables";t.mkdir();atomic_json(a.output/"RUNNING.json",{"status":"RUNNING","started_at":now()})
 hsha=sha256(a.h5ad);asha=sha256(a.assignments)
 if hsha!=INPUT_H5AD_SHA or asha!=ASSIGNMENTS_SHA:raise RuntimeError("frozen input hash mismatch")
 ass=pd.read_csv(a.assignments,sep="\t");
 if len(ass)!=229032 or ass.exact_cell_id.nunique()!=len(ass):raise RuntimeError("assignment identity")
 sel=ass[ass.state_v2.isin(STATES)].copy();reg=pd.read_csv(a.signature_registry,sep="\t");genes=pd.read_csv(a.signature_genes,sep="\t")
 backed=ad.read_h5ad(a.h5ad,backed="r")
 if "counts" not in backed.layers or "X_scVI" not in backed.obsm or not backed.obs_names.is_unique:raise RuntimeError("H5AD structure")
 idx=backed.obs_names.get_indexer(sel.exact_cell_id.astype(str));
 if (idx<0).any():raise RuntimeError("assignment cells missing")
 obj=backed[idx].to_memory();backed.file.close();obj.obs=sel.set_index("exact_cell_id").loc[obj.obs_names].copy()
 cells=[];patient_rows=[];state_rows=[];graph_rows=[];qc_rows=[];roots=[]
 for study in sorted(obj.obs.study_v1.unique()):
  print(f"[{now()}] {study}",flush=True);sub=obj[obj.obs.study_v1==study].copy();counts=sub.obs.state_v2.value_counts();
  if sub.n_obs<200 or counts.get("Monocyte",0)<50 or counts.get("SPP1-state Macrophage",0)<50:raise RuntimeError(f"study coverage {study}")
  sc.pp.neighbors(sub,n_neighbors=30,n_pcs=30,use_rep="X_scVI",metric="euclidean",random_state=SEED);sc.tl.diffmap(sub,n_comps=20,random_state=SEED)
  labels_n,comp=connected_components(sub.obsp["connectivities"],directed=False);latent=np.asarray(sub.obsm["X_scVI"],float);mono=np.flatnonzero(sub.obs.state_v2.astype(str).eq("Monocyte").to_numpy());center=np.median(latent[mono],axis=0);root=int(mono[np.argmin(((latent[mono]-center)**2).sum(axis=1))]);root_comp=int(comp[root]);root_fraction=float((comp==root_comp).mean())
  if root_fraction<.95:raise RuntimeError(f"root component {study} {root_fraction}")
  sub.uns["iroot"]=root;sc.tl.dpt(sub,n_dcs=15);pt=np.asarray(sub.obs.dpt_pseudotime,float)
  if not np.isfinite(pt).all():raise RuntimeError(f"nonfinite DPT {study}")
  scores,qc=score(sub,genes,reg);qc.insert(0,"study",study);qc_rows.append(qc);frame=sub.obs[["patient","study_v1","state_v2"]].copy();frame.insert(0,"cell_id",sub.obs_names);frame["dpt_pseudotime"]=pt
  for c in scores.columns:frame[c]=scores[c].to_numpy()
  cells.append(frame);roots.append({"study":study,"root_cell":str(sub.obs_names[root]),"root_state":str(sub.obs.state_v2.iloc[root]),"n_cells":sub.n_obs,"n_components":labels_n,"root_component_fraction":root_fraction})
  for state,g in frame.groupby("state_v2",observed=True):state_rows.append({"study":study,"state":state,"n_cells":len(g),"n_patients":g.patient.nunique(),"median_pseudotime":float(g.dpt_pseudotime.median()),"median_monocyte_program":float(g.monocyte_program.median()),"median_spp1_state_no_spp1":float(g.spp1_state_no_spp1.median()),"median_spp1_vs_c1qc":float(g.spp1_vs_c1qc.median())})
  conn=sub.obsp["connectivities"].tocsr();labs=sub.obs.state_v2.astype(str).to_numpy()
  for source in STATES:
   sm=labs==source;den=float(conn[sm].sum())
   for target in STATES:graph_rows.append({"study":study,"source":source,"target":target,"edge_weight":float(conn[sm][:,labs==target].sum()),"source_normalized_fraction":float(conn[sm][:,labs==target].sum()/den) if den else math.nan})
  endpoints=[x[0] for x in PRIMARY]+SUPPORTING
  for patient,g in frame.groupby("patient",sort=True):
   eligible=len(g)>=50 and g.state_v2.nunique()>=2 and np.ptp(g.dpt_pseudotime)>0
   for endpoint in endpoints:patient_rows.append({"study":study,"patient":patient,"endpoint":endpoint,"role":"primary" if endpoint in dict(PRIMARY) else "supporting","expected_direction":dict(PRIMARY).get(endpoint,0),"n_cells":len(g),"n_states":g.state_v2.nunique(),"eligible":eligible and np.ptp(g[endpoint])>0,"spearman_rho":rho(g.dpt_pseudotime,g[endpoint]) if eligible else math.nan})
 patient=pd.DataFrame(patient_rows);rng=np.random.default_rng(SEED);agg=[];loo=[]
 for study in sorted(patient.study.unique()):
  for endpoint,direction in PRIMARY:
   f=patient[(patient.study==study)&(patient.endpoint==endpoint)&patient.eligible&patient.spearman_rho.notna()];v=f.spearman_rho.to_numpy();n=len(v)
   if n<5:agg.append({"study":study,"endpoint":endpoint,"expected_direction":direction,"eligible_patients":n,"median_rho":math.nan,"expected_direction_patients":int((v*direction>0).sum()),"opposite_direction_patients":int((v*direction<0).sum()),"wilcoxon_p":math.nan,"sign_test_p":math.nan,"bootstrap_ci_low":math.nan,"bootstrap_ci_high":math.nan,"loo_expected_direction_fraction":math.nan,"decision":"NOT_TESTED_LT5_PATIENTS"});continue
   try:wp=float(wilcoxon(v).pvalue)
   except ValueError:wp=1.0
   boot=np.median(v[rng.integers(0,n,size=(BOOTSTRAPS,n))],axis=1);lv=[]
   for left in f.patient:
    med=float(np.median(f.loc[f.patient!=left,"spearman_rho"]));lv.append(med*direction);loo.append({"study":study,"endpoint":endpoint,"left_out_patient":left,"median_rho":med,"expected_direction":direction})
   agg.append({"study":study,"endpoint":endpoint,"expected_direction":direction,"eligible_patients":n,"median_rho":float(np.median(v)),"expected_direction_patients":int((v*direction>0).sum()),"opposite_direction_patients":int((v*direction<0).sum()),"wilcoxon_p":wp,"sign_test_p":float(binomtest(int((v*direction>0).sum()),n,.5).pvalue),"bootstrap_ci_low":float(np.quantile(boot,.025)),"bootstrap_ci_high":float(np.quantile(boot,.975)),"loo_expected_direction_fraction":float((np.asarray(lv)>0).mean()),"decision":"PENDING"})
 aggregate=pd.DataFrame(agg);aggregate["wilcoxon_fdr_bh_within_study"]=aggregate.groupby("study")["wilcoxon_p"].transform(lambda x:bh(x.to_numpy()))
 test=aggregate.decision.eq("PENDING");need=np.ceil(aggregate.eligible_patients*2/3);support=(aggregate.median_rho*aggregate.expected_direction>=.1)&(aggregate.expected_direction_patients>=need)&(aggregate.wilcoxon_fdr_bh_within_study<.05)&(aggregate.bootstrap_ci_low*aggregate.expected_direction>0)&(aggregate.bootstrap_ci_high*aggregate.expected_direction>0)&(aggregate.loo_expected_direction_fraction==1)
 aggregate.loc[test,"decision"]=np.where(support[test],"STUDY_CONTINUITY_SUPPORT",np.where(aggregate.loc[test,"expected_direction_patients"]>=need[test],"DIRECTIONAL_NOT_FDR",np.where((aggregate.loc[test,"expected_direction_patients"]>=2)&(aggregate.loc[test,"opposite_direction_patients"]>=2),"HETEROGENEOUS","NOT_SUPPORTED")))
 synth=[]
 for endpoint,direction in PRIMARY:
  f=aggregate[(aggregate.endpoint==endpoint)&aggregate.median_rho.notna()];synth.append({"endpoint":endpoint,"expected_direction":direction,"studies_tested":len(f),"expected_direction_studies":int((f.median_rho*direction>0).sum()),"median_study_median_rho":float(f.median_rho.median()) if len(f) else math.nan,"studies_with_support":int(f.decision.eq("STUDY_CONTINUITY_SUPPORT").sum()),"interpretation":"DESCRIPTIVE_DIRECTION_SYNTHESIS_ONLY"})
 outputs={"myeloid_trajectory_cell_scores_v1.tsv.gz":pd.concat(cells,ignore_index=True),"patient_pseudotime_associations_v1.tsv":patient,"study_primary_inference_v1.tsv":aggregate,"leave_one_patient_out_v1.tsv":pd.DataFrame(loo),"cross_study_direction_synthesis_v1.tsv":pd.DataFrame(synth),"state_trajectory_summary_v1.tsv":pd.DataFrame(state_rows),"graph_state_connectivity_v1.tsv":pd.DataFrame(graph_rows),"signature_gene_coverage_v1.tsv":pd.concat(qc_rows,ignore_index=True),"trajectory_root_qc_v1.tsv":pd.DataFrame(roots)}
 for n,f in outputs.items():atomic_tsv(f,t/n)
 atomic_json(a.output/"software_versions_v1.json",{"python":platform.python_version(),"anndata":ad.__version__,"scanpy":sc.__version__,"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"script_sha256":sha256(Path(__file__))});manifest=[]
 for f in sorted(x for x in a.output.rglob("*") if x.is_file() and x.name!="RUNNING.json"):manifest.append({"file":f.relative_to(a.output).as_posix(),"size_bytes":f.stat().st_size,"sha256":sha256(f)})
 atomic_tsv(pd.DataFrame(manifest),a.output/"OUTPUT_MANIFEST_v1.tsv");status={"schema_version":"phase21_monocyte_spp1_trajectory_v1","status":"PASS_ANALYSIS_COMPLETE","completed_at":now(),"studies":5,"patients":int(sel.patient.nunique()),"cells":len(sel),"study_separate_graphs":True,"primary_endpoints":[x[0] for x in PRIMARY],"patient_as_biological_unit":True,"cross_study_cell_pooling":False,"lineage_conversion_claim":False,"input_h5ad_sha256":hsha,"assignments_sha256":asha,"registry_sha256":sha256(a.analysis_registry)};atomic_json(a.output/"PHASE21_MONOCYTE_SPP1_TRAJECTORY_STATUS_v1.json",status);(a.output/"RUNNING.json").unlink();print(json.dumps(status,indent=2))
if __name__=="__main__":main()
