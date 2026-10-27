from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,h5py
from anndata.io import read_elem,sparse_dataset
from scipy.sparse import load_npz
from scipy.stats import rankdata,wilcoxon
from statsmodels.stats.multitest import multipletests
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase182_IL10_LGALS3_receiver_audit_20260916_v1"
src=Path("/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad")
mf=r/"part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz";gf=r/"part2_rebuild_v1/08_regulons/input_v1/GENES.tsv";xf=r/"part2_rebuild_v1/08_regulons/input_v1/GENES_BY_CELLS_RAW.npz";tf=r/"part2_rebuild_v1/07_trajectory/monocle3_v3/ALL_CELL_PSEUDOTIME.tsv.gz"
files=[src,mf,gf,xf,tf,d/"STAT3_TARGETS_WITHOUT_LGALS3.tsv"];panel=["IL10","IL10RA","IL10RB","IL6","IL6R","IL6ST","SOCS3","ADORA2B","LGALS3","SPP1"]
m=pd.read_csv(mf,sep="\t",index_col=0,dtype={"candidate":str},low_memory=False)
with h5py.File(src,"r") as h:
 obs=read_elem(h["obs"]);var=read_elem(h["var"]);assert len(obs)==258120 and obs.sample_id.nunique()==75
 assert set(panel)<=set(var.index)
 xp=sparse_dataset(h["layers/counts"]);pi=var.index.get_indexer(panel);raw=np.zeros((len(obs),len(panel)),dtype=np.float64)
 for start in range(0,len(obs),5000):
  block=xp[start:start+5000,:].tocsr();raw[start:start+len(block.indptr)-1]=block[:,pi].toarray()
obs["source_class"]=obs.cell_type_broad_v1.astype(str)
obs.loc[obs.source_class.isin(["Macrophage","Monocyte"]),"source_class"]+="_parent_not_in_refined"
obs.loc[obs.malignant_conservative_v1.astype(bool),"source_class"]="Malignant_conservative"
for cid,mm in m.groupby("candidate"):
 label="C4" if cid=="4" else "Monocyte_refined" if cid in ["1","Monocyte"] else "OtherMacrophage_refined"
 assert mm.index.isin(obs.index).all();obs.loc[mm.index,"source_class"]=label
small=obs[["study","sample_id","curated_patient_id","source_tissue","source_class","total_counts_raw"]].copy()
small.to_csv(d/"SOURCE_CELL_METADATA.tsv.gz",sep="\t",index_label="cell_id",compression="gzip")
np.savez_compressed(d/"SOURCE_PANEL_CACHE.npz",counts=raw,genes=np.array(panel),cell_ids=obs.index.to_numpy(str))
tum=obs.source_tissue.astype(str).str.lower()=="tumor"
out=[]
for (pat,cl),o in obs[tum].groupby(["curated_patient_id","source_class"],observed=True):
 ix=obs.index.get_indexer(o.index);ct=raw[ix];tot=o.total_counts_raw.sum()
 patient_idx=obs.index.get_indexer(obs[tum&(obs.curated_patient_id==pat)].index)
 for j,g in enumerate(panel):
  denominator=raw[patient_idx,j].sum()
  out.append(dict(patient=pat,study=str(o.study.iloc[0]),source_class=cl,gene=g,cells=len(o),eligible=len(o)>=20,detected=int((ct[:,j]>0).sum()),detection=float((ct[:,j]>0).mean()),counts=float(ct[:,j].sum()),log2CPM=float(np.log2(ct[:,j].sum()/tot*1e6+1)),observed_RNA_share=float(ct[:,j].sum()/denominator) if denominator else np.nan))
T=pd.DataFrame(out);T.to_csv(d/"SOURCE_PATIENT_PANEL.tsv",sep="\t",index=False)
s=[]
for (cl,g),p in T[T.eligible].groupby(["source_class","gene"]):
 s.append(dict(source_class=cl,gene=g,n_patients=len(p),patients_detected=int((p.detected>0).sum()),patients_ge10_detected=int((p.detected>=10).sum()),median_detection=p.detection.median(),median_log2CPM=p.log2CPM.median(),median_observed_RNA_share=p.observed_RNA_share.median()))
SS=pd.DataFrame(s);SS.to_csv(d/"SOURCE_PANEL_SUMMARY.tsv",sep="\t",index=False)
print("SOURCE_COMPLETE",len(T),flush=True)
# Refined receiver coverage; total macrophage pool includes C4 and other mature macrophages.
mac=m[(~m.candidate.isin(["1","Monocyte"]))&(m.source_tissue.astype(str).str.lower()=="tumor")]
cov=[]
for pat,mm in mac.groupby("curated_patient_id"):
 ix=obs.index.get_indexer(mm.index);rr=raw[ix]
 if len(mm)<20:continue
 for g in panel:
  a=rr[:,panel.index(g)]>0;cov.append(dict(patient=pat,feature=g,cells=len(mm),detected=int(a.sum()),fraction=a.mean(),association_eligible=int(a.sum())>=10))
 for label,a,b in [("IL10RA_AND_IL10RB","IL10RA","IL10RB"),("IL6R_AND_IL6ST","IL6R","IL6ST")]:
  hit=(rr[:,panel.index(a)]>0)&(rr[:,panel.index(b)]>0);cov.append(dict(patient=pat,feature=label,cells=len(mm),detected=int(hit.sum()),fraction=hit.mean(),association_eligible=int(hit.sum())>=10))
CV=pd.DataFrame(cov);CV.to_csv(d/"RECEIVER_PATIENT_COVERAGE.tsv",sep="\t",index=False)
cs=[]
for g,p in CV.groupby("feature"):
 cs.append(dict(feature=g,n_patients=len(p),patients_any=int((p.detected>0).sum()),patients_ge10_cells=int(p.association_eligible.sum()),median_detection=p.fraction.median(),RNA_gate_pass=int(p.association_eligible.sum())>=10))
CS=pd.DataFrame(cs);CS.to_csv(d/"RECEIVER_COVERAGE_SUMMARY.tsv",sep="\t",index=False)
print("RECEIVER\n"+CS.to_csv(sep="\t",index=False),flush=True)
# Fixed scores and three predefined edges.
K=pd.read_csv(d/"STAT3_TARGETS_WITHOUT_LGALS3.tsv",sep="\t");genes=pd.read_csv(gf,sep="\t").gene.tolist();X=load_npz(xf).tocsr();target=sorted(set(K.gene));assert not set(target)&{"LGALS3","SPP1","TYMP","STAT3"}
E=np.log1p(X[[genes.index(g) for g in target]].toarray().T/m.total_counts_raw.to_numpy()[:,None]*1e4)
rr=raw[obs.index.get_indexer(m.index)];genexpr=np.log1p(rr/m.total_counts_raw.to_numpy()[:,None]*1e4)
pt=pd.read_csv(tf,sep="\t",index_col=0,dtype={"partition":str},low_memory=False);pt=pt[pt["mode"]=="primary"]
def score(A,sign,kind):
 if kind=="rankmean":A=rankdata(A,axis=0)
 sd=A.std(0,ddof=1);ok=sd>1e-12
 return np.round(((A[:,ok]-A[:,ok].mean(0))/sd[ok]*sign[ok]).mean(1),12),int(ok.sum())
checks=[]
def partial(x,y,cv):
 C=np.column_stack([np.ones(len(x))]+[rankdata(a) for a in cv]);rx=rankdata(x);ry=rankdata(y)
 a=rx-C@np.linalg.lstsq(C,rx,rcond=None)[0];b=ry-C@np.linalg.lstsq(C,ry,rcond=None)[0]
 if np.std(a)<1e-12 or np.std(b)<1e-12:return np.nan
 rho=np.corrcoef(a,b)[0,1]
 if np.linalg.matrix_rank(C)==C.shape[1]:
  Q,_=np.linalg.qr(C);aa=rx-Q@(Q.T@rx);bb=ry-Q@(Q.T@ry);checks.append(dict(check="QR_partial",passed=abs(rho-np.corrcoef(aa,bb)[0,1])<1e-10))
 return rho
effects=[];elig=[];paired=[];cache={}
for pat,mm in mac.groupby("curated_patient_id"):
 ix=m.index.get_indexer(mm.index);l=genexpr[ix,panel.index("LGALS3")];s=genexpr[ix,panel.index("SPP1")];dep=np.log1p(mm.total_counts_raw.to_numpy());ng=mm.n_genes_raw.to_numpy();c4=(mm.candidate=="4").to_numpy(float)
 if sum(c4==1)>=20 and sum(c4==0)>=20:
  lc=rr[ix,panel.index("LGALS3")];tot=mm.total_counts_raw.to_numpy();delta=np.log2(lc[c4==1].sum()/tot[c4==1].sum()*1e6+1)-np.log2(lc[c4==0].sum()/tot[c4==0].sum()*1e6+1)
  paired.append(dict(patient=pat,study=str(mm.study.iloc[0]),delta=delta,C4_cells=int(sum(c4==1)),other_cells=int(sum(c4==0))))
 if len(mm)<50:continue
 for res,kk in K.groupby("resource"):
  ixg=[target.index(g) for g in kk.gene]
  for scoring in ["zmean","rankmean"]:
   st,nvar=score(E[ix][:,ixg],kk.sign.to_numpy(),scoring)
   for edge,x,y,detok in [("LGALS3_SPP1",l,s,(l>0).sum()>=10 and (s>0).sum()>=10),("STAT3_LGALS3",st,l,(l>0).sum()>=10),("STAT3_SPP1",st,s,(s>0).sum()>=10)]:
    ok=detok and nvar>=(20 if res=="collectri" else 10)
    elig.append(dict(patient=pat,resource=res,scoring=scoring,edge=edge,cells=len(mm),LGALS3_detected=int((l>0).sum()),SPP1_detected=int((s>0).sum()),variable_targets=nvar,eligible=ok))
    if not ok:continue
    base=dict(patient=pat,study=str(mm.study.iloc[0]),resource=res,scoring=scoring,edge=edge)
    for adj,cv in [("raw",[]),("depth",[dep]),("depth_genes",[dep,ng]),("depth_C4",[dep,c4])]:
     effects.append({**base,"adjustment":adj,"rho":partial(x,y,cv)})
    q=pt.reindex(mm.index);vals=[]
    for part in q.partition.dropna().unique():
     ii=(q.partition==part).to_numpy()&np.isfinite(q.pt_primary.to_numpy())
     if sum(ii)<50 or q.pt_primary.to_numpy()[ii].size<10 or len(np.unique(q.pt_primary.to_numpy()[ii]))<10:continue
     if (y[ii]>0).sum()<10:continue
     if edge=="LGALS3_SPP1" and (x[ii]>0).sum()<10:continue
     vals.append(partial(x[ii],y[ii],[dep[ii],q.pt_primary.to_numpy()[ii]]))
    if vals:effects.append({**base,"adjustment":"depth_pseudotime","rho":float(np.nanmedian(vals))})
   cache["__".join([res,scoring,pat])]=np.column_stack([st,l,s,dep,ng,c4])
F=pd.DataFrame(effects);F.to_csv(d/"PATIENT_EDGE_EFFECTS.tsv",sep="\t",index=False)
pd.DataFrame(elig).to_csv(d/"EDGE_ELIGIBILITY.tsv",sep="\t",index=False);pd.DataFrame(paired).to_csv(d/"C4_LGALS3_PATIENTS.tsv",sep="\t",index=False);np.savez_compressed(d/"SCORE_EDGE_CACHE.npz",**cache)
su=[];los=[]
for key,p in F.groupby(["resource","scoring","adjustment","edge"]):
 v=p.rho.dropna().to_numpy();n=len(v);rng=np.random.default_rng(182);bs=np.median(v[rng.integers(0,n,(5000,n))],axis=1)
 su.append(dict(zip(["resource","scoring","adjustment","edge"],key))|dict(n=n,median_rho=np.median(v),positive=int(sum(v>0)),negative=int(sum(v<0)),P=wilcoxon(v,method="approx",correction=True).pvalue if n>=10 else np.nan,CI_low=np.quantile(bs,.025),CI_high=np.quantile(bs,.975)))
 for st in p.study.unique():
  z=p[p.study!=st].rho.dropna();los.append(dict(zip(["resource","scoring","adjustment","edge"],key))|dict(omitted_study=st,n=len(z),median_rho=z.median(),positive=int(sum(z>0))))
SU=pd.DataFrame(su)
for _,ix in SU.groupby(["resource","scoring","adjustment"]).groups.items():
 good=SU.loc[ix,"P"].dropna();SU.loc[good.index,"FDR_three_edges"]=multipletests(good,method="fdr_bh")[1]
SU.to_csv(d/"EDGE_SUMMARY.tsv",sep="\t",index=False);pd.DataFrame(los).to_csv(d/"LEAVE_STUDY.tsv",sep="\t",index=False)
v=pd.DataFrame(paired).delta.to_numpy();cc=dict(n=len(v),median_log2CPM_delta=float(np.median(v)),positive=int(sum(v>0)),negative=int(sum(v<0)),P=float(wilcoxon(v,method="approx",correction=True).pvalue))
(d/"C4_LGALS3_SUMMARY.json").write_text(json.dumps(cc,indent=2))
pd.DataFrame(checks).to_csv(d/"NUMERIC_VALIDATION.tsv",sep="\t",index=False);assert all(x["passed"] for x in checks)
manifest=[]
for file in files:
 hh=hashlib.sha256()
 with file.open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):hh.update(b)
 manifest.append(dict(path=str(file),bytes=file.stat().st_size,sha256=hh.hexdigest()))
(d/"INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2));assert manifest[0]["sha256"]=="e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7"
(d/"COMPUTE_STATUS.json").write_text(json.dumps(dict(status="COMPUTED",checks=len(checks)),indent=2))
print("EDGES\n"+SU[(SU.resource=="collectri")&(SU.scoring=="zmean")].to_csv(sep="\t",index=False),flush=True)
print("C4",cc,flush=True)
print("IL10_SOURCES\n"+SS[SS.gene=="IL10"].sort_values("median_detection",ascending=False).to_csv(sep="\t",index=False),flush=True)
