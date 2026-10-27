from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,h5py
from scipy.stats import rankdata,wilcoxon
from statsmodels.stats.multitest import multipletests
from anndata.io import sparse_dataset,read_elem
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase182_IL10_LGALS3_receiver_audit_20260916_v1"
checks=[]
def ck(n,v,e=0):checks.append(dict(check=n,passed=bool(v),error=float(e)))
cache=np.load(d/"SCORE_EDGE_CACHE.npz");F=pd.read_csv(d/"PATIENT_EDGE_EFFECTS.tsv",sep="\t")
for key in cache.files:
 res,sc,pat=key.split("__",2);a=cache[key];cols={"LGALS3_SPP1":(1,2),"STAT3_LGALS3":(0,1),"STAT3_SPP1":(0,2)}
 for edge,(x,y) in cols.items():
  row=F[(F.resource==res)&(F.scoring==sc)&(F.patient==pat)&(F.edge==edge)&(F.adjustment=="depth")]
  if row.empty:continue
  C=np.column_stack([np.ones(len(a)),rankdata(a[:,3])]);Q,_=np.linalg.qr(C);xx=rankdata(a[:,x]);yy=rankdata(a[:,y]);xx-=Q@(Q.T@xx);yy-=Q@(Q.T@yy)
  rr=np.dot(xx,yy)/np.sqrt(np.dot(xx,xx)*np.dot(yy,yy));err=abs(rr-row.rho.iloc[0]);ck("cache_"+key+edge,err<1e-10,err)
SU=pd.read_csv(d/"EDGE_SUMMARY.tsv",sep="\t")
for i,row in SU.iterrows():
 p=F[(F.resource==row.resource)&(F.scoring==row.scoring)&(F.adjustment==row.adjustment)&(F.edge==row.edge)].rho.dropna().to_numpy()
 ck("summary_"+str(i),len(p)==row.n and np.isclose(np.median(p),row.median_rho) and sum(p>0)==row.positive and sum(p<0)==row.negative)
for _,ix in SU.groupby(["resource","scoring","adjustment"]).groups.items():
 pv=SU.loc[ix,"P"].dropna();q=multipletests(pv,method="fdr_bh")[1];ck("FDR_"+str(ix[0]),np.allclose(q,SU.loc[pv.index,"FDR_three_edges"]))
# Replicate original STAT3-SPP1 patient effects: no target change was required.
old=pd.read_csv(r/"phase179_STAT3_external_activity_validation_20260916_v1/run_v2/SC_PATIENT_EFFECTS.tsv",sep="\t")
new=F[F.edge=="STAT3_SPP1"].merge(old,on=["patient","study","resource","scoring","adjustment"],suffixes=("_new","_old"),validate="one_to_one")
ck("old_edge_row_count",len(new)==len(old))
ck("old_edge_unchanged",np.allclose(new.rho_new,new.rho_old,atol=1e-12))
# Raw panel spot checks directly against immutable H5AD.
z=np.load(d/"SOURCE_PANEL_CACHE.npz");raw=z["counts"];genes=z["genes"].tolist()
src=Path(json.loads((d/"INPUT_MANIFEST.json").read_text())[0]["path"])
rng=np.random.default_rng(262182);ix=np.sort(rng.choice(len(raw),80,replace=False))
with h5py.File(src,"r") as h:
 var=read_elem(h["var"]);xx=sparse_dataset(h["layers/counts"])[ix,:].tocsr()[:,var.index.get_indexer(genes)].toarray()
 for j,i in enumerate(ix):ck("H5AD_panel_"+str(i),np.array_equal(xx[j],raw[i]))
obs=pd.read_csv(d/"SOURCE_CELL_METADATA.tsv.gz",sep="\t",index_col=0);T=pd.read_csv(d/"SOURCE_PATIENT_PANEL.tsv",sep="\t")
for (pat,g),p in T.groupby(["patient","gene"]):
 ii=(obs.curated_patient_id==pat)&(obs.source_tissue.str.lower()=="tumor");total=raw[ii.to_numpy(),genes.index(g)].sum()
 ck("RNA_total_"+pat+g,np.isclose(p.counts.sum(),total))
 if total>0:ck("RNA_share_"+pat+g,abs(p.observed_RNA_share.sum()-1)<1e-10)
for item in json.loads((d/"INPUT_MANIFEST.json").read_text()):
 h=hashlib.sha256()
 with Path(item["path"]).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 ck("hash_"+item["path"],h.hexdigest()==item["sha256"])
pd.DataFrame(checks).to_csv(d/"INDEPENDENT_VALIDATION.tsv",sep="\t",index=False);assert all(c["passed"] for c in checks)
numeric=pd.read_csv(d/"NUMERIC_VALIDATION.tsv",sep="\t");assert numeric.passed.all()
status=dict(independent_checks=len(checks),numeric_checks=len(numeric),total_checks=len(checks)+len(numeric),all_pass=True,max_error=max(c["error"] for c in checks))
(d/"VALIDATION_STATUS.json").write_text(json.dumps(status,indent=2))
print(json.dumps(status));print("SENSITIVITY\n"+SU[(SU.adjustment=="depth")].to_csv(sep="\t",index=False))
print("LOSO\n"+pd.read_csv(d/"LEAVE_STUDY.tsv",sep="\t").query("resource=='collectri' and scoring=='zmean' and adjustment=='depth'").to_csv(sep="\t",index=False))
