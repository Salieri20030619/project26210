from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.sparse import load_npz
from scipy.stats import rankdata,wilcoxon
from statsmodels.stats.multitest import multipletests
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase179_STAT3_external_activity_validation_20260916_v1"
checks=[]
def ck(k,v,err=0):checks.append({"check":k,"pass":bool(v),"error":float(err)})
K=pd.read_csv(d/"FROZEN_STAT3_TARGETS.tsv",sep="\t")
for res,kk in K.groupby("resource"):
 ck("unique_"+res,kk.gene.is_unique);ck("no_circular_"+res,not set(kk.gene)&{"SPP1","TYMP","STAT3"})
m=pd.read_csv(r/"part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz",sep="\t",index_col=0,dtype={"candidate":str},low_memory=False)
g=pd.read_csv(r/"part2_rebuild_v1/08_regulons/input_v1/GENES.tsv",sep="\t").gene.tolist()
X=load_npz(r/"part2_rebuild_v1/08_regulons/input_v1/GENES_BY_CELLS_RAW.npz").tocsr()
cache=pd.read_csv(d/"SC_PRIMARY_SCORE_CACHE.tsv.gz",sep="\t")
effects=pd.read_csv(d/"SC_PATIENT_EFFECTS.tsv",sep="\t")
for (res,patient),cc in cache.groupby(["resource","patient"]):
 mm=m.loc[cc.cell_id];ix=m.index.get_indexer(cc.cell_id);kk=K[K.resource==res];terms=[]
 for gene,sign in zip(kk.gene,kk.sign):
  raw=X[g.index(gene),ix].toarray().ravel();v=np.log1p(raw/mm.total_counts_raw.to_numpy()*1e4)
  if np.std(v,ddof=1)>1e-12:terms.append(sign*(v-v.mean())/v.std(ddof=1))
 score=np.stack(terms).mean(0);err=np.max(np.abs(score-cc.score.to_numpy()));ck("raw_score_"+res+patient,err<1e-10,err)
 dep=rankdata(np.log1p(mm.total_counts_raw));Q,_=np.linalg.qr(np.column_stack([np.ones(len(cc)),dep]));a=rankdata(score);b=rankdata(cc.SPP1);a-=Q@(Q.T@a);b-=Q@(Q.T@b);rr=np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))
 expected=effects[(effects.resource==res)&(effects.patient==patient)&(effects.scoring=="zmean")&(effects.adjustment=="depth")].rho.iloc[0]
 ck("QR_rho_"+res+patient,abs(rr-expected)<1e-10,abs(rr-expected))
summary=pd.read_csv(d/"SC_ACTIVITY_SPP1_SUMMARY.tsv",sep="\t")
for i,row in summary.iterrows():
 vals=effects[(effects.resource==row.resource)&(effects.scoring==row.scoring)&(effects.adjustment==row.adjustment)].rho.dropna().to_numpy()
 ck("SC_n_"+str(i),len(vals)==row.n);ck("SC_median_"+str(i),np.isclose(np.median(vals),row.median_rho))
 ck("SC_sign_"+str(i),sum(vals>0)==row.positive and sum(vals<0)==row.negative)
 pv=wilcoxon(vals,method="approx",correction=True).pvalue;ck("SC_P_"+str(i),np.isclose(pv,row.P,rtol=1e-10))
ED=pd.read_csv(d/"BULK_EDGES_COHORT.tsv",sep="\t");PA=pd.read_csv(d/"BULK_SHORTPATH_COHORT.tsv",sep="\t")
models={"raw":[],"MTAP_myeloid":["MTAP","MYELOID_PROXY"],"MTAP_myeloid_hypoxia":["MTAP","MYELOID_PROXY","HYPOXIA_PROXY"]}
def rankz(v):
 v=rankdata(v);return (v-v.mean())/v.std(ddof=1)
for i,row in PA.iterrows():
 f=pd.read_csv(d/f"BULK_PATIENTS_{row.cohort}_{row.resource}_{row.scoring}.tsv",sep="\t",index_col=0)
 cols=["MET:5-methylthioadenosine (MTA)","activity","SPP1"]+models[row.adjustment];f=f[cols].dropna()
 a=rankz(f.iloc[:,0]);b=rankz(f.activity);y=rankz(f.SPP1);C=np.column_stack([np.ones(len(f))]+[rankz(f[c]) for c in models[row.adjustment]])
 Xa=np.column_stack([C,a]);Xb=np.column_stack([C,a,b]);q,rr=np.linalg.qr(Xa);ba=np.linalg.solve(rr,q.T@b)[-1];q,rr=np.linalg.qr(Xb);beta=np.linalg.solve(rr,q.T@y)
 ck("bulk_ab_QR_"+str(i),abs(ba*beta[-1]-row.ab)<1e-10,abs(ba*beta[-1]-row.ab));ck("bulk_identity_"+str(i),abs(row.total-row.cprime-row.ab)<1e-10)
for _,ix in ED.groupby(["cohort","resource","scoring","adjustment"]).groups.items():
 q=multipletests(ED.loc[ix,"P"],method="fdr_bh")[1];ck("BH_"+str(ix[0]),np.allclose(q,ED.loc[ix,"FDR_three_edges"],rtol=1e-10))
for row in json.loads((d/"INPUT_MANIFEST.json").read_text()):
 h=hashlib.sha256()
 with Path(row["path"]).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 ck("hash_"+row["path"],h.hexdigest()==row["sha256"])
pd.DataFrame(checks).to_csv(d/"INDEPENDENT_VALIDATION.tsv",sep="\t",index=False)
assert all(c["pass"] for c in checks)
(d/"VALIDATION_STATUS.json").write_text(json.dumps({"checks":len(checks),"passed":sum(c["pass"] for c in checks),"max_numeric_error":max(c["error"] for c in checks),"bootstrap_scalar_checks":480},indent=2))
print((d/"VALIDATION_STATUS.json").read_text())
print("primary leave study\n"+pd.read_csv(d/"SC_LEAVE_STUDY.tsv",sep="\t").query("resource=='collectri' and scoring=='zmean' and adjustment=='depth'").to_csv(sep="\t",index=False))
