from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime,shutil
from scipy.stats import rankdata,norm
r=Path("/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1");d=r/"phase178_MTA_SPP1_candidate_mechanism_20260916_v1"
pub=Path("/mnt/BioDisk/wangmingyuan/project26210_phase98_public_paired_metabolomics_v1")
models={"raw":[],"MTAP_myeloid":["MTAP","MYELOID_PROXY"],"MTAP_myeloid_hypoxia":["MTAP","MYELOID_PROXY","HYPOXIA_PROXY"]}
B=2000;cohorts=["RC18","RC20","CPTAC","CPTAC_val"]
plan={"primary":"MTA x ADORA2B interaction coefficient on bulk SPP1; MTAP_myeloid model","hypothesis":"Higher receptor RNA background associated with stronger positive MTA-SPP1 relationship","models":models,"bootstrap":B,"unit":"independent patients within cohort","FDR":"single prespecified primary pooled interaction; no FDR. Raw and hypoxia are sensitivity, not additional selected primary tests","not_mediation":True,"not_activity_measurement":True,"cell_source":"bulk; not macrophage-specific","created_at":datetime.datetime.now().isoformat()}
(d/"RECEPTOR_INTERACTION_PLAN_BEFORE_FIT.json").write_text(json.dumps(plan,indent=2))
mapping=pd.read_csv(pub/"MasterMapping_ccRCC.csv",dtype=str).fillna("")
files=[pub/"MasterMapping_ccRCC.csv"];rows=[];draws={};checks=[]
def z(v):v=rankdata(v);return (v-v.mean())/v.std(ddof=1)
def design(f,cov):
 a=z(f["MET:5-methylthioadenosine (MTA)"]);b=z(f.ADORA2B)
 return np.column_stack([np.ones(len(f)),a,b,a*b]+[z(f[c]) for c in cov]),z(f.SPP1)
def fit(X,y):
 beta=np.linalg.lstsq(X,y,rcond=None)[0];inv=np.linalg.inv(X.T@X)
 residual=y-X@beta;hat=np.einsum("ni,ij,nj->n",X,inv,X);u=residual/(1-hat)
 v=inv@((X*u[:,None]).T@(X*u[:,None]))@inv;se=np.sqrt(np.diag(v))
 return beta,se
for ci,c in enumerate(cohorts):
 fp=r/"phase162_MTA_polyamine_part1_reassessment_20260915_v1"/("PUBLIC_PATIENTS_"+c+".tsv")
 rp=pub/"downloads_v1"/("matched_tpm_"+c+".csv");files += [fp,rp]
 f=pd.read_csv(fp,sep="\t",index_col=0)
 h=pd.read_csv(rp,nrows=0).columns[0];rna=pd.read_csv(rp,usecols=[h,"ADORA2B"],index_col=0)
 m=mapping[(mapping.Dataset==c)&(mapping.TN=="Tumor")].copy()
 m["patient"]=m.CommonID.where(m.CommonID!="",m.RNAID.str.replace(r"-T$","",regex=True)).str.replace(r"^RC(?:18|20):","",regex=True)
 assert not m.patient.duplicated().any()
 v=np.log2(rna.loc[m.RNAID,"ADORA2B"].to_numpy()+1);g=pd.DataFrame({"ADORA2B":v},index=m.patient)
 f=f.join(g,validate="one_to_one");assert f.ADORA2B.notna().all()
 f.to_csv(d/("RECEPTOR_PATIENTS_"+c+".tsv"),sep="\t")
 for ai,(adj,cov) in enumerate(models.items()):
  ff=f[["MET:5-methylthioadenosine (MTA)","ADORA2B","SPP1"]+cov].dropna()
  X,y=design(ff,cov);beta,se=fit(X,y);assert np.linalg.matrix_rank(X)==X.shape[1]
  rng=np.random.default_rng(262178+ci*100+ai)
  ii=rng.integers(0,len(ff),(B,len(ff)))
  arrays=[]
  for col in ff.columns:
   a=rankdata(ff[col].to_numpy()[ii],axis=1);a=(a-a.mean(1)[:,None])/a.std(1,ddof=1)[:,None];arrays.append(a)
  a,b,yb=arrays[:3];xb=np.stack([np.ones_like(a),a,b,a*b]+arrays[3:],axis=2)
  gram=np.einsum("bni,bnj->bij",xb,xb);rhs=np.einsum("bni,bn->bi",xb,yb)
  bb=np.einsum("bij,bj->bi",np.linalg.pinv(gram),rhs)[:,3];assert np.isfinite(bb).all()
  for j in range(10):
   xx,yy=design(ff.iloc[ii[j]],cov);direct=np.linalg.lstsq(xx,yy,rcond=None)[0][3]
   checks.append({"check":f"bootstrap_{c}_{adj}_{j}","pass":bool(np.isclose(bb[j],direct,atol=1e-10))})
  q,rr=np.linalg.qr(X);checkbeta=np.linalg.solve(rr,q.T@y)
  checks.append({"check":f"QR_{c}_{adj}","pass":bool(np.allclose(beta,checkbeta,atol=1e-10))})
  rows.append({"cohort":c,"adjustment":adj,"n":len(ff),"MTA_beta_at_mean_receptor":beta[1],"receptor_beta":beta[2],"interaction_beta":beta[3],"HC3_SE":se[3],"P_HC3_normal":2*norm.sf(abs(beta[3]/se[3])),"bootstrap_low":np.quantile(bb,.025),"bootstrap_high":np.quantile(bb,.975),"bootstrap_positive_fraction":np.mean(bb>0),"bootstrap_valid":len(bb)})
  draws[c+"__"+adj]=bb
R=pd.DataFrame(rows);R.to_csv(d/"RECEPTOR_INTERACTION_COHORT.tsv",sep="\t",index=False)
def meta(t):
 v=t.HC3_SE.to_numpy()**2;b=t.interaction_beta.to_numpy();w=1/v;mu=np.average(b,weights=w);Q=np.sum(w*(b-mu)**2);df=len(b)-1
 tau=max(0,(Q-df)/(w.sum()-(w*w).sum()/w.sum())) if len(b)>1 else 0
 ww=1/(v+tau);mu=np.average(b,weights=ww);se=1/np.sqrt(ww.sum())
 return {"cohorts":len(t),"n":int(t.n.sum()),"positive_cohorts":int((b>0).sum()),"random_interaction_beta":mu,"random_SE":se,"P":2*norm.sf(abs(mu/se)),"CI_low":mu-1.96*se,"CI_high":mu+1.96*se,"I2":max(0,(Q-df)/Q)*100 if Q>0 else 0,"tau2":tau}
M=[];L=[]
for adj,t in R.groupby("adjustment"):
 out=meta(t);v=t.HC3_SE.to_numpy()**2;w=1/(v+out["tau2"]);w=w/w.sum()
 boot=np.sum(np.stack([draws[c+"__"+adj] for c in t.cohort])*w[:,None],axis=0)
 M.append({"adjustment":adj,**out,"fixed_meta_weight_bootstrap_low":np.quantile(boot,.025),"fixed_meta_weight_bootstrap_high":np.quantile(boot,.975)})
 for c in t.cohort:L.append({"adjustment":adj,"omitted":c,**meta(t[t.cohort!=c])})
pd.DataFrame(M).to_csv(d/"RECEPTOR_INTERACTION_META.tsv",sep="\t",index=False);pd.DataFrame(L).to_csv(d/"RECEPTOR_INTERACTION_LEAVE_COHORT.tsv",sep="\t",index=False)
np.savez_compressed(d/"RECEPTOR_INTERACTION_BOOTSTRAP.npz",**draws)
manifest=[]
for p in files:
 hh=hashlib.sha256()
 with p.open("rb") as f:
  for block in iter(lambda:f.read(8*1024*1024),b""):hh.update(block)
 manifest.append({"path":str(p),"sha256":hh.hexdigest(),"bytes":p.stat().st_size})
(d/"INPUT_MANIFEST.json").write_text(json.dumps(manifest,indent=2))
pd.DataFrame(checks).to_csv(d/"INTERACTION_VALIDATION.tsv",sep="\t",index=False);assert all(x["pass"] for x in checks)
n=r/"part2_rebuild_v1/08_regulons/network_v2_20260915"
reg=pd.read_csv(n/"REGULON_REGISTRY.tsv",sep="\t")
gate={"STAT3_regulons":reg.loc[reg.TF=="STAT3","regulon"].tolist(),"STAT3_AUC_columns":[c for c in pd.read_csv(n/"AUC_ALL_CELLS.csv",nrows=0).columns if c.startswith("STAT3(")],"interpretation":"Absence of inferred regulon is measurement gap, not biological absence; no substitution with STAT3 RNA"}
(d/"STAT3_ACTIVITY_GATE.json").write_text(json.dumps(gate,indent=2))
print("COHORT\n"+R.to_csv(sep="\t",index=False));print("META\n"+pd.DataFrame(M).to_csv(sep="\t",index=False));print("VALIDATION",len(checks),gate)
