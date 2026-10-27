from pathlib import Path
import json, hashlib, datetime, sys
import numpy as np, pandas as pd
from scipy import sparse, stats
P=Path(__file__).resolve().parents[1]; R=P.parent
F=R/'phase229_spatial_F1_partition_eligibility_20260921_v1'
A=R/'phase223_spatial_cohort_reinclusion_20260920_v1'
S=R/'phase221_spatial_S1_identity_eligibility_freeze_20260920_v1'
def save(name,obj):
 (P/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def tab(name,d):d.to_csv(P/name,sep='\t',index=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stamp():return datetime.datetime.now().astimezone().isoformat()
def meta(p):
 p=Path(p);s=p.stat();return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino)
def fit(x,y,cov=None,mask=None,min_n=50):
 x=np.asarray(x);y=np.asarray(y)
 if cov is None:cov=np.empty((len(x),0))
 mask=np.ones(len(x),bool) if mask is None else mask.copy()
 mask &= np.isfinite(x)&np.isfinite(y)&np.isfinite(cov).all(axis=1)
 a=np.column_stack([x,y,cov])[mask];n=len(a)
 if n<min_n:return np.nan,n,'INSUFFICIENT_OBSERVATIONS'
 if min(len(np.unique(a[:,0])),len(np.unique(a[:,1])))<3:return np.nan,n,'INSUFFICIENT_VARIATION'
 rk=np.column_stack([stats.rankdata(a[:,j]) for j in range(a.shape[1])])
 c=rk[:,2:];sd=c.std(axis=0)
 if (sd==0).any():return np.nan,n,'CONSTANT_COVARIATE'
 D=np.column_stack([np.ones(n),(c-c.mean(axis=0))/sd]) if c.shape[1] else np.ones((n,1))
 if np.linalg.matrix_rank(D)<D.shape[1] or np.linalg.cond(D)>=30:return np.nan,n,'RANK_OR_CONDITION'
 e=rk[:,:2]-D@np.linalg.lstsq(D,rk[:,:2],rcond=None)[0]
 if (e.std(axis=0)<1e-10).any():return np.nan,n,'RESIDUAL_CONSTANT'
 v=float(np.corrcoef(e.T)[0,1])
 if c.shape[1]:
  q=np.linalg.qr(D,mode='reduced')[0];ee=rk[:,:2]-q@(q.T@rk[:,:2])
  assert abs(float(np.corrcoef(ee.T)[0,1])-v)<1e-10
 return v,n,''
def infer(v,seed=262234):
 v=np.asarray(v);v=v[np.isfinite(v)];n=len(v);nz=int((v!=0).sum())
 o=dict(n=n,positive=int((v>0).sum()),negative=int((v<0).sum()),zero=n-nz,median=float(np.median(v)) if n else np.nan,P=np.nan,P_wilcoxon=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=10:
  o['P']=float(stats.binomtest(o['positive'],nz,.5).pvalue) if nz else 1.
  o['P_wilcoxon']=float(stats.wilcoxon(v,method='auto').pvalue) if nz else 1.
  rng=np.random.default_rng(seed);boot=np.median(v[rng.integers(n,size=(10000,n))],axis=1)
  o['CI_low'],o['CI_high']=np.quantile(boot,[.025,.975])
 return o
KEYS=['family','endpoint','exposure','target','layer','scale','qualification','model','effect_type']
def summarize(sec,prefix):
 pats=[]
 for key,g in sec.groupby(['dataset','patient']+KEYS,dropna=False,sort=True):
  v=g.effect.dropna().to_numpy()
  effect=(float(np.tanh(np.arctanh(np.clip(v,-.999999999,.999999999)).mean())) if key[-1]=='rho' else float(np.mean(v))) if len(v) else np.nan
  pats.append(dict(zip(['dataset','patient']+KEYS,key),effect=effect,n_sections=len(v),sections=';'.join(g.loc[g.effect.notna(),'sample'])))
 pat=pd.DataFrame(pats);tab('tables/'+prefix+'_PATIENT_EFFECTS.tsv',pat)
 summary=[]
 for i,(key,g) in enumerate(pat.groupby(['dataset']+KEYS,dropna=False,sort=True)):
  summary.append(dict(zip(['dataset']+KEYS,key),**infer(g.effect,262234+i)))
 su=pd.DataFrame(summary);su['FDR']=np.nan;su['FDR_wilcoxon']=np.nan
 # Family pooling includes all registered parallel variables and every scale/model within layer.
 for key,g in su.groupby(['dataset','family','layer'],dropna=False):
  su.loc[g.index,'family_size_registered']=len(g)
  su.loc[g.index,'family_size_estimable']=int(g.P.notna().sum())
  for pcol,qcol in [('P','FDR'),('P_wilcoxon','FDR_wilcoxon')]:
   q=stats.false_discovery_control(g[pcol].fillna(1).to_numpy())
   su.loc[g.index,qcol]=np.where(g[pcol].notna(),q,np.nan)
 su['support']=np.where(su.P.isna(),'DESCRIPTIVE_OR_NOT_ESTIMABLE',np.where(su.P<.05,np.where(su.FDR<.05,'FDR_SUPPORTED','NOMINAL_ONLY'),'NOT_SUPPORTED'))
 tab('tables/'+prefix+'_COHORT_SUMMARY.tsv',su)
 return pat,su
