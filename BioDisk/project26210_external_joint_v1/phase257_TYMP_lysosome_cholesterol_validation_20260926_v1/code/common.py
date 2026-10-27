from pathlib import Path
import sys,json,datetime,hashlib
import numpy as np,pandas as pd
from scipy import stats
P=Path(__file__).resolve().parents[1];ROOT=P.parent
sys.path.insert(0,str(P/'code'))
from phase256_input_helpers import gate,fingerprint
CHECKS=[]
def stamp():return datetime.datetime.now().astimezone().isoformat()
def out(name,d):d.to_csv(P/name,sep='\t',index=False)
def jwrite(name,x):(P/name).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda a:a.item() if hasattr(a,'item') else str(a)),encoding='utf-8')
def correlation(x,y,covs=None,min_n=50):
 x=np.asarray(x,float);y=np.asarray(y,float)
 covs={} if covs is None else covs
 keep=np.isfinite(x)&np.isfinite(y)
 for c in covs.values():keep &= np.isfinite(c)
 n=int(keep.sum());info=dict(n_observations=n,rho=np.nan,reason='',covariates=';'.join(covs),dropped_constant='',condition=np.nan,residual_df=np.nan)
 if n<min_n:info['reason']='INSUFFICIENT_OBSERVATIONS';return info,None
 x=x[keep];y=y[keep]
 if len(np.unique(x))<3 or len(np.unique(y))<3:info['reason']='INSUFFICIENT_VARIATION';return info,None
 rx=stats.rankdata(x);ry=stats.rankdata(y);cols=[];drop=[]
 for key,c in covs.items():
  v=stats.rankdata(np.asarray(c)[keep]);sd=v.std()
  if sd<1e-10:drop.append(key)
  else:cols.append((v-v.mean())/sd)
 C=np.column_stack([np.ones(n)]+cols);rank=np.linalg.matrix_rank(C);cond=float(np.linalg.cond(C));df=n-rank-2
 info.update(dropped_constant=';'.join(drop),condition=cond,residual_df=df)
 if rank!=C.shape[1] or cond>=30 or df<3:info['reason']='RANK_CONDITION_OR_DF';return info,None
 Q=np.linalg.qr(C,mode='reduced')[0]
 a=rx-Q@(Q.T@rx);b=ry-Q@(Q.T@ry)
 if min(a.std(),b.std())<1e-10:info['reason']='CONSTANT_RESIDUAL';return info,None
 rho=float(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))
 aa=rx-C@np.linalg.lstsq(C,rx,rcond=None)[0];bb=ry-C@np.linalg.lstsq(C,ry,rcond=None)[0]
 ref=float(stats.pearsonr(aa,bb).statistic);err=abs(ref-rho);assert err<1e-9
 CHECKS.append(err);info['rho']=rho
 return info,dict(a=a,b=b,Q=Q,keep=keep)
def cohort_summary(v,key):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v);pos=int((v>1e-12).sum());neg=int((v< -1e-12).sum());nz=pos+neg
 r=dict(patient_n=n,positive=pos,negative=neg,zero=n-nz,median_rho=float(np.median(v)) if n else np.nan,P=np.nan,CI_low=np.nan,CI_high=np.nan)
 if n>=10:
  r['P']=float(stats.binomtest(pos,nz,.5).pvalue) if nz else 1.
  seed=int(hashlib.sha256(key.encode()).hexdigest()[:8],16);rng=np.random.default_rng(seed)
  bs=np.median(v[rng.integers(n,size=(10000,n))],axis=1)
  r['CI_low'],r['CI_high']=np.quantile(bs,[.025,.975])
 return r
def correction(d,groups,pcol='P',qcol='FDR'):
 d[qcol]=np.nan
 for k,ix in d.groupby(groups,dropna=False).groups.items():
  ix=np.asarray(list(ix));pv=d.loc[ix,pcol];ok=pv.notna()
  d.loc[ix,'registered_family_size']=len(ix)
  if ok.any():
   q=stats.false_discovery_control(pv.fillna(1).to_numpy())
   d.loc[ix[ok],qcol]=q[ok]
 return d
