from pathlib import Path
import warnings
import numpy as np,pandas as pd
from scipy.stats import wilcoxon,spearmanr,rankdata,t
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
D=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1'
def bh(p):
 p=np.asarray(p,float);out=np.full(p.shape,np.nan);ok=np.isfinite(p);v=p[ok];idx=np.argsort(v)
 if len(v):
  q=np.minimum.accumulate((v[idx]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];rev=np.empty(len(v));rev[idx]=np.minimum(q,1);out[ok]=rev
 return out
def summary(v):
 v=np.asarray(v,float);v=v[np.isfinite(v)];n=len(v)
 with warnings.catch_warnings():
  warnings.simplefilter('ignore')
  p=np.nan if n<3 else (float(wilcoxon(v,method='approx',correction=True).pvalue) if np.any(v) else 1.)
 return dict(n=n,median=float(np.median(v)) if n else np.nan,positive=int((v>0).sum()),negative=int((v<0).sum()),P=p)
def partial(x,y,covs=()):
 x=rankdata(x).astype(float);y=rankdata(y).astype(float)
 C=np.column_stack([np.ones(len(x))]+[rankdata(v) for v in covs])
 x-=C@np.linalg.lstsq(C,x,rcond=None)[0];y-=C@np.linalg.lstsq(C,y,rcond=None)[0]
 if np.std(x)<1e-10 or np.std(y)<1e-10:return np.nan
 return float(np.corrcoef(x,y)[0,1])
def paired_pvals(eff,test):
 p=np.full(eff.shape[1],np.nan)
 if eff.shape[0]>=3 and test.any():
  with warnings.catch_warnings():
   warnings.simplefilter('ignore');p[test]=wilcoxon(eff[:,test],axis=0,zero_method='wilcox',correction=True,method='approx').pvalue
 return p
def add_fdr(frame,keys,p='P',name='FDR'):
 frame[name]=np.nan
 for _,ix in frame.groupby(keys,dropna=False).groups.items():frame.loc[ix,name]=bh(frame.loc[ix,p])
 return frame
