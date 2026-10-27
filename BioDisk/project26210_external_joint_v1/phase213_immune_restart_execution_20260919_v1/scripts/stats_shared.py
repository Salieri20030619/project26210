import numpy as np,pandas as pd,math,itertools,collections
from scipy.stats import rankdata,spearmanr

def groups(study):return [np.flatnonzero(np.asarray(study)==s) for s in sorted(set(study))]
def residual_ranks(x,study,cov=None):
 r=rankdata(x).astype(float)
 Z=np.column_stack([np.ones(len(x))]+[(np.asarray(study)==s).astype(float) for s in sorted(set(study))[1:]]+([] if cov is None else [rankdata(cov).astype(float)]))
 return r-Z@np.linalg.lstsq(Z,r,rcond=None)[0],np.linalg.matrix_rank(Z)
def effect(x,y,study,cov=None):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if len(x)<3 or not np.isfinite(x).all() or not np.isfinite(y).all():return np.nan
 a,k=residual_ranks(x,study,cov);b,_=residual_ranks(y,study,cov);den=np.linalg.norm(a)*np.linalg.norm(b)
 if den<1e-10 or len(x)-k<2:return np.nan
 return float(np.dot(a,b)/den)
def unique_perms(values):
 c=collections.Counter(values);keys=sorted(c);n=len(values);tmp=[]
 def rec():
  if len(tmp)==n:yield tuple(tmp);return
  for k in keys:
   if c[k]:
    c[k]-=1;tmp.append(k);yield from rec();tmp.pop();c[k]+=1
 yield from rec()
def main_test(x,y,study,seed=262212,B=99999):
 n=len(x);a,k=residual_ranks(x,study);b,_=residual_ranks(y,study);obs=effect(x,y,study);gs=groups(study);nperm=1
 for g in gs:
  v,c=np.unique(b[g],return_counts=True);nperm*=math.factorial(len(g))//math.prod(math.factorial(int(z)) for z in c)
 if not np.isfinite(obs):return {'rho':np.nan,'p':np.nan,'method':'not_estimable','permutations':0,'possible_permutations':str(nperm),'study_design_rank':int(k),'residual_df':n-int(k)-1},np.array([])
 den=np.linalg.norm(a)*np.linalg.norm(b)
 if nperm<=B:
  options=[list(unique_perms(b[g].tolist())) for g in gs];dist=np.empty(nperm)
  for j,pp in enumerate(itertools.product(*options)):
   v=b.copy()
   for g,p in zip(gs,pp):v[g]=p
   dist[j]=np.dot(a,v)/den
  p=float(np.mean(np.abs(dist)>=abs(obs)-1e-12));method='exact_study_blocked'
 else:
  rng=np.random.default_rng(seed);dist=np.empty(B)
  for start in range(0,B,10000):
   m=min(10000,B-start);v=np.tile(b,(m,1))
   for g in gs:
    if len(g)>1:v[:,g]=b[g][np.argsort(rng.random((m,len(g))),axis=1)]
   dist[start:start+m]=v@a/den
  p=(1+int((np.abs(dist)>=abs(obs)-1e-12).sum()))/(B+1);method='MC_study_blocked_plus1'
 return {'rho':obs,'p':p,'method':method,'permutations':len(dist),'possible_permutations':str(nperm),'study_design_rank':int(k),'residual_df':n-int(k)-1},dist

def bootstrap_effect(x,y,study,B=5000,seed=263212):
 x=np.asarray(x);y=np.asarray(y);gs=groups(study);rng=np.random.default_rng(seed);ix=np.tile(np.arange(len(x)),(B,1))
 for g in gs:ix[:,g]=rng.choice(g,size=(B,len(g)),replace=True)
 a=rankdata(x[ix],axis=1).astype(float);b=rankdata(y[ix],axis=1).astype(float)
 for g in gs:
  a[:,g]-=a[:,g].mean(axis=1,keepdims=True);b[:,g]-=b[:,g].mean(axis=1,keepdims=True)
 den=np.sqrt((a*a).sum(axis=1)*(b*b).sum(axis=1));v=np.full(B,np.nan);good=den>1e-10;v[good]=(a[good]*b[good]).sum(axis=1)/den[good]
 ci=np.quantile(v[good],[.025,.975]) if good.mean()>=.8 else [np.nan,np.nan]
 return float(ci[0]),float(ci[1]),int(good.sum()),v

def bh(values):
 p=np.asarray(values,float);order=np.argsort(p);q=p[order]*len(p)/np.arange(1,len(p)+1);q=np.minimum.accumulate(q[::-1])[::-1];out=np.empty(len(p));out[order]=np.clip(q,0,1);return out
