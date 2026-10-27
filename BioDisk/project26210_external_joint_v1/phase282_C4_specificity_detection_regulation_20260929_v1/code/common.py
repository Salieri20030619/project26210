from pathlib import Path
import json,hashlib,datetime,importlib.util,sys,subprocess
import numpy as np,pandas as pd,scipy
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent
PRE=R/'phase280_TYMP_source_niche_function_20260928_v1'
M=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1'
spec=importlib.util.spec_from_file_location('previous',PRE/'code/source.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);b.P=P
infer=b.infer;sha=b.sha;gate=b.gate;bh=b.bh
CONTEXT=['Interferon Alpha Response','Interferon Gamma Response','Hypoxia']
GENES=['SPP1','CTSB','CTSD','LGALS3']
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def read(f,**kw):return pd.read_csv(f,sep='\t',**kw)
def tab(d,name):pd.DataFrame(d).to_csv(P/'tables'/name,sep='\t',index=False)
def data():
 a=read(PRE/'inputs/CELL_CACHE.tsv.gz',dtype={'candidate':str})
 c=read(PRE/'inputs/CONTEXT_CELL_CACHE.tsv.gz')
 z=a.merge(c,on=['cell_id','patient'],validate='one_to_one')
 assert len(z)==20426 and z.cell_id.is_unique
 return z
def covariates(z,model,other=False,extra=None):
 cols=[z.log_depth.to_numpy(),z.n_genes.to_numpy()]
 if model=='context':cols += [z[c].to_numpy() for c in CONTEXT]
 if other:cols += [v for v in pd.get_dummies(z.candidate,dtype=float,drop_first=True).to_numpy().T]
 if extra is not None:cols.append(np.asarray(extra))
 return np.column_stack(cols)
def correlation(x,y,cov=None,min_n=50,binary=False):
 x=np.asarray(x,float);y=np.asarray(y,float);cov=np.empty((len(x),0)) if cov is None else np.asarray(cov,float)
 keep=np.isfinite(x)&np.isfinite(y)&np.isfinite(cov).all(1);x=x[keep];y=y[keep];cov=cov[keep];n=len(x)
 if n<min_n:return np.nan,n,'INSUFFICIENT_CELLS'
 if min(len(np.unique(x)),len(np.unique(y)))<(2 if binary else 3):return np.nan,n,'INSUFFICIENT_VARIATION'
 v=np.column_stack([x,y]) if binary else np.column_stack([stats.rankdata(x),stats.rankdata(y)])
 c=np.column_stack([stats.rankdata(v) for v in cov.T]) if cov.shape[1] else cov
 sd=c.std(0)
 if (sd==0).any():return np.nan,n,'CONSTANT_COVARIATE'
 D=np.column_stack([np.ones(n),(c-c.mean(0))/sd]) if c.shape[1] else np.ones((n,1))
 if np.linalg.matrix_rank(D)<D.shape[1] or np.linalg.cond(D)>=30:return np.nan,n,'RANK_OR_CONDITION'
 residual=v-D@np.linalg.lstsq(D,v,rcond=None)[0]
 if (residual.std(0)<1e-10).any():return np.nan,n,'RESIDUAL_CONSTANT'
 return float(np.corrcoef(residual.T)[0,1]),n,''
def summarize(e,keys,familykeys,name):
 rows=[]
 for key,z in e.groupby(keys,dropna=False):
  if not isinstance(key,tuple):key=(key,)
  rows.append(dict(zip(keys,key),**infer(z.effect,26212820+len(rows))))
 s=bh(pd.DataFrame(rows),familykeys);tab(s,name)
 return s
def study_descriptions(e,keys,prefix):
 rows=[];leave=[]
 for key,g in e.groupby(keys,dropna=False):
  if not isinstance(key,tuple):key=(key,)
  for st in sorted(g.study.unique()):
   for out,sel in [(rows,g.study.eq(st)),(leave,g.study.ne(st))]:
    v=g.loc[sel,'effect'].dropna();out.append(dict(zip(keys,key),study=st,n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=float(v.median()) if len(v) else np.nan))
 tab(rows,prefix+'_STUDY.tsv');tab(leave,prefix+'_LEAVE_STUDY.tsv')
