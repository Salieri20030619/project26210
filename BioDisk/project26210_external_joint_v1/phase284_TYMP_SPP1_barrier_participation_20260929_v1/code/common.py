from pathlib import Path
import json,hashlib,datetime,importlib.util,sys,subprocess
import numpy as np,pandas as pd,scipy
from scipy import sparse,stats
P=Path(__file__).resolve().parents[1];R=P.parent
PRE=R/'phase280_TYMP_source_niche_function_20260928_v1'
M=R/'phase200_union_C4_C10_C15_reassessment_20260918_v1'
spec=importlib.util.spec_from_file_location('base280',PRE/'code/source.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b);b.P=P
fit=b.fit;infer=b.infer;cross=b.cross;bh=b.bh;sha=b.sha;gate=b.gate
def read(f,**kw):return pd.read_csv(f,sep='\t',float_precision='round_trip',**kw)
def tab(d,name):pd.DataFrame(d).to_csv(P/'tables'/name,sep='\t',index=False)
def dump(f,o):Path(f).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
def summary(e,keys,families,filename):
 rows=[]
 for key,d in e.groupby(keys,dropna=False):
  if not isinstance(key,tuple):key=(key,)
  rows.append(dict(zip(keys,key),**infer(d.effect,26212840+len(rows))))
 s=bh(pd.DataFrame(rows),families);tab(s,filename);return s
def describe_study(e,keys,prefix):
 rows=[];leave=[]
 for key,d in e.groupby(keys,dropna=False):
  if not isinstance(key,tuple):key=(key,)
  for st in sorted(d.study.unique()):
   for out,flag in [(rows,d.study.eq(st)),(leave,d.study.ne(st))]:
    v=d.loc[flag,'effect'].dropna()
    out.append(dict(zip(keys,key),study=st,n=len(v),positive=int((v>0).sum()),negative=int((v<0).sum()),median=v.median() if len(v) else np.nan))
 tab(rows,prefix+'_STUDY.tsv');tab(leave,prefix+'_LEAVE_STUDY.tsv')
