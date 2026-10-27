from pathlib import Path
import pandas as pd,numpy as np,json,subprocess,hashlib
from scipy import stats
P=Path(__file__).resolve().parents[1]
with (P/'logs/memory_gate.log').open('a') as f:subprocess.run(['/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],stdout=f,check=True)
freeze=dict(protocol_sha256=hashlib.sha256((P/'INTERACTION_ADDENDUM_BEFORE_TESTING.md').read_bytes()).hexdigest(),input_sha256=hashlib.sha256((P/'tables/FROZEN_PATIENT_INPUTS.tsv').read_bytes()).hexdigest())
(P/'qc/INTERACTION_FREEZE.json').write_text(json.dumps(freeze,indent=2))
d=pd.read_csv(P/'tables/FROZEN_PATIENT_INPUTS.tsv',sep='\t')
def design(z,source):
 def zrank(v):
  r=stats.rankdata(v);return (r-r.mean())/r.std(ddof=0)
 a=zrank(z[source]);b=zrank(z.macro_SPP1);y=zrank(z.collagen_maturation)
 D=np.column_stack([np.ones(len(z))]+[(z.study==s).astype(float) for s in sorted(z.study.unique())[1:]]+[a,b])
 w=a*b;X=np.column_stack([D,w]);return D,w,X,y
out=[];leave=[];checks=[]
for source in ['malignant_TYMP','macro_TYMP']:
 D,w,X,y=design(d,source);n=len(y);rank=np.linalg.matrix_rank(X);df=n-rank
 beta=np.linalg.lstsq(X,y,rcond=None)[0];u=w-D@np.linalg.lstsq(D,w,rcond=None)[0];v=y-D@np.linalg.lstsq(D,y,rcond=None)[0]
 independent=u@v/(u@u);assert np.isclose(beta[-1],independent,atol=1e-12)
 denom=np.sqrt((v-(u@v/(u@u))*u)@(v-(u@v/(u@u))*u)/df/(u@u));t=beta[-1]/denom
 p=np.nan;B=99999
 if n>=10 and df>=5:
  rng=np.random.default_rng(29601);idx=np.tile(np.arange(n),(B,1))
  for st,g in d.groupby('study'):
   ix=g.index.to_numpy()
   if len(ix)>1:idx[:,ix]=ix[np.argsort(rng.random((B,len(ix))),axis=1)]
  V=v[idx];Q=np.linalg.qr(D,mode='reduced')[0];V-=V@Q@Q.T
  bp=V@u/(u@u);err=V-bp[:,None]*u
  se=np.sqrt((err*err).sum(1)/df/(u@u));tt=bp/se
  p=(np.count_nonzero(np.abs(tt)>=abs(t)-1e-12)+1)/(B+1)
 out.append(dict(source=source,n=n,design_rank=rank,residual_df=df,TYMP_main_beta=beta[-3],SPP1_main_beta=beta[-2],interaction_beta=beta[-1],interaction_t=t,P=p,SPP1_slope_TYMP_minus1=beta[-2]-beta[-1],SPP1_slope_TYMP_mean=beta[-2],SPP1_slope_TYMP_plus1=beta[-2]+beta[-1],method='exploratory_FreedmanLane_study_blocked_approximate',permutations=B))
 checks.append(dict(source=source,independent_beta=independent,full_model_beta=beta[-1],error=abs(beta[-1]-independent)))
 for typ in ['patient','study']:
  for val in d[typ].unique():
   z=d[d[typ]!=val];dd,ww,xx,yy=design(z,source)
   ok=np.linalg.matrix_rank(xx)==xx.shape[1] and len(z)>xx.shape[1]
   bb=np.linalg.lstsq(xx,yy,rcond=None)[0][-1] if ok else np.nan
   leave.append(dict(source=source,omission_type=typ,omitted=val,n=len(z),interaction_beta=bb,status='DESCRIPTIVE' if ok else 'INSUFFICIENT_DF'))
o=pd.DataFrame(out);o['FDR']=stats.false_discovery_control(o.P.to_numpy());o['family']='TYMP_SPP1_interaction_two_sources';o['family_size']=2
o.to_csv(P/'tables/INTERACTION_EXPLORATORY.tsv',sep='\t',index=False)
pd.DataFrame(leave).to_csv(P/'tables/INTERACTION_LEAVE_OUT.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(P/'qc/INTERACTION_INDEPENDENT_CHECKS.tsv',sep='\t',index=False)
print(o.to_string(index=False))
