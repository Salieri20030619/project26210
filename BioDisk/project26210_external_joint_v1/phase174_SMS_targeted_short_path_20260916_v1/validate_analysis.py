from pathlib import Path
import json,hashlib,math,datetime
import numpy as np,pandas as pd
from scipy.stats import rankdata,t,norm,false_discovery_control
O=Path(__file__).resolve().parent;R=O.parent
checks=[];booterrs=[]
def ck(scope,name,val):
 checks.append(dict(scope=scope,check=name,passed=bool(val)))
 if not val:raise RuntimeError((scope,name))
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
D=pd.read_csv(O/'MALIGNANT_PATIENT_SMS_EFFECTS.tsv',sep='\t')
S=pd.read_csv(O/'MALIGNANT_SMS_SUMMARY.tsv',sep='\t')
with np.load(O/'MALIGNANT_DESIGN_CACHE.npz') as Z:
 for r in D.itertuples(index=False):
  C=Z[r.cache_key+'_C'];x=Z[r.cache_key+'_x'];y=Z[r.cache_key+'_y']
  ck(r.cache_key,'biological_patient_eligibility',r.cells>=50 and r.TYMP_detected>=10 and r.SMS_detected>=10)
  if r.eligible:
   q,_=np.linalg.qr(C,mode='reduced');xx=x-q@(q.T@x);yy=y-q@(q.T@y);rr=np.corrcoef(xx,yy)[0,1]
   ck(r.cache_key,'QR_partial_correlation',abs(rr-r.rho)<1e-10)
for r in S.itertuples(index=False):
 a=D[(D.model==r.model)&D.eligible];v=a.rho.to_numpy()
 ck(r.model,'patient_summary',len(a)==r.n_patients and a.patient.nunique()==r.n_patients and (v>0).sum()==r.positive and abs(np.median(v)-r.median_rho)<1e-12)
 # exact signed-rank distribution for n=16, no zero/tied absolute effects here
 ranks=rankdata(abs(v));ck(r.model,'unique_nonzero_signed_rank_inputs',(v!=0).all() and len(set(abs(v)))==len(v))
 obs=float(ranks[v>0].sum());n=len(v);weights=np.arange(1,n+1);dp=np.array([1],dtype=np.int64)
 for j in weights:
  new=np.zeros(len(dp)+j,dtype=np.int64);new[:len(dp)]+=dp;new[j:]+=dp;dp=new
 tail=min(obs,float(weights.sum()-obs));pp=min(1.,2*dp[:int(tail)+1].sum()/2**n)
 ck(r.model,'independent_exact_Wilcoxon',abs(pp-r.P)<1e-12)
new=S[~S.model.isin(['raw','depth','depth_cycle'])]
ck('malignant','BH4',np.allclose(false_discovery_control(new.P.to_numpy()),new.supplement_FDR_4))
F=pd.read_csv(O/'PUBLIC_SHORT_PATH_COHORT.tsv',sep='\t');M=pd.read_csv(O/'PUBLIC_SHORT_PATH_META.tsv',sep='\t');E=pd.read_csv(O/'PUBLIC_CONDITIONAL_EDGES_META.tsv',sep='\t')
def qrfit(arr):
 z=rankdata(arr,axis=0);z=(z-z.mean(0))/z.std(0,ddof=1)
 a,x,y=z[:,:3].T;C=np.column_stack([np.ones(len(arr)),z[:,3:]])
 def reg(d,v):
  q,r=np.linalg.qr(d,mode='reduced');return np.linalg.solve(r,q.T@v)
 aa=reg(np.column_stack([C,a]),x)[-1];bb=reg(np.column_stack([C,a,x]),y);tot=reg(np.column_stack([C,a]),y)[-1]
 return aa,float(bb[-1]),float(aa*bb[-1]),float(bb[-2]),float(tot)
with np.load(O/'PUBLIC_INPUT_CACHE.npz') as Z,np.load(O/'PUBLIC_AB_BOOTSTRAPS.npz') as B:
 for r in F.itertuples(index=False):
  key=r.cohort+'__'+r.adjustment;arr=Z[key+'_values'];bb=B[key];ids=Z[key+'_ids']
  aa,b,ab,cp,total=qrfit(arr)
  ck(key,'unique_patients',len(ids)==r.n and len(set(ids))==r.n)
  ck(key,'QR_coefficients_and_identity',max(abs(aa-r.a),abs(b-r.b),abs(ab-r.ab),abs(cp-r.cprime),abs(total-r.total),abs(total-cp-ab))<1e-10)
  ck(key,'5000_bootstrap_and_intervals',len(bb)==5000 and np.isfinite(bb).all() and np.allclose(np.quantile(bb,[.025,.975]),[r.ab_low,r.ab_high],atol=1e-12))
  for j,ii in enumerate(Z[key+'_indices']):
   b0=qrfit(arr[ii])[2];err=abs(b0-bb[j]);booterrs.append(err);ck(key+':'+str(j),'independent_bootstrap_rerank_QR',err<1e-9)
  for edge,cols,k in [('a',[0,1]+list(range(3,arr.shape[1])),arr.shape[1]-3),('b',[1,2,0]+list(range(3,arr.shape[1])),arr.shape[1]-2)]:
   c=np.corrcoef(rankdata(arr[:,cols],axis=0),rowvar=False);pr=np.linalg.inv(c);rho=-pr[0,1]/np.sqrt(pr[0,0]*pr[1,1])
   ck(key+edge,'precision_partial_correlation',abs(rho-getattr(r,edge+'_rho'))<1e-10)
   pp=2*t.sf(abs(rho)*np.sqrt((len(arr)-k-2)/(1-rho*rho)),len(arr)-k-2)
   ck(key+edge,'edge_P_value',abs(pp-getattr(r,edge+'_P'))<1e-10)
 for r in M.itertuples(index=False):
  g=F[F.adjustment==r.adjustment]
  w=1/g.ab_bootstrap_variance.to_numpy();eff=g.ab.to_numpy();fixed=np.dot(w,eff)/w.sum();q=np.dot(w,(eff-fixed)**2)
  tau=max(0,(q-3)/(w.sum()-np.dot(w,w)/w.sum()));rw=1/(1/w+tau);mu=np.dot(rw,eff)/rw.sum();se=1/np.sqrt(rw.sum());pp=math.erfc(abs(mu/se)/np.sqrt(2))
  ck(r.adjustment,'independent_meta_ab',abs(mu-r.estimate)<1e-12 and abs(se-r.SE)<1e-12 and abs(pp-r.P)<1e-12)
  arr=np.stack([B[x.cohort+'__'+r.adjustment] for x in g.itertuples(index=False)])
  pooled=np.dot(rw/rw.sum(),arr)
  ck(r.adjustment,'pooled_bootstrap',np.allclose(pooled,B['POOLED__'+r.adjustment],atol=1e-12) and np.allclose(np.quantile(pooled,[.025,.975]),[r.bootstrap_low,r.bootstrap_high],atol=1e-12))
ix=M.adjustment!='MTAP_myeloid';ck('meta','sensitivity_BH3',np.allclose(false_discovery_control(M.loc[ix,'P'].to_numpy()),M.loc[ix,'sensitivity_FDR3']))
hist=pd.read_csv(R/'phase162_MTA_polyamine_part1_reassessment_20260915_v1/PUBLIC_META_AXIS_ASSOCIATIONS.tsv',sep='\t')
EH=E.copy()
for i,r in EH.iterrows():
 if r.edge=='a' and r.adjustment=='raw':hh=hist[(hist.gene=='SMS')&(hist.adjustment=='raw')&(hist.edge=='a')&hist.outcome.str.contains('MTA')]
 elif r.edge=='b' and r.adjustment in ['raw','MTAP','MTAP_myeloid']:
  adj='myeloid' if r.adjustment=='MTAP_myeloid' else 'raw'
  edge='b' if r.adjustment=='raw' else 'b_adjust_MTAP'
  hh=hist[(hist.gene=='SMS')&(hist.adjustment==adj)&(hist.edge==edge)&hist.outcome.str.contains('MTA')]
 else:continue
 ck('historical','same_edge_rho',len(hh)==1 and abs(float(hh.iloc[0].random_rho)-r.rho)<1e-10)
 EH.loc[i,'historical_FDR6']=float(hh.iloc[0].random_FDR)
EH.to_csv(O/'PUBLIC_EDGE_EVIDENCE_HISTORY.tsv',sep='\t',index=False)
for row in pd.read_csv(O/'INPUT_MANIFEST.tsv',sep='\t').itertuples(index=False):
 p=Path(row.path);ck(p.name,'input_identity_unchanged',p.stat().st_size==row.size and p.stat().st_mtime_ns==row.mtime_ns and sha(p)==row.sha256)
freeze=json.loads((O/'PLAN_FREEZE.json').read_text());ck('plan','freeze_hash',sha(O/'PLAN_BEFORE_ANALYSIS_CN.md')==freeze['sha256'])
pd.DataFrame(checks).to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
out=dict(checks=len(checks),passed=sum(x['passed'] for x in checks),max_bootstrap_QR_error=max(booterrs),independent_bootstrap_draws=len(booterrs),validated_at=datetime.datetime.now().astimezone().isoformat())
(O/'VALIDATION_STATUS.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out));print(S.to_string(index=False))
print('LEAVE COHORT');print(pd.read_csv(O/'PUBLIC_LEAVE_ONE_COHORT.tsv',sep='\t').query("adjustment=='MTAP_myeloid'").to_string(index=False))
print('MALIGNANT study primary');print(pd.read_csv(O/'MALIGNANT_STUDY_SENSITIVITY.tsv',sep='\t').query("model=='depth_cycle_hypoxia'").to_string(index=False))
print('PRIMARY_PATIENT_JACKKNIFE');print(pd.read_csv(O/'PUBLIC_LEAVE_ONE_PATIENT.tsv',sep='\t').groupby('cohort').ab.agg(['min','max']).to_string())
