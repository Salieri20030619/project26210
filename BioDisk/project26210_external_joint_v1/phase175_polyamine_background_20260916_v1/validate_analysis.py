from pathlib import Path
import json,math,hashlib,datetime
import numpy as np,pandas as pd
from scipy.stats import rankdata,false_discovery_control
O=Path(__file__).resolve().parent;R=O.parent
checks=[];errs=[]
def ck(scope,name,v):
 checks.append(dict(scope=scope,check=name,passed=bool(v)))
 if not v:raise RuntimeError((scope,name))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
D=pd.read_csv(O/'MALIGNANT_PROGRAM_EFFECTS.tsv',sep='\t');S=pd.read_csv(O/'MALIGNANT_PROGRAM_SUMMARY.tsv',sep='\t')
with np.load(O/'MALIGNANT_CACHE.npz') as Z:
 for r in D.itertuples(index=False):
  ck(r.key,'frozen_technical_gate',r.cells>=50 and r.eligible and r.rank==r.columns and r.condition_number<30)
  C=Z[r.key+'_C'];x=Z[r.key+'_x'];y=Z[r.key+'_y'];Q,_=np.linalg.qr(C,mode='reduced')
  rr=np.corrcoef(x-Q@(Q.T@x),y-Q@(Q.T@y))[0,1]
  ck(r.key,'independent_QR_rho',abs(rr-r.rho)<1e-10)
for r in S.itertuples(index=False):
 a=D[(D.module==r.module)&(D.score==r.score)&(D.model==r.model)&D.eligible]
 v=a.rho.to_numpy();n=len(v)
 ck(r.module,'patient_n_effect_sign',a.patient.nunique()==r.n==n and (v>0).sum()==r.positive and abs(np.median(v)-r.median_rho)<1e-12)
 ranks=rankdata(abs(v));ck(r.module,'untied_nonzero_effects',len(set(abs(v)))==n and (v!=0).all())
 obs=ranks[v>0].sum();tail=int(min(obs,n*(n+1)/2-obs));dp=np.array([1],dtype=np.int64)
 for j in range(1,n+1):
  temp=np.zeros(len(dp)+j,dtype=np.int64);temp[:len(dp)]+=dp;temp[j:]+=dp;dp=temp
 pp=min(1.,2*dp[:tail+1].sum()/2**n);ck(r.module,'exact_signed_rank_P',abs(pp-r.P)<1e-12)
for key,a in S.groupby(['score','model']):ck(str(key),'BH2',np.allclose(false_discovery_control(a.P.to_numpy()),a.FDR2,atol=1e-12))
F=pd.read_csv(O/'PUBLIC_PROGRAM_SHORT_COHORT.tsv',sep='\t');M=pd.read_csv(O/'PUBLIC_PROGRAM_SHORT_META.tsv',sep='\t');E=pd.read_csv(O/'PUBLIC_PROGRAM_EDGES_META.tsv',sep='\t')
def qfit(a):
 v=rankdata(a,axis=0);v=(v-v.mean(0))/v.std(0,ddof=1);A,X,Y=v[:,:3].T;C=np.column_stack([np.ones(len(v)),v[:,3:]])
 def reg(x,y):
  q,r=np.linalg.qr(x,mode='reduced');return np.linalg.solve(r,q.T@y)
 aa=reg(np.column_stack([C,A]),X)[-1];yy=reg(np.column_stack([C,A,X]),Y);tt=reg(np.column_stack([C,A]),Y)[-1]
 return aa,yy[-1],aa*yy[-1],yy[-2],tt
with np.load(O/'PUBLIC_INPUT_CACHE.npz') as Z,np.load(O/'PUBLIC_BOOTSTRAPS.npz') as B:
 for r in F.itertuples(index=False):
  key=r.cohort+'__'+r.module+'__'+r.adjustment;a=Z[key+'_values'];bb=B[key]
  aa,b,ab,cp,tt=qfit(a)
  ck(key,'QR_path_identity',max(abs(aa-r.a),abs(b-r.b),abs(ab-r.ab),abs(cp-r.cprime),abs(tt-r.total),abs(tt-cp-ab))<1e-10)
  ck(key,'5000_bootstrap_CI',len(bb)==5000 and np.isfinite(bb).all() and np.allclose(np.quantile(bb,[.025,.975]),[r.ab_low,r.ab_high]))
  for i,ii in enumerate(Z[key+'_indices']):
   err=abs(qfit(a[ii])[2]-bb[i]);errs.append(err);ck(key,'bootstrap_QR',err<1e-9)
  for edge,cols in [('a',[0,1]+list(range(3,a.shape[1]))),('b',[1,2,0]+list(range(3,a.shape[1])))]:
   C=np.corrcoef(rankdata(a[:,cols],axis=0),rowvar=False);pr=np.linalg.inv(C);rho=-pr[0,1]/np.sqrt(pr[0,0]*pr[1,1]);ck(key,edge+'_precision',abs(rho-getattr(r,edge+'_rho'))<1e-10)
 for r in M.itertuples(index=False):
  g=F[(F.module==r.module)&(F.adjustment==r.adjustment)];v=g.ab_var.to_numpy();a=g.ab.to_numpy();w=1/v
  q=np.dot(w,(a-np.dot(w,a)/w.sum())**2);tau=max(0,(q-3)/(w.sum()-np.dot(w,w)/w.sum()));ww=1/(v+tau)
  mu=np.dot(ww,a)/ww.sum();se=1/np.sqrt(ww.sum());pv=math.erfc(abs(mu/se)/np.sqrt(2))
  ck(r.module+r.adjustment,'random_effect_meta',abs(mu-r.estimate)<1e-12 and abs(pv-r.P)<1e-12)
  combined=np.dot(ww/ww.sum(),np.stack([B[x.cohort+'__'+r.module+'__'+r.adjustment] for x in g.itertuples(index=False)]))
  ck(r.module+r.adjustment,'combined_bootstrap',np.allclose(combined,B['POOLED__'+r.module+'__'+r.adjustment]) and np.allclose(np.quantile(combined,[.025,.975]),[r.bootstrap_low,r.bootstrap_high]))
for key,a in M.groupby('adjustment'):ck(str(key),'short_BH2',np.allclose(false_discovery_control(a.P.to_numpy()),a.FDR2))
for key,a in E.groupby(['adjustment','edge']):ck(str(key),'edge_BH2',np.allclose(false_discovery_control(a.P.to_numpy()),a.FDR2))
for r in pd.read_csv(O/'INPUT_MANIFEST.tsv',sep='\t').itertuples(index=False):
 p=Path(r.path);ck(p.name,'source_unchanged',p.stat().st_size==r.bytes and p.stat().st_mtime_ns==r.mtime_ns and sha(p)==r.sha256)
ck('freeze','unchanged',sha(O/'PLAN_BEFORE_ANALYSIS_CN.md')==json.loads((O/'PLAN_FREEZE.json').read_text())['sha256'])
pd.DataFrame(checks).to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
out=dict(checks=len(checks),passed=sum(x['passed'] for x in checks),independent_bootstrap_draws=len(errs),max_bootstrap_QR_error=max(errs),completed_at=datetime.datetime.now().astimezone().isoformat())
(O/'VALIDATION_STATUS.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
g=pd.read_csv(O/'GENE_LEAVE_ONE_DIAGNOSTIC.tsv',sep='\t');print('DROP_GENE',g[g.model=='depth_cycle_hypoxia'].groupby('omitted_gene').agg(n=('patient','size'),median=('rho','median'),positive=('rho',lambda x:(x>0).sum())).to_string())
s=pd.read_csv(O/'MALIGNANT_STUDY_SENSITIVITY.tsv',sep='\t');print('STUDIES',s.query("model=='depth_cycle_hypoxia' and score=='gene_z_mean'").to_string(index=False))
j=pd.read_csv(O/'MATCHED_CORE_VS_SMS.tsv',sep='\t');print('MATCHED',j.groupby('model').agg(n=('patient','size'),core=('module_rho','median'),SMS=('SMS_rho','median'),paired_delta=('difference','median')).to_string())
print('LOO',pd.read_csv(O/'PUBLIC_LEAVE_ONE_COHORT.tsv',sep='\t').query("adjustment=='MTAP_myeloid'").to_string(index=False))
