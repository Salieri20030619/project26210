from pathlib import Path
import math,json,hashlib,datetime
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
O=Path(__file__).resolve().parent
d=pd.read_csv(O/'SECTION_EFFECTS.tsv',sep='\t');p=pd.read_csv(O/'PATIENT_EFFECTS.tsv',sep='\t');s=pd.read_csv(O/'PATIENT_SUMMARY.tsv',sep='\t')
checks=[];recomputed=[]
def ck(scope,name,value):
 checks.append(dict(scope=scope,check=name,passed=bool(value)))
 if not value:raise RuntimeError(str((scope,name)))
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
ck('ALL','expected_section_grid',len(d)==780 and not d.duplicated(['section','role','candidate','model']).any())
ck('ALL','primary_source_candidate_sections',len(d[(d.model=='extended_primary')&(d.role=='source')&d.eligible])==96)
ck('ALL','primary_receiver_candidate_sections',len(d[(d.model=='extended_primary')&(d.role=='receiver_link')&d.eligible])==80)
for gsm,a in d.groupby('section'):
 with np.load(O/'design_cache'/(gsm+'.npz')) as cache:
  for row in a.itertuples(index=False):
   X=cache[row.cache_key+'_X'];y=cache[row.cache_key+'_y'];ids=cache[row.cache_key+'_spots']
   ck(gsm+':'+row.cache_key,'cache_dimensions',len(y)==row.spots and len(ids)==row.spots and X.shape==(row.spots,row.columns) and len(np.unique(ids))==row.spots)
   if not row.eligible:continue
   C=X[:,np.r_[0,2:X.shape[1]]]
   Q,_=np.linalg.qr(C,mode='reduced')
   rx=X[:,1]-Q@(Q.T@X[:,1]);ry=y-Q@(Q.T@y)
   r=float(np.corrcoef(rx,ry)[0,1])
   ck(gsm+':'+row.cache_key,'QR_partial_correlation_agreement',abs(r-row.rho)<1e-10)
   ck(gsm+':'+row.cache_key,'technical_gate',row.spots>=50 and row.rank==row.columns and row.condition_number<30)
   recomputed.append(dict(section=gsm,cache_key=row.cache_key,rho_qr=r,abs_difference=abs(r-row.rho)))
for row in p.itertuples(index=False):
 a=d[(d.section.isin(row.section_ids.split(';')))&(d.role==row.role)&(d.candidate==row.candidate)&(d.model==row.model)]
 ck(row.patient_id+':'+row.model,'patient_source_consistency',len(a)==row.sections and a.eligible.all() and (a.patient_id==row.patient_id).all())
 v=a.rho.to_numpy()
 if len(v)==1:r=float(v[0])
 else:
  ck(row.patient_id,'only_c2_has_two_sections',row.patient_id.endswith(':c_2') and len(v)==2)
  odds=np.sqrt(np.prod((1+v)/(1-v)));r=float((odds-1)/(odds+1))
 ck(row.patient_id+':'+row.model,'patient_equal_weight_aggregation',abs(r-row.rho)<1e-10)
for row in s.itertuples(index=False):
 a=p[(p.role==row.role)&(p.candidate==row.candidate)&(p.model==row.model)&(p.stratum==row.stratum)&(p.policy==row.policy)]
 ck('SUMMARY','unique_patient_n',len(a)==row.n_patients and a.patient_id.nunique()==row.n_patients and a.sections.sum()==row.n_sections)
 if len(a):ck('SUMMARY','effect_and_sign_counts',np.isclose(a.rho.median(),row.median_rho,atol=1e-12) and (a.rho>0).sum()==row.positive and (a.rho<0).sum()==row.negative)
 if np.isfinite(row.P):
  k=min(row.positive,row.negative);n=row.test_n
  exact=min(1.0,2*sum(math.comb(n,j) for j in range(k+1))/(2**n))
  ck('SUMMARY','exact_sign_test_and_biological_n',n>=10 and abs(exact-row.P)<1e-12)
 else:ck('SUMMARY','planned_no_inference',row.n_patients<10 or row.stratum.startswith('cohort_') or row.stratum.startswith('exclude_cohort_') or row.model.startswith('component_'))
for key,a in s.groupby(['model','stratum','policy','family']):
 N=2 if key[-1]=='MTAP_background' else 3
 vals=np.r_[a.P.fillna(1).to_numpy(),np.ones(N-len(a))]
 fdr=false_discovery_control(vals,method='bh')[:len(a)]
 m=a.P.notna().to_numpy()
 ck('FDR:'+str(key),'independent_scipy_BH',np.allclose(fdr[m],a.FDR.to_numpy()[m],atol=1e-12))
for row in pd.read_csv(O/'INPUT_MANIFEST.tsv',sep='\t').itertuples(index=False):
 path=Path(row.path)
 ck(path.name,'input_size_mtime_unchanged',path.stat().st_size==row.bytes and path.stat().st_mtime_ns==row.mtime_ns)
 ck(path.name,'input_sha256_unchanged',sha(path)==row.sha256)
freeze=json.loads((O/'PLAN_FREEZE.json').read_text())
ck('PLAN','frozen_plan_hash_unchanged',sha(O/'FROZEN_ANALYSIS_PLAN.md')==freeze['sha256'])
pd.DataFrame(recomputed).to_csv(O/'QR_RECOMPUTED_EFFECTS.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'VALIDATION.tsv',sep='\t',index=False)
status=dict(checks=len(checks),passed=sum(x['passed'] for x in checks),max_qr_error=max(x['abs_difference'] for x in recomputed),validated_effects=len(recomputed),input_files=101,completed_at=datetime.datetime.now().astimezone().isoformat())
(O/'VALIDATION_STATUS.json').write_text(json.dumps(status,indent=2))
print(json.dumps(status))
print('SENSITIVITIES')
print(s[(s.stratum=='all')&(s.policy=='equal_fisher_mean')][['role','candidate','model','n_patients','positive','median_rho','P','FDR']].to_string(index=False))
print('PREP')
print(s[s.stratum.str.startswith('preparation_')][['role','candidate','stratum','n_patients','positive','median_rho','P','FDR']].to_string(index=False))
print('COHORT')
print(s[s.stratum.str.startswith('cohort_')][['role','candidate','stratum','n_patients','positive','median_rho']].to_string(index=False))
