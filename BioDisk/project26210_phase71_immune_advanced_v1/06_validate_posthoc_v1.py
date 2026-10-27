#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase71_immune_advanced_v1');RES=ROOT/'posthoc_state_results_v1';C=[]
def add(i,o,d):C.append({'check_id':i,'status':'PASS' if o else 'FAIL','detail':str(d)})
def read(n):return pd.read_csv(RES/n,sep='\t')
def bh(p):p=np.array(p,float);n=len(p);o=np.argsort(p);r=p[o];q=np.minimum.accumulate((r*n/np.arange(1,n+1))[::-1])[::-1];z=np.empty(n);z[o]=np.minimum(q,1);return z
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
s=read('state_adjusted_summary_all_v1.tsv');st=read('state_stratified_summary_all_v1.tsv');l=read('response_loo_summary_all_v1.tsv');sc=read('state_adjusted_cross_cohort_v1.tsv');lr=read('response_loo_robustness_v1.tsv');status=read('state_audit_status_all_v1.tsv')
add('status',len(status)==3 and status.status.eq('PASS_STATE_CONFOUNDING_AUDIT_V1').all(),len(status));add('dimensions',len(s)==18 and len(l)==36,len(s));qok=True
for d,g in s.groupby('dataset'):qok&=np.allclose(g.q_value,bh(g.p_value));
for d,g in st.groupby('dataset'):qok&=np.allclose(g.q_value,bh(g.p_value));
for d,g in l.groupby('dataset'):qok&=np.allclose(g.q_value,bh(g.p_value));
add('bh',qok,qok);dec=True
for r in s.itertuples(index=False):dec&=(r.decision==('ASSOCIATED' if r.expected=='two_sided' and r.q_value<.05 else ('SUPPORTED' if r.expected!='two_sided' and r.q_value<.05 and r.median_rho>0 and r.positive_fraction>=2/3 else 'NOT_SUPPORTED')))
add('state_decisions',dec,dec);add('loo_robust_rows',len(lr)==3,len(lr));manifest=read('FILE_MANIFEST_SHA256_v1.tsv');hok=all((RES/r.file).is_file() and sha(RES/r.file)==r.sha256 for r in manifest.itertuples(index=False));add('hashes',hok,len(manifest))
out=ROOT/'posthoc_independent_validation_v1';out.mkdir(exist_ok=True);pd.DataFrame(C).to_csv(out/'validation_checks_v1.tsv',sep='\t',index=False,lineterminator='\n');x={'checks':len(C),'passed':sum(z['status']=='PASS' for z in C),'failed':sum(z['status']=='FAIL' for z in C)};(out/'validation_summary_v1.json').write_text(json.dumps(x,indent=2)+'\n');
if x['failed']==0:(out/'POSTHOC_INDEPENDENT_VALIDATION_PASS_v1').write_text('PASS\n')
print(json.dumps(x));raise SystemExit(0 if x['failed']==0 else 1)
