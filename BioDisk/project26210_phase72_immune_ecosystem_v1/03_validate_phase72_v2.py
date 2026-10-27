#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
ROOT=Path('/mnt/BioDisk/wangmingyuan/project26210_phase72_immune_ecosystem_v1');RES=ROOT/'results_v1';C=[]
def add(i,o,d):C.append({'check_id':i,'status':'PASS' if o else 'FAIL','detail':str(d)})
def read(n):return pd.read_csv(RES/n,sep='\t')
def bh(p):return stats.false_discovery_control(np.asarray(p,float),method='bh')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
s=read('phase72_primary_relations_all_v1.tsv');x=read('phase72_cross_cohort_decisions_v1.tsv');status=read('phase72_dataset_status_all_v1.tsv');cov=read('phase72_program_coverage_all_v1.tsv');tri=read('PHASE70_PHASE72_TRIANGULATION_v1.tsv')
add('status',len(status)==3 and status.status.eq('PASS_PHASE72_DATASET_V1').all(),len(status));add('dimensions',len(s)==24 and len(x)==9,(len(s),len(x)));add('coverage',len(cov)==18 and cov.n_present.eq(cov.n_frozen).all(),len(cov));qok=True;dec=True
for d,g in s.groupby('dataset'):qok&=np.allclose(g.q_value,bh(g.p_value))
for r in s.itertuples(index=False):dec&=(r.decision==('COHORT_SUPPORTED' if r.q_value<.05 and ((r.expected=='positive' and r.rho>0)or(r.expected=='negative' and r.rho<0)) else ('DIRECTION_CONFLICT' if r.q_value<.05 else 'NOT_SUPPORTED')))
add('bh',qok,qok);add('decisions',dec,dec);crossok=True
for r in x.itertuples(index=False):
 ds=[d for d in ('GSE264586','GSE294109') if hasattr(r,f'{d}_decision') and pd.notna(getattr(r,f'{d}_decision'))];n=sum(getattr(r,f'{d}_decision')=='COHORT_SUPPORTED' for d in ds);expected=('TWO_EXTERNAL_SUPPORTED' if n==2 else ('ONE_EXTERNAL_SUPPORTED' if n==1 else 'NOT_SUPPORTED')) if len(ds)==2 else ('GSE294109_SUPPORTED' if n==1 else 'NOT_SUPPORTED');crossok&=r.external_decision==expected
add('cross_decisions',crossok,crossok);add('triangulation',len(tri)==5,len(tri));manifest=read('FILE_MANIFEST_SHA256_v1.tsv');hok=all((RES/r.file).is_file() and sha(RES/r.file)==r.sha256 for r in manifest.itertuples(index=False));add('hashes',hok,len(manifest))
out=ROOT/'independent_validation_v2';out.mkdir(exist_ok=True);pd.DataFrame(C).to_csv(out/'validation_checks_v2.tsv',sep='\t',index=False,lineterminator='\n');p={'checks':len(C),'passed':sum(z['status']=='PASS' for z in C),'failed':sum(z['status']=='FAIL' for z in C)};(out/'validation_summary_v2.json').write_text(json.dumps(p,indent=2)+'\n');
if p['failed']==0:(out/'PHASE72_INDEPENDENT_VALIDATION_PASS_v2').write_text('PASS\n')
print(json.dumps(p));raise SystemExit(0 if p['failed']==0 else 1)
