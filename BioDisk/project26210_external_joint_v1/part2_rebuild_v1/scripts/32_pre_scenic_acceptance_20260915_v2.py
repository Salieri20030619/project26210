from pathlib import Path
import json,hashlib,datetime,sys
import numpy as np,pandas as pd
from scipy.stats import wilcoxon
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1')
O=R/'provenance/pre_scenic_acceptance_20260915_v2'
O.mkdir(exist_ok=False)
# Reuse only upstream root/identity checks from existing validator; no writes to old results.
src=(R/'scripts/29_validate_remaining_results_v2.py').read_text()
start=src.index("m=pd.read_csv")
end=src.index("for filename,keys,qcol")
checks=[]
def ck(name,value): checks.append({'check':name,'pass':bool(value)})
exec(src[start:end])
def close(a,b): return np.allclose(a,b,atol=1e-10,rtol=1e-7,equal_nan=True)
def bh(v):
 from scipy.stats import false_discovery_control
 v=np.asarray(v,float); out=np.full(len(v),np.nan); mask=np.isfinite(v)
 if mask.any():out[mask]=false_discovery_control(v[mask],method='bh')
 return out
# Raw pseudobulk reaggregation from the independent functional validator, redirect outputs.
s=(R/'scripts/19b_functional_review_validation_json_fixed.py').read_text()
s=s[s.index("genes=pd.read_csv"):s.index("assert all(x['pass']")]
# The original resets checks; preserve earlier checks.
old=checks.copy()
F=R/'06_functional_programs/run_v1'
env=dict(globals());env.update(O=F,M=R/'04_markers_annotation/run_v4_macrophage_only')
from scipy import sparse
from scipy.stats import hypergeom
env.update(sparse=sparse,hypergeom=hypergeom)
exec(s,env); checks=old+env['checks']
for fn,keys in [('ALL_PATIENT_PROGRAM_COMPARISONS.tsv.gz',['cluster','library','version']),('ALL_ORA_RESULTS.tsv.gz',['cluster','library','version'])]:
 df=pd.read_csv(F/fn,sep='\t')
 qcols=[c for c in df if 'FDR' in c]
 ck(fn+'_has_FDR',len(qcols)==1)
 if len(qcols)==1:
  for fam,g in df.groupby(keys):ck(fn+str(fam)+'_BH',close(bh(g.P),g[qcols[0]]))
D=R/'07_trajectory/dynamics_v2'
for rawfile,patientfile,keys,values in [
 ('PARTITION_PATIENT_PAIRED_PSEUDOTIME.tsv','PATIENT_PAIRED_PSEUDOTIME.tsv',['mode','scope','pt_version','study','patient','endpoint'],['effect']),
 ('PARTITION_PATIENT_DYNAMIC_TRENDS.tsv.gz','PATIENT_DYNAMIC_TRENDS.tsv',['mode','scope','pt_version','study','patient','family','feature'],['rho','depth_adjusted_rho'])]:
 raw=pd.read_csv(D/rawfile,sep='\t'); pat=pd.read_csv(D/patientfile,sep='\t')
 agg=raw.groupby(keys)[values].median().sort_index(); given=pat.set_index(keys)[values].sort_index()
 ck(patientfile+'_unique',given.index.is_unique);ck(patientfile+'_partition_collapse',agg.index.equals(given.index) and close(agg.values,given.values))
for sf,pf,keys,values,family in [
 ('PSEUDOTIME_SUMMARY.tsv','PATIENT_PAIRED_PSEUDOTIME.tsv',['mode','scope','pt_version','endpoint'],['effect'],['mode','scope','pt_version']),
 ('DYNAMIC_TREND_SUMMARY.tsv','PATIENT_DYNAMIC_TRENDS.tsv',['mode','scope','pt_version','family','feature'],['rho','depth_adjusted_rho'],['mode','scope','pt_version','family','adjustment'])]:
 summary=pd.read_csv(D/sf,sep='\t');pat=pd.read_csv(D/pf,sep='\t')
 for key,g in pat.groupby(keys):
  subset=summary.copy()
  for k,v in zip(keys,key):subset=subset[subset[k]==v]
  for value in values:
   row=subset[subset.adjustment==value].iloc[0] if len(values)>1 else subset.iloc[0]
   a=g[value].to_numpy();a=a[np.isfinite(a)]
   p=wilcoxon(a,method='approx',correction=True).pvalue if len(a)>=3 and np.any(a) else (1. if len(a)>=3 else np.nan)
   label=sf+str(key)+value
   ck(label+'_summary',len(a)==row.n_patients and (a>0).sum()==row.positive and (a<0).sum()==row.negative and close(np.median(a) if len(a) else np.nan,row.median_effect) and close(p,row.P))
   ck(label+'_small_n',row.n_patients>=10 or row.reporting_tier in ['small_n_exploratory','descriptive'])
 for fam,g in summary.groupby(family):ck(sf+str(fam)+'_BH',close(bh(g.P),g.FDR))
# Independently recompute every eligible same-partition C4 pseudotime contrast.
raw=pd.read_csv(D/'PARTITION_PATIENT_PAIRED_PSEUDOTIME.tsv',sep='\t',dtype={'partition':str})
for i,row in raw.iterrows():
 z=t[(t['mode']==row['mode'])&(t.study==row.study)&(t.curated_patient_id==row.patient)&(t.partition==row.partition)]
 if row.scope=='tumor_only':z=z[z.source_tissue.str.lower()=='tumor']
 z=z[np.isfinite(z[row.pt_version])]
 a=z.loc[z.candidate=='4',row.pt_version]
 b=z.loc[z.candidate=='Monocyte',row.pt_version] if row.endpoint=='C4_minus_Monocyte' else z.loc[~z.candidate.isin(['4','1','Monocyte']),row.pt_version]
 ck('partition_contrast_'+str(i),len(a)>=20 and len(b)>=20 and len(a)==row.C4_cells and len(b)==row.comparator_cells and close(a.median()-b.median(),row.effect))
for mode,g in t.groupby('mode'):
 expected=m.index if mode=='primary' else m.index[~m.candidate.isin(['7','11','15'])]
 ck(mode+'_exact_cells',g.index.is_unique and set(g.index)==set(expected))
ck('ten_models',t.groupby(['mode','study']).ngroups==10)
coverage=pd.read_csv(D/'PROGRAM_COVERAGE.tsv',sep='\t')
ck('programs_exclude_endpoints',all(not {'SPP1','TYMP'}&set(s.split(';')) for s in coverage.gene_symbols))
for directory in ['06_functional_programs','07_trajectory','08_regulons']:
 ck(directory+'_no_new_figures',not any(p.suffix.lower() in ['.png','.pdf','.svg','.jpg'] for p in (R/directory).rglob('*')))
(O/'CHECKS.json').write_text(json.dumps(checks,indent=2))
status={'status':'PASS' if all(x['pass'] for x in checks) else 'FAIL','checks':len(checks),'passed':sum(x['pass'] for x in checks),'failures':[x for x in checks if not x['pass']],'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Upstream functional, trajectory and patient statistics; excludes incomplete SCENIC','python':sys.version,'figures_created':False}
(O/'STATUS.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
assert status['status']=='PASS'

