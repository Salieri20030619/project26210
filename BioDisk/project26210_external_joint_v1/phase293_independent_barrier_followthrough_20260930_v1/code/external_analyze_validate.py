from pathlib import Path
import json,hashlib,re,subprocess
import pandas as pd,numpy as np
from scipy import stats
P=Path(__file__).resolve().parents[1]
q=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(P/'logs/external_stats_memory.log').write_text(q.stdout+q.stderr);assert q.returncode==0
f=json.loads((P/'qc/SOURCEPAIR_FREEZE.json').read_text())
for n,h in f['files'].items():assert hashlib.sha256((P/n).read_bytes()).hexdigest()==h
d=pd.read_csv(P/'tables/EXTERNAL_SOURCEPAIR_FREEZE.tsv',sep='\t')
out=[]
for key,z in d.groupby(['dataset','stratum','receiver','source','endpoint']):
 assert z.patient.is_unique
 n=len(z);row=dict(zip(['dataset','stratum','receiver','source','endpoint'],key))
 for model in ['unadjusted','receiver_depth']:
  D=np.ones((n,1)) if model=='unadjusted' else np.column_stack([np.ones(n),stats.rankdata(z.receiver_depth)])
  rank=np.linalg.matrix_rank(D);rho=np.nan
  if n>=3 and n-rank>=2 and rank==D.shape[1]:
   a=stats.rankdata(z.x);b=stats.rankdata(z.y)
   a-=D@np.linalg.lstsq(D,a,rcond=None)[0];b-=D@np.linalg.lstsq(D,b,rcond=None)[0]
   den=np.linalg.norm(a)*np.linalg.norm(b)
   if den>1e-10:rho=float(a@b/den)
  assert n<10,'Unexpected inferential sample size: implement prespecified permutation before proceeding'
  out.append(dict(**row,model=model,n=n,rho=rho,P=np.nan,FDR=np.nan,status='EXPLORATORY_DESCRIPTIVE_N_LT10' if np.isfinite(rho) else 'NOT_ESTIMABLE_N_OR_DF',residual_df=n-rank-1))
a=pd.DataFrame(out);a.to_csv(P/'tables/EXTERNAL_ASSOCIATIONS.tsv',sep='\t',index=False)
# Independently reconstruct qualification predicates and aggregate patient frozen inputs.
lib=pd.read_csv(P/'tables/EXTERNAL_LIBRARY_INPUTS.tsv',sep='\t',float_precision='round_trip');checks=[]
def chk(n,v):checks.append(dict(check=n,passed=bool(v)));assert v,n
for row in d.itertuples():
 z=lib[(lib.dataset==row.dataset)&(lib.patient==row.patient)&(lib.stratum==row.stratum)&(lib.receiver==row.receiver)&lib.library.isin(row.libraries.split(';'))]
 chk('library_count_'+str(len(checks)),len(z)==row.n_libraries)
 cnt='tumor_like_cells' if row.source=='tumor_like_TYMP' else 'macro_cells';cut=50 if row.source=='tumor_like_TYMP' else 100
 chk('source_gate_'+str(len(checks)),z[cnt].ge(cut).all())
 chk('receiver_gate_'+str(len(checks)),z.receiver_cells.ge(20 if row.receiver=='strict_fibro' else 50).all())
 chk('x_mean_'+str(len(checks)),np.isclose(z[row.source].mean(),row.x,atol=1e-12))
 chk('y_mean_'+str(len(checks)),np.isclose(z[row.endpoint].mean(),row.y,atol=1e-12))
for row in a.itertuples():
 z=d[(d.dataset==row.dataset)&(d.stratum==row.stratum)&(d.receiver==row.receiver)&(d.source==row.source)&(d.endpoint==row.endpoint)]
 if row.model=='unadjusted' and np.isfinite(row.rho):
  rr=z[['x','y']].rank().corr().iloc[0,1];chk('independent_rho_'+str(len(checks)),abs(rr-row.rho)<1e-12)
 elif np.isfinite(row.rho):
  xyz=z[['x','y','receiver_depth']].rank().corr();rxy=xyz.iloc[0,1];rxz=xyz.iloc[0,2];ryz=xyz.iloc[1,2]
  rr=(rxy-rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))
  chk('independent_partial_'+str(len(checks)),abs(rr-row.rho)<1e-12)
chk('no_inference_n_lt10',a.P.isna().all() and a.FDR.isna().all())
# Source metadata extraction excludes expression tables; preserve public sample identity evidence.
geo=P/'resources/GSE264586_metadata.txt'
if geo.exists():
 rows=[]
 for block in re.split(r'\^SAMPLE = ',geo.read_text(errors='replace'))[1:]:
  rec=dict(GSM=block.splitlines()[0])
  for k in ['title','source_name_ch1','characteristics_ch1']:
   rec[k]=';'.join(v.split(' = ',1)[1] for v in block.splitlines() if v.startswith('!Sample_'+k+' = '))
  rows.append(rec)
 pd.DataFrame(rows).to_csv(P/'tables/GSE264586_GEO_SAMPLE_IDENTITY.tsv',sep='\t',index=False)
 print(pd.DataFrame(rows)[['GSM','title','source_name_ch1','characteristics_ch1']].head(8).to_string(index=False))
pd.DataFrame(checks).to_csv(P/'qc/EXTERNAL_STATS_CHECKS.tsv',sep='\t',index=False)
(P/'qc/EXTERNAL_STATS_VALIDATED.json').write_text(json.dumps(dict(status='PASS',checks=len(checks),rows=len(a),n_max=int(a.n.max()),inferential_models=0),indent=2))
print(a[(a.receiver=='strict_fibro')&a.model.eq('unadjusted')].to_string(index=False))
print('EXTERNAL_DESCRIPTIVE_COMPLETE')
