from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import rankdata
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase187_IL6_project_receiver_validation_20260917_v1'
checks=[]
def qr(x,y,cv):
 C=np.column_stack([np.ones(len(x))]+[rankdata(np.round(v,12)) for v in cv]);Q,R=np.linalg.qr(C)
 assert np.linalg.matrix_rank(C)==C.shape[1]
 a=rankdata(np.round(x,12));b=rankdata(np.round(y,12));a=a-Q@(Q.T@a);b=b-Q@(Q.T@b)
 return np.corrcoef(a,b)[0,1]
for row in pd.read_csv(d/'PUBLIC_COHORT_EDGES.tsv',sep='\t').to_dict('records'):
 f=pd.read_csv(d/f"PUBLIC_PATIENTS_{row['cohort']}.tsv",sep='\t',index_col=0)
 cv=['MTAP','MYELOID_PROXY']+(['HYPOXIA_PROXY'] if row['adjustment']=='hypoxia' else [])
 x,y=('MET:5-methylthioadenosine (MTA)','IL6') if row['edge']=='MTA_IL6' else ('IL6','SPP1')
 if row['edge']=='IL6_SPP1':cv+=['MET:5-methylthioadenosine (MTA)']
 rho=qr(f[x].to_numpy(),f[y].to_numpy(),[f[c].to_numpy() for c in cv]);checks.append({'type':'bulk_QR','cohort':row['cohort'],'edge':row['edge'],'adjustment':row['adjustment'],'passed':abs(rho-row['rho'])<1e-10})
for row in pd.read_csv(d/'IL6_STAT3_PUBLIC_COHORT.tsv',sep='\t').to_dict('records'):
 f=pd.read_csv(d/f"PUBLIC_PATIENTS_{row['cohort']}.tsv",sep='\t',index_col=0)
 score=pd.read_csv(r/'phase179_STAT3_external_activity_validation_20260916_v1/run_v2'/f"BULK_PATIENTS_{row['cohort']}_{row['resource']}_zmean.tsv",sep='\t',index_col=0).loc[f.index,'activity']
 cv=['MTAP','MYELOID_PROXY','MET:5-methylthioadenosine (MTA)']+(['HYPOXIA_PROXY'] if row['adjustment']=='hypoxia' else [])
 rho=qr(f.IL6.to_numpy(),score.to_numpy(),[f[c].to_numpy() for c in cv]);checks.append({'type':'IL6_STAT3_QR','cohort':row['cohort'],'edge':row['resource'],'adjustment':row['adjustment'],'passed':abs(rho-row['rho'])<1e-10})
assert all(c['passed'] for c in checks)
pd.DataFrame(checks).to_csv(d/'INDEPENDENT_QR_VALIDATION.tsv',sep='\t',index=False)
# Patient-independent source qualification recomputation.
o=r/'phase182_IL10_LGALS3_receiver_audit_20260916_v1'
m=pd.read_csv(o/'SOURCE_CELL_METADATA.tsv.gz',sep='\t',index_col=0);z=np.load(o/'SOURCE_PANEL_CACHE.npz')
mask=(m.source_class=='Fibroblast')&(m.source_tissue.str.lower()=='tumor')
fib=m[mask];counts=z['counts'][mask.to_numpy(),list(z['genes']).index('IL6')]
q={'all_cells':len(m),'libraries':int(m.sample_id.nunique()),'tumor_fibroblast_cells':len(fib),'tumor_fibroblast_patients':int(fib.curated_patient_id.nunique()),'tumor_fibroblast_IL6_counts':float(counts.sum()),'tumor_fibroblast_IL6_detected_cells':int(sum(counts>0)),'inference':'not_estimable'}
(d/'FIBROBLAST_SOURCE_GATE.json').write_text(json.dumps(q,indent=2))
pd.DataFrame({'file':[str(p) for p in [o/'SOURCE_PANEL_CACHE.npz',o/'SOURCE_CELL_METADATA.tsv.gz',o/'SCORE_EDGE_CACHE.npz']],'sha256':[hashlib.sha256(p.read_bytes()).hexdigest() for p in [o/'SOURCE_PANEL_CACHE.npz',o/'SOURCE_CELL_METADATA.tsv.gz',o/'SCORE_EDGE_CACHE.npz']]}).to_csv(d/'INPUT_HASHES.tsv',sep='\t',index=False)
print('QR',len(checks),'PASS',q)
