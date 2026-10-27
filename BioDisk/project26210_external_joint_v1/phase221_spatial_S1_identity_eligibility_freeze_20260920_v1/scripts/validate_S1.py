from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import sparse
P=Path(__file__).resolve().parents[1];R=P.parent;B=R.parent
assert (P/'EXIT_STATUS').read_text().strip()=='0'
e=pd.read_csv(P/'S2_SECTION_ELIGIBILITY.tsv',sep='\t');patients=pd.read_csv(P/'S2_PATIENT_ELIGIBILITY.tsv',sep='\t');mp=pd.read_csv(P/'FROZEN_SECTION_PATIENT_MAP.tsv',sep='\t')
assert len(e)==48 and e.section.nunique()==24 and mp.patient_id.nunique()==23 and len(patients)==46
assert e.eligible.all() and (e.source_spots>=50).all() and (e.target_supported_spots>=20).all() and (e.design_rank==7).all() and (e.condition<30).all()
checks=[]
for sid in mp.section:
 f=P/'frozen_inputs'/(sid+'_SPOT_INPUTS.tsv.gz');x=pd.read_csv(f,sep='\t',index_col='spot_id');a=pd.read_csv(B/'project26210_phase78_fine_state_spatial_v1/run_v1/sections_v1'/sid/'fine_state_abundance_v1.tsv.gz',sep='\t',index_col='spot_id');assert a.index.equals(x.index)
 macro=a['q05__SPP1-state Macrophage']+a['q05__C1QC-state Macrophage']+a['q05__Other Macrophage'];ratio=a['q05__SPP1-state Macrophage']/macro
 assert np.allclose(ratio,x.SPP1_share_within_mature_macro,rtol=1e-12,atol=1e-12)
 src=(x.total_macro_fraction>=.1)&(x.degree>=3);assert np.array_equal(src,x.source_eligible)
 G=sparse.load_npz(P/'frozen_inputs'/(sid+'_adjacency.npz')).tocsr();inds=np.flatnonzero(src);select=inds[[0,len(inds)//2,-1]]
 for i in select:
  nb=G.indices[G.indptr[i]:G.indptr[i+1]]
  for lab in ['endothelial','pericyte']:
   assert abs(x[lab+'_fraction'].iloc[nb].mean()-x['neighbor_'+lab+'_fraction'].iloc[i])<1e-12
 assert int(src.sum())==int(e[e.section==sid].source_spots.iloc[0])
 checks.append({'section':sid,'spots':len(x),'source_spots':int(src.sum()),'independent_neighbor_checks':6})
assert sum(x['spots'] for x in checks)==72952
for endpoint,g in patients.groupby('endpoint'):
 double=g[g.n_eligible_sections==2];assert len(double)==1 and double.patient.iloc[0]=='GSE175540:BIONIKK:c_2'
 for x in g.itertuples():assert set(x.eligible_sections.split(';'))==set(e[(e.endpoint==endpoint)&(e.patient==x.patient)&e.eligible].section)
for name,h in json.loads((P/'FREEZE_SHA256.json').read_text()).items():assert hashlib.sha256((P/name).read_bytes()).hexdigest()==h
for f in json.loads((P/'INPUT_MANIFEST.json').read_text()):
 q=Path(f['path']);assert q.stat().st_size==f['bytes']
 if f['sha256']:assert hashlib.sha256(q.read_bytes()).hexdigest()==f['sha256']
for f in json.loads((P/'audit/RDS_BEFORE.json').read_text()):
 st=Path(f['path']).stat();assert (st.st_size,st.st_mtime_ns,st.st_ino)==(f['bytes'],f['mtime_ns'],f['inode'])
result={'status':'PASS','sections_checked':24,'spots':72952,'endpoint_section_rows':48,'endpoint_patient_rows':46,'independent_macrophage_share_checks':24,'independent_neighbor_checks':144,'c2_duplicate_resolved_per_endpoint':True,'frozen_hashes_match':True,'input_files_unchanged':True,'RDS_unchanged':333,'spatial_correlations_P_FDR_calculated':False,'checks':checks}
(P/'S1_FINAL_VALIDATION.json').write_text(json.dumps(result,indent=2));print('S1_VALIDATED_24_SECTIONS_23_PATIENTS')
