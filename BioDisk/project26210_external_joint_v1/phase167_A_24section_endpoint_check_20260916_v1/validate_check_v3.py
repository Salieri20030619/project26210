from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.sparse import load_npz
from scipy.sparse.csgraph import connected_components
O=Path(__file__).resolve().parent;ROOT=O.parent
d=pd.read_csv(O/'DESIGN_ENDPOINT_ELIGIBILITY.tsv',sep='\t');g=pd.read_csv(O/'GRAPH_QC.tsv',sep='\t')
checks=[];comps=[]
def ck(s,k,v):
 checks.append(dict(section=s,check=k,passed=bool(v)))
 if not v:raise RuntimeError(str((s,k)))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
ck('ALL','complete_unique_576_grid',len(d)==576 and not d.duplicated(['section','role','tier','candidate','transform']).any())
ck('ALL','24_sections',g.section.nunique()==24 and d.section.nunique()==24)
ck('ALL','eligibility_recomputed',np.array_equal(d.eligible,(d.spots>=50)&(d['rank']==d['columns'])&(d.condition_number<30)&(d.outcome_variance>0)))
old=ROOT/'phase166_A_technical_pilot_20260915_v1/strict_region_v2'
for s in g.itertuples(index=False):
 A=load_npz(O/(s.section+'_adjacency.npz'));q=pd.read_csv(O/(s.section+'_SPOT_TECHNICAL_QC.tsv.gz'),sep='\t')
 ck(s.section,'saved_degree',np.array_equal(np.asarray(A.sum(axis=1)).ravel(),q.degree))
 ck(s.section,'symmetric_loopfree',(A!=A.T).nnz==0 and not A.diagonal().any())
 eligible=q.source_region.to_numpy()&(q.eligible_macrophage_neighbors.to_numpy()>=3)
 ck(s.section,'receiver_membership',np.array_equal(eligible,q.eligible_receiver_link))
 ck(s.section,'region_counts',q.source_region.sum()==s.source_spots and eligible.sum()==s.receiver_link_spots)
 nc,lab=connected_components(A)
 ck(s.section,'components',nc==s.components and np.array_equal(lab,q.component))
 for k in range(nc):
  m=lab==k
  comps.append(dict(section=s.section,preparation=s.preparation,component=k,spots=int(m.sum()),source_spots=int(q.loc[m,'source_region'].sum()),receiver_link_spots=int(q.loc[m,'eligible_receiver_link'].sum())))
 if (old/(s.section+'_adjacency.npz')).exists():
  ck(s.section,'pilot_graph_exact',(A!=load_npz(old/(s.section+'_adjacency.npz'))).nnz==0)
  ck(s.section,'pilot_spot_qc_exact',q.equals(pd.read_csv(old/(s.section+'_SPOT_TECHNICAL_QC.tsv.gz'),sep='\t')))
pilot=pd.read_csv(old/'DESIGN_IDENTIFIABILITY.tsv',sep='\t')
new=d[d['transform']=='linear'].merge(pilot,on=['section','role','tier','candidate'],suffixes=('_new','_old'))
ck('PILOT','designs_reproduced_float32_aware_rtol1e5',len(new)==24 and np.array_equal(new['rank_new'],new['rank_old']) and np.allclose(new.condition_number_new,new.condition_number_old,rtol=1e-5,atol=1e-8) and np.array_equal(new.passes_frozen_design_rule,new.passes_pilot_rule))
for x in pd.read_csv(O/'INPUT_MANIFEST.tsv',sep='\t').itertuples(index=False):
 p=Path(x.path)
 ck(p.name,'input_size_mtime_unchanged',p.stat().st_size==x.bytes and p.stat().st_mtime_ns==x.mtime_ns)
 ck(p.name,'input_sha256_unchanged',sha(p)==x.sha256)
summary=d.groupby(['section','preparation','role']).eligible.all().unstack('role').reset_index().rename(columns={'source':'source_all_designs_eligible','receiver_link':'receiver_all_designs_eligible'})
summary=summary.merge(g[['section','source_spots','receiver_link_spots','components','largest_component_fraction']],on='section')
summary.to_csv(O/'SECTION_ENDPOINT_SUMMARY.tsv',sep='\t',index=False)
pd.DataFrame(comps).to_csv(O/'CONNECTED_COMPONENT_QC.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(O/'RELOAD_VALIDATION.tsv',sep='\t',index=False)
print(summary.to_string(index=False))
print('VALIDATION',len(checks),sum(x['passed'] for x in checks))
print('MAX_COND',d.loc[d.eligible].groupby(['role','transform']).condition_number.max().to_dict())
print(pd.DataFrame(comps).query("section=='GSM5924048'").to_string(index=False))
