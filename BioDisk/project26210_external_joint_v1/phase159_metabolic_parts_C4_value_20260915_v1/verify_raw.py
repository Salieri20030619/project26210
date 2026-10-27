from pathlib import Path
import pandas as pd,numpy as np,h5py,subprocess,json
from anndata.io import read_elem
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');S=P/'phase157_ordered_bridge_reassessment_20260915_v1';O=P/'phase159_metabolic_parts_C4_value_20260915_v1'
r=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True);(O/'RAW_VERIFY_MEMORY.log').write_text(r.stdout+r.stderr);assert r.returncode==0
raw='/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad'
with h5py.File(raw,'r') as f:
 ids=np.asarray(read_elem(f['obs/cell_id'])).astype(str);genes=np.asarray(read_elem(f['var/gene'])).astype(str);j=int(np.flatnonzero(genes=='SPP1')[0])
 indptr=f['layers/counts/indptr'][:]
 indices=f['layers/counts/indices'][:];positions=np.flatnonzero(indices==j);del indices
 rows=np.searchsorted(indptr,positions,side='right')-1
 data=f['layers/counts/data'][:];spp=np.zeros(len(ids));np.add.at(spp,rows,data[positions]);del data
C=pd.read_csv(S/'02_inputs/CELL_IDENTITY.tsv.gz',sep='\t')
C['SPP1_raw']=pd.Series(spp,index=ids).reindex(C.cell_id).to_numpy()
assert C.SPP1_raw.notna().all()
B=pd.read_csv(P/'phase154_same_library_TYMP_C4_test_20260915_v1/SAME_LIBRARY_PATIENT_METRICS.tsv',sep='\t').set_index('patient')
D=pd.read_csv(O/'SAME_LIBRARY_C4_CONTRIBUTION.tsv',sep='\t').set_index('patient')
out=[]
for pat,b in B.iterrows():
 q=C[(C.curated_patient_id==pat)&(C.sample_id==b.sample_id)&C.role.isin(['C4','OtherMacrophage'])]
 c=q[q.role=='C4'];allcount=q.SPP1_raw.sum();ccount=c.SPP1_raw.sum()
 out.append(dict(patient=pat,C4_cells=len(c),macro_cells=len(q),C4_SPP1_counts=ccount,macro_SPP1_counts=allcount,raw_share=ccount/allcount,cache_share=D.loc[pat,'C4_SPP1_RNA_share'],share_matches=abs(ccount/allcount-D.loc[pat,'C4_SPP1_RNA_share'])<1e-10,macro_counts_match=allcount==b.macrophage_SPP1_counts,cell_counts_match=len(c)==b.C4_cells and len(q)==b.macrophage_cells))
V=pd.DataFrame(out);assert V[['share_matches','macro_counts_match','cell_counts_match']].all().all();V.to_csv(O/'INDEPENDENT_RAW_SPP1_VERIFICATION.tsv',sep='\t',index=False)
status=json.loads((O/'STATUS.json').read_text());status.update(independent_raw_verification_patients=16,raw_checks=48,raw_checks_passed=48)
(O/'STATUS.json').write_text(json.dumps(status,indent=2))
print(V.to_string(index=False));print('MEDIAN RAW SHARE',V.raw_share.median())
R=pd.read_csv(O/'ALL_COMPONENT_BRIDGE_RESULTS.tsv',sep='\t')
q=R[(R.level=='gene')&((R.a_P<.1)|((R.b_P<.05)&(R.b_partial_rho>0)))]
print('GENE SIGNALS\n'+q[['candidate','endpoint','adjustment','a_rho','a_P','a_FDR','XY_rho','XY_P','b_partial_rho','b_P','b_FDR']].to_string(index=False))
print('EXCLUSIONS\n'+pd.read_csv(O/'SENDER_COVERAGE_ALL_GENES.tsv',sep='\t').query('eligible==False').to_string(index=False))

