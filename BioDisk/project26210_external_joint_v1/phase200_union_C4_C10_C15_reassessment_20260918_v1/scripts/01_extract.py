from pathlib import Path
import json,hashlib,datetime,subprocess
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1')
d=r/'phase200_union_C4_C10_C15_reassessment_20260918_v1'
src=r/'part2_rebuild_v1'
gate=subprocess.run(['bash','/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh'],capture_output=True,text=True)
(d/'qc/MEMORY_GATE.log').write_text(gate.stdout+gate.stderr);assert gate.returncode==0
mp=r/'phase192_sections28_212_figures_20260917_v1/inputs/MACRO_DISPLAY_METADATA.tsv.gz'
m=pd.read_csv(mp,sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
assert len(m)==22444 and m.index.is_unique
m['is_union']=m.candidate.isin(['4','10','15'])
assert m.is_union.sum()==4502
dominant=m.loc[m.candidate=='10','curated_patient_id'].value_counts().index[0]
h=src/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad'
a=ad.read_h5ad(h,backed='r')
positions=a.obs_names.get_indexer(m.index);assert (positions>=0).all()
order=np.argsort(positions);x=a.X[positions[order]].tocsr()[np.argsort(order)].astype(np.float64)
genes=a.var_names.astype(str).to_numpy();a.file.close()
lib=np.asarray(x.sum(1)).ravel()
assert np.allclose(lib,m.total_counts_raw.to_numpy(),atol=1e-5,rtol=0)
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
m['library_counts']=lib
m.to_csv(d/'inputs/MACRO_METADATA.tsv.gz',sep='\t',index_label='cell_id')
sparse.save_npz(d/'inputs/MACRO_COUNTS.npz',x)
pd.DataFrame({'gene':genes}).to_csv(d/'inputs/GENES.tsv',sep='\t',index=False)
keys=['curated_patient_id','study','sample_id','source_tissue','candidate']
m['group_key']=m[keys].astype(str).agg('|'.join,axis=1)
codes,labels=pd.factorize(m.group_key,sort=True)
G=sparse.csr_matrix((np.ones(len(m)),(codes,np.arange(len(m)))),shape=(len(labels),len(m)))
pb=(G@x).tocsr();det=(G@(x>0)).tocsr()
norm=x.multiply(10000/lib[:,None]).tocsr();norm.data=np.log1p(norm.data)
meansum=(G@norm).tocsr()
group=m.groupby('group_key',sort=True).agg(patient=('curated_patient_id','first'),study=('study','first'),sample_id=('sample_id','first'),tissue=('source_tissue','first'),cluster=('candidate','first'),cells=('candidate','size')).reset_index()
assert group.group_key.tolist()==list(labels)
group.to_csv(d/'inputs/GROUP_INDEX.tsv',sep='\t',index=False)
sparse.save_npz(d/'inputs/GROUP_COUNTS.npz',pb);sparse.save_npz(d/'inputs/GROUP_DETECTED.npz',det);sparse.save_npz(d/'inputs/GROUP_LOGSUM.npz',meansum)
old=src/'04_markers_annotation/run_v4_macrophage_only'
og=pd.read_csv(old/'PSEUDOBULK_GENE_INDEX.tsv',sep='\t').gene.to_numpy()
assert np.array_equal(genes,og)
opb=sparse.load_npz(old/'PATIENT_CLUSTER_PSEUDOBULK_COUNTS.npz')
ori=pd.read_csv(old/'PSEUDOBULK_ROW_INDEX.tsv',sep='\t',dtype={'cluster':str})
for i,row in ori.iterrows():
 ix=(group.patient==row.patient)&(group.cluster==row.cluster)
 got=np.asarray(pb[ix.to_numpy()].sum(0)).ravel()
 assert np.array_equal(got,opb[i].toarray().ravel())
 assert group.loc[ix,'cells'].sum()==row.n_cells
composition=[]
for name,cs in {'U':['4','10','15'],'C4':['4'],'C10':['10'],'C15':['15']}.items():
 z=m[m.candidate.isin(cs)];pc=z.curated_patient_id.value_counts();sc=z.study.value_counts()
 composition.append(dict(group=name,cells=len(z),patients=len(pc),studies=len(sc),dominant_patient_fraction=pc.iloc[0]/len(z),dominant_study_fraction=sc.iloc[0]/len(z)))
pd.DataFrame(composition).to_csv(d/'core/COMPOSITION.tsv',sep='\t',index=False)
(d/'inputs/FROZEN_UNION.json').write_text(json.dumps({'members':['4','10','15'],'cells':4502,'all_mature_macrophages':22444,'dominant_C10_patient':str(dominant),'source_metadata':str(mp),'metadata_sha256':hashlib.sha256(mp.read_bytes()).hexdigest(),'raw_source':str(h),'raw_source_bytes':h.stat().st_size,'full_old_patient_cluster_counts_equal':len(ori),'created':datetime.datetime.now().isoformat()},indent=2))
print('EXTRACTION_VALIDATED',x.shape,'groups',len(group),'old_count_rows',len(ori),flush=True)
