from pathlib import Path
import json,datetime
import numpy as np,pandas as pd,anndata as ad
from scipy.io import mmwrite
import gzip
R=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/part2_rebuild_v1');O=R/'07_trajectory/input_v1'
assert json.loads((R/'06_functional_programs/run_v1/STATUS.json').read_text())['status']=='COMPLETE_FUNCTIONAL_NUMERICAL'
O.mkdir(exist_ok=False)
a=ad.read_h5ad(R/'01_input_audit/run_v1/myeloid_raw_counts_v1.h5ad');cl=pd.read_csv(R/'03_macrophage_clustering/run_v4_exploratory/CELL_CLUSTER_ASSIGNMENTS.tsv.gz',sep='\t',index_col=0,dtype={'cluster':str},low_memory=False)
macro=cl[~cl.cluster.isin(['2','13','14','16'])].copy();mono=a.obs[a.obs.cell_type_broad_v1.astype(str)=='Monocyte'].copy();assert len(mono)==8215 and len(macro)==25416
ids=list(macro.index)+list(mono.index);assert len(set(ids))==len(ids)
b=a[ids].copy();m=b.obs.copy();m['candidate']=pd.Series(macro.cluster,index=macro.index).reindex(m.index).fillna('Monocyte');m['origin_label']=np.where(m.candidate=='Monocyte','original_Monocyte','purified_macrophage_catalogue');x=b.X.tocsr();tot=np.asarray(x.sum(axis=1)).ravel();norm=x.astype(float).multiply(10000/tot[:,None]).tocsr();norm.data=np.log1p(norm.data)
panels=json.loads((R/'02_myeloid_qc/run_v2/STATUS.json').read_text())['panels'];panels['DC']=['CD1C','FCER1A','CLEC10A','CD1E'];panels['classical_mono']=['FCN1','S100A8','S100A9','VCAN','LYZ']
for name,gs in panels.items():
 idx=b.var_names.get_indexer([g for g in gs if g in b.var_names]);m[name+'_score']=np.asarray(norm[:,idx].mean(axis=1)).ravel();m[name+'_detected']=np.asarray((x[:,idx]>0).sum(axis=1)).ravel()
m['low_complexity']=(m.n_genes_raw<500)|(m.total_counts_raw<1000);m['nonmyeloid_dominant']=False
for name in panels:
 if name not in ['myeloid_backbone','classical_mono']:m['nonmyeloid_dominant']|=(m[name+'_detected']>=3)&(m[name+'_score']>m.myeloid_backbone_score)
# C4 and existing mature catalogue membership stay frozen; screening additional mono/C1 boundary inputs.
m['additional_input']=m.candidate.isin(['Monocyte','1']);m['keep_primary']=~(m.additional_input&(m.low_complexity|m.nonmyeloid_dominant|(m.myeloid_backbone_detected<3)))
m['root_eligible']=m.keep_primary&(m.candidate=='Monocyte')&(m.classical_mono_detected>=2)&(m.myeloid_backbone_detected>=3)&~m.nonmyeloid_dominant
m['is_C4']=m.candidate=='4';m['is_cycling_or_stress']=m.candidate.isin(['7','11','15']);m.to_csv(O/'CELL_INPUT_AUDIT.tsv.gz',sep='\t',index_label='cell_id')
kept=m[m.keep_primary].copy();assert (kept.candidate=='4').sum()==2674;b=b[kept.index].copy();b.obs=kept;b.write_h5ad(O/'TRAJECTORY_RAW_COUNTS.h5ad',compression='gzip')
kept.to_csv(O/'CELL_METADATA.tsv.gz',sep='\t',index_label='cell_id')
for field in ['study','curated_patient_id','sample_id','source_tissue']:pd.crosstab(kept[field],kept.candidate).to_csv(O/('COUNTS_BY_'+field+'.tsv'),sep='\t')
qc=kept.groupby('study',observed=True).agg(cells=('candidate','size'),patients=('curated_patient_id','nunique'),libraries=('sample_id','nunique'),root_cells=('root_eligible','sum'),C4_cells=('is_C4','sum'));qc.to_csv(O/'STUDY_INPUT_QC.tsv',sep='\t')
stat={'status':'COMPLETE_INPUT_QC','source_candidates':len(m),'retained_cells':len(kept),'excluded_additional_inputs':int((~m.keep_primary).sum()),'patients':kept.curated_patient_id.nunique(),'libraries':kept.sample_id.nunique(),'C4_frozen_cells':2674,'root_eligible':int(kept.root_eligible.sum()),'figures_created':False,'primary_graph':'Per study; original Monocyte roots within each partition; rootless partitions remain unrooted. C1 is boundary retained but not automatic root.','panels':panels};(O/'STATUS.json').write_text(json.dumps(stat,indent=2));print(json.dumps(stat),flush=True);print(qc.to_csv(sep='\t'),flush=True)
