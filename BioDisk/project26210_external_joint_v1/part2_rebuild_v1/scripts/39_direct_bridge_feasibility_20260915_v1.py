from pathlib import Path
import pandas as pd,anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');R=P/'part2_rebuild_v1'
O=P/'phase152_TYMP_SPP1_direct_bridge_plan_20260915_v1';O.mkdir(exist_ok=False)
src=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
a=ad.read_h5ad(src,backed='r');obs=a.obs.copy();a.file.close()
mal=obs[(obs.eligible_for_subclone_selection_v1==True)&(obs.source_tissue.astype(str).str.lower()=='tumor')]
m=pd.read_csv(R/'08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str},low_memory=False)
mac=m[(~m.candidate.isin(['1','Monocyte']))&(m.source_tissue.str.lower()=='tumor')]
keys=['study','curated_patient_id','sample_id']
x=mal.groupby(keys,observed=True).size().rename('malignant_cells').reset_index()
y=mac.groupby(keys,observed=True).size().rename('macrophage_cells').reset_index()
z=x.merge(y,on=keys,how='outer').fillna({'malignant_cells':0,'macrophage_cells':0})
z['eligible_same_library']=(z.malignant_cells>=20)&(z.macrophage_cells>=50)
z.to_csv(O/'SAME_LIBRARY_FEASIBILITY.tsv',sep='\t',index=False)
(O/'FEASIBILITY_STATUS.json').write_text(json.dumps(dict(status='METADATA_FEASIBILITY_COMPLETE',eligible_libraries=int(z.eligible_same_library.sum()),eligible_patients=int(z[z.eligible_same_library].curated_patient_id.nunique()),source=str(src),input_matrix_read=False,model_rerun=False,figures_created=False),indent=2))
print('COUNTS',len(mal),len(mac),'eligible_libraries',int(z.eligible_same_library.sum()),'eligible_patients',z[z.eligible_same_library].curated_patient_id.nunique())
print(z[z.eligible_same_library].groupby('study').agg(libraries=('sample_id','nunique'),patients=('curated_patient_id','nunique')).to_string())
old=pd.read_csv(P/'phase151_TYMP_SPP1_threeTF_role_20260915_v1/TUMOUR_TYMP_MACROPHAGE_MATCHED_PATIENTS.tsv',sep='\t')
print('PREVIOUS_17',old.groupby('study').size().to_string())
print('MALIGNANT_METADATA_COLUMNS',mal.columns.tolist())
for base in [P/'phase107_evidence_escalation_audit_local_v1',P/'phase110_targeted_supplement_planning_local_v1']:
 print('REFERENCE',base)
 for f in base.glob('*.tsv'):print(f.name)

