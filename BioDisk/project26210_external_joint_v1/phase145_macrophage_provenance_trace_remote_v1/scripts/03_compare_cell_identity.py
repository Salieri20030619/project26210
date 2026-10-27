from pathlib import Path
import pandas as pd,numpy as np,anndata as ad,json
P=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1');O=P/'identity';O.mkdir()
a=pd.read_csv(P/'legacy_object/metadata.tsv.gz',sep='\t',low_memory=False)
p='/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad'
b=ad.read_h5ad(p,backed='r');obs=b.obs.copy();obs.index=obs.index.astype(str);b.file.close()
ref=pd.read_csv('/mnt/BioDisk/wangmingyuan/project26210_part4_sensitivity_pilot_v1/reference_v2/reference_v2_cell_assignments.tsv.gz',sep='\t',low_memory=False)
print('KEYS',json.dumps({'old_id_examples':a.cell_id.head(2).tolist(),'current_id_examples':obs.index[:2].tolist(),'ref_columns':list(ref.columns),'current_columns':list(obs.columns)},ensure_ascii=False))
common=set(a.cell_id)&set(obs.index);print('EXACT_ID_INTERSECTION',len(common))
pd.crosstab(a.macro_cluster_final,a.macro_subtype).to_csv(O/'cluster_to_subtype.tsv',sep='\t')
checks={c:int((a[c].astype(str)==a.macro_cluster_final.astype(str)).sum()) for c in a if c.startswith('RNA_snn_res.')}
by_study=pd.crosstab(a.macro_subtype,a.study);by_study.to_csv(O/'legacy_subtype_by_study.tsv',sep='\t')
by_sample=pd.crosstab(a.macro_subtype,a.sample_id);by_sample.to_csv(O/'legacy_subtype_by_library.tsv',sep='\t')
# Record exact matching only; derive alternative keys only if evidence later requires it.
j=a.set_index('cell_id').join(obs,how='inner',lsuffix='_old',rsuffix='_current')
j.to_csv(O/'old_current_exact_join.tsv.gz',sep='\t',compression='gzip')
if len(j):
 for c in ['cell_type_broad_v1','curated_patient_id','source_tissue']:
  if c in j:pd.crosstab(j.macro_subtype,j[c]).to_csv(O/('legacy_subtype_by_current_'+c+'.tsv'),sep='\t')
 new_col='state_v2';r=ref.set_index('exact_cell_id');jr=j.join(r[['state_v2']],how='left')
 pd.crosstab(jr.macro_subtype,jr.state_v2.fillna('NOT_IN_REFERENCE')).to_csv(O/'legacy_subtype_by_current_reference_state.tsv',sep='\t')
 jr[['macro_subtype','macro_cluster_final','state_v2','curated_patient_id','cell_type_broad_v1']].to_csv(O/'cell_label_crosswalk.tsv.gz',sep='\t',compression='gzip')
state={'old_cells':len(a),'old_libraries':a.sample_id.nunique(),'old_studies':a.study.nunique(),'current_cells':len(obs),'exact_overlap':len(common),'old_cluster_resolution_matches':checks,'old_subtypes':a.macro_subtype.value_counts().to_dict(),'current_patients_among_join':int(j.curated_patient_id.nunique()) if 'curated_patient_id' in j else None}
(O/'IDENTITY_STATUS.json').write_text(json.dumps(state,ensure_ascii=False,indent=2));print('STATUS',json.dumps(state,ensure_ascii=False));print('STUDY_TABLE',by_study.to_string());print('MARKERS',pd.read_csv(P/'legacy_object/MARKER_DESCRIPTIVE_MEANS.tsv',sep='\t').query("gene in ['SPP1','C1QC','GPNMB']").to_string(index=False))
