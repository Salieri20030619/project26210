from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase210_immune_identity_pairing_check_20260919_v1';M=R/'local_records_migration_20260908_150455/project';v1=M/'phase7_annotation_local_v1/runs/20260810_234538_376/tables';v2=M/'phase7_annotation_local_v2/runs/20260810_234538_376/tables'
a=pd.read_csv(v2/'final_cluster_annotations.tsv',sep='\t');sel=a[a.dataset.isin(['GSE264586','GSE294109','GSE207493'])&a.broad_annotation.isin(['T_cell','NK_cell'])];sel.to_csv(O/'EXTERNAL_TNK_CLUSTER_REVIEW.tsv',sep='\t',index=False)
rows=[]
for ds in ['GSE264586','GSE294109','GSE207493']:
 x=pd.read_csv(v1/f'{ds}_canonical_marker_expression.tsv',sep='\t');print(ds,'MARKER_COLS',list(x.columns));print(sel[sel.dataset==ds][['cluster','cluster_cells','fine_annotation','annotation_state','confidence','review_note']].to_string(index=False))
 x.to_csv(O/f'{ds}_EXISTING_CANONICAL_MARKERS.tsv',sep='\t',index=False)
for ds in ['GSE264586','GSE294109']:
 p=H/f'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2/{ds}_phase67_metadata_v1.tsv.gz';m=pd.read_csv(p,sep='\t')
 for k,g in m.groupby(['patient_key','sample_id','timepoint'],dropna=False):
  c=g.phase7_fine_annotation.value_counts();nmac=int(c.get('Macrophage_C1QC',0)+c.get('Macrophage_SPP1',0));cd8=int(c.get('CD8_cytotoxic',0));nk=int(c.get('NK_cell',0));cd4=int(c.get('CD4_naive_memory',0)+c.get('Treg',0));treg_available=ds=='GSE294109'
  rows.append(dict(dataset=ds,patient=k[0],sample_id=k[1],timepoint=k[2],n_mature_candidate=nmac,n_CD8_label=cd8,n_NK_label=nk,n_CD4_candidate=cd4,n_Treg_label=int(c.get('Treg',0)),CD8_count_ready=nmac>=50 and cd8>=20,NK_count_ready=nmac>=50 and nk>=20,Treg_CD4_count_ready=nmac>=50 and cd4>=50 and treg_available,Treg_annotation_available=treg_available,C4_available=False,qualification='Existing external annotation/count screen; no new expression validation or tumour provenance recuration'))
e=pd.DataFrame(rows);e.to_csv(O/'EXTERNAL_PATIENT_SAMPLE_READINESS.tsv',sep='\t',index=False)
for ds,g in e.groupby('dataset'):
 for ep in ['CD8','NK','Treg_CD4']:
  sub=g[g[ep+'_count_ready']];print('EXTERNAL_READY',ds,ep,'patients',sub.patient.nunique(),'sample_timepoints',len(sub))
# Independent consistency checks of saved tables and raw marker predicates.
d=pd.read_csv(O/'CELL_IDENTITY_AUDIT.tsv.gz',sep='\t',index_col=0);c=pd.read_csv(O/'CELL_PANEL_RAW_COUNTS.tsv.gz',sep='\t',index_col=0)
assert d.index.equals(c.index)
checks={'unique_TNK_IDs':d.index.is_unique and len(d)==86444,'T_backbone_detection':np.array_equal((c[['CD3D','CD3E','TRAC']]>0).sum(axis=1),d.T_backbone_detected),'Treg_marker_rule':np.array_equal(d.loc[d.legacy_candidate=='Treg','marker_supported'],(((c.CD3D>0).astype(int)+(c.CD3E>0).astype(int)+(c.TRAC>0).astype(int)>=2)&(c.FOXP3>0)&((c.IL2RA>0)|(c.CTLA4>0))).loc[d.legacy_candidate=='Treg']),'strict_requires_broad_compatible':not ((d.strict_identity_supported)&~d.broad_compatible).any(),'no_shared_NK_adaptors_in_myeloid_flag':'FCER1G' not in c.columns and 'TYROBP' not in c.columns}
t=pd.read_csv(O/'ALL_LIBRARY_IDENTITY_COUNTS.tsv',sep='\t');q=pd.read_csv(O/'SAME_TUMOR_LIBRARY_READINESS.tsv',sep='\t');checks['all75libraries']=len(t)==75;checks['C4_membership_total']=t.C4.sum()==2674;checks['mature_total']=t.mature_macrophage.sum()==22444;checks['macrophage_minimum']=not ((q.candidate_count_ready)&(q.n_mature<50)).any();checks['receiver_minimum']=np.array_equal(q.candidate_count_ready,(q.n_mature>=50)&(q.n_receiver_candidate>=q.receiver_threshold))
s=pd.read_csv(O/'READINESS_SUMMARY.tsv',sep='\t')
for _,row in s.iterrows():
 x=q[(q.endpoint==row.endpoint)&(q.receiver_threshold==row.threshold)];x=x[x.candidate_count_ready if row.identity_version=='candidate' else x.strict_count_ready];assert row.patients==x.patient.nunique() and row.libraries==len(x)
checks['patient_not_library_counts']=True
# Link actual gene coverage to candidate same-library readiness, without association testing.
g=pd.read_csv(O/'PROGRAM_EXPRESSION_COVERAGE.tsv',sep='\t');g['endpoint']=g.candidate.map({'CD8_candidate':'CD8','NK_candidate':'NK'});merged=g.merge(q[q.receiver_threshold==20],on=['study','patient','library','endpoint'],validate='many_to_one');summary=[]
for (ep,program),x in merged.groupby(['endpoint','program']):
 for label,mask in [('count_only',x.candidate_count_ready),('count_plus_coverage',x.candidate_count_ready&x.coverage_ready)]:
  z=x[mask];summary.append(dict(endpoint=ep,program=program,qualification=label,libraries=len(z),patients=z.patient.nunique(),studies=z.study.nunique()))
pd.DataFrame(summary).to_csv(O/'PROGRAM_READINESS_SUMMARY.tsv',sep='\t',index=False)
checks['passed']=all(checks.values());assert checks['passed'];(O/'VALIDATION.json').write_text(json.dumps(checks,indent=2));print('VALIDATION',checks);print(pd.DataFrame(summary).to_string(index=False))
print('TUMOR_TREG_MARKERS',d[(d.source_tissue=='Tumor')&(d.legacy_candidate=='Treg')].groupby('cell_type_broad_v1')[['marker_supported','strict_identity_supported']].agg(['sum','count']).to_string())
