from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,h5py
H=Path.home();R=H/'project26210_external_joint_v1';O=R/'phase211_immune_receiver_refreeze_20260919_v1';F=O/'freeze_v2';checks={}
d=pd.read_csv(F/'DISCOVERY_RECEIVER_CELL_FREEZE.tsv.gz',sep='\t',index_col='cell_id');parent=pd.read_csv(R/'phase210_immune_identity_pairing_check_20260919_v1/CELL_IDENTITY_AUDIT.tsv.gz',sep='\t',index_col='cell_id');checks['discovery_exact_full_ID_set']=d.index.is_unique and set(d.index)==set(parent.index) and len(d)==86444
checks['metadata_identity_unchanged']=d[['study','curated_patient_id','sample_id','source_tissue']].equals(parent.loc[d.index,['study','curated_patient_id','sample_id','source_tissue']])
c4path=R/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz';c4=pd.read_csv(c4path,sep='\t').cell_id;prior=json.loads((R/'phase210_immune_identity_pairing_check_20260919_v1/INPUT_MANIFEST.json').read_text());oldhash=next(x['sha256'] for x in prior if x['path']==str(c4path));checks['C4_unchanged_SHA256']=hashlib.sha256(c4path.read_bytes()).hexdigest()==oldhash;checks['C4_disjoint_receivers']=not set(c4)&set(d.index)
z=np.load(O/'RAW_PANEL_CACHE.npz');genes=list(z['genes']);raw=pd.DataFrame(z['counts'],index=z['cell_id'],columns=genes);checks['cache_exact_order']=np.array_equal(z['cell_id'],d.index)
# Independent marker summaries by frozen identity, separate from original cluster review.
markers=[]
for (tissue,label),g in d.groupby(['source_tissue','receiver_label_v1']):
 for gene in ['CD3D','CD3E','TRAC','CD4','IL7R','CD8A','CD8B','FOXP3','IKZF2','IL2RA','KLRD1','TRDC','TRGC1','TRGC2']:
  x=raw.loc[g.index,gene];markers.append(dict(tissue=tissue,label=label,gene=gene,cells=len(g),detected=int((x>0).sum()),detection_fraction=float((x>0).mean())))
pd.DataFrame(markers).to_csv(F/'FROZEN_IDENTITY_MARKER_SUMMARY.tsv',sep='\t',index=False)
treg=d[d.receiver_label_v1=='CD4_Treg'];checks['Treg_complete_cluster_not_detection_selected']=len(treg)==3956 and set(treg.index)==set(d.index[d.receiver_cluster==14]) and (raw.loc[treg.index,'FOXP3']==0).any()
# Fresh source gene-count identity spot check with stratified deterministic cell selection.
SRC=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad');probe=[]
for _,g in d.groupby('receiver_label_v1'):probe.extend(g.sample(n=min(3,len(g)),random_state=211).index)
with h5py.File(SRC,'r') as f:
 ids=pd.Index(f['obs/cell_id'].asstr()[:]);source_genes=f['var/gene'].asstr()[:];ptr=f['layers/counts/indptr'][:];x=f['layers/counts'];gc=['CD3D','CD8B','CD4','FOXP3','KLRD1','TRDC','PDCD1'];ix={gene:int(np.flatnonzero(source_genes==gene)[0]) for gene in gc}
 for cid in probe:
  i=ids.get_loc(cid);lo,hi=int(ptr[i]),int(ptr[i+1]);idx=x['indices'][lo:hi];val=x['data'][lo:hi];lookup=dict(zip(idx,val))
  for gene in gc:assert float(lookup.get(ix[gene],0))==raw.loc[cid,gene]
checks['fresh_source_marker_counts']=True
# Reconstruct counts from ID assignments independently and compare saved eligibility rows.
mac=pd.read_csv(R/'part2_rebuild_v1/08_regulons/input_v1/CELL_METADATA.tsv.gz',sep='\t',index_col=0,dtype={'candidate':str});mac=mac[~mac.candidate.isin(['1','Monocyte'])];e=pd.read_csv(F/'DISCOVERY_ENDPOINT_LIBRARY_FREEZE.tsv',sep='\t');pf=pd.read_csv(F/'DISCOVERY_PATIENT_ENDPOINT_FREEZE.tsv',sep='\t');definition=json.loads((F/'FREEZE_DEFINITION.json').read_text());programs=definition['programs'];countchecks=0
for row in e.itertuples():
 g=d[(d.study==row.study)&(d.curated_patient_id==row.patient)&(d.sample_id==row.library)&(d.source_tissue==row.tissue)];m=mac[(mac.study==row.study)&(mac.curated_patient_id==row.patient)&(mac.sample_id==row.library)&(mac.source_tissue==row.tissue)];assert len(m)==row.n_mature and (m.candidate=='4').sum()==row.n_C4
 labels=['CD4_Tcon','CD4_Treg'] if row.endpoint=='Treg_CD4' else (['CD8_T'] if row.endpoint.startswith('CD8') else ['NK']);receiver=g[g.receiver_label_v1.isin(labels)];assert len(receiver)==row.n_receiver
 u=~g.receiver_label_v1.isin(['CD4_Tcon','CD4_Treg','CD8_T','NK','gdT_enriched']);assert abs(u.mean()-row.unresolved_fraction)<1e-12
 cover=True
 if row.endpoint in programs:
  detection=(raw.loc[receiver.index,programs[row.endpoint]]>0).sum();n=int((detection>=5).sum());assert n==row.covered_genes;cover=n>=3
 valid=row.tissue=='Tumor' and len(m)>=50 and len(receiver)>=(50 if row.endpoint=='Treg_CD4' else 20) and u.mean()<=0.2 and cover
 assert bool(row.eligible)==valid;assert bool(row.C4_intrinsic_eligible)==(valid and (m.candidate=='4').sum()>=20);countchecks+=1
checks['library_counts_and_eligibility_recomputed']=True
for row in pf.itertuples():
 g=e[(e.study==row.study)&(e.patient==row.patient)&(e.endpoint==row.endpoint)&e.eligible];s=set(str(row.eligible_libraries).split(';')) if row.eligible else set();assert s==set(g.library) and row.n_libraries==len(g)
checks['patient_library_lists_and_no_pseudoreplication']=not pf.duplicated(['study','patient','endpoint']).any()
# External frozen labels must have exact source IDs and preserve patient/treatment/tissue fields.
ee=pd.read_csv(F/'EXTERNAL_ENDPOINT_SAMPLE_FREEZE.tsv',sep='\t');ext_checks=[]
for ds in ['GSE264586','GSE294109']:
 a=pd.read_csv(F/f'{ds}_RECEIVER_CELL_FREEZE.tsv.gz',sep='\t',index_col='cell_id');src=pd.read_csv(O/f'external/{ds}_CELL_MAPPING.tsv.gz',sep='\t',index_col='cell_id');counts=pd.read_csv(O/f'external/{ds}_TNK_PANEL_RAW_COUNTS.tsv.gz',sep='\t',index_col='cell_id');assert a.index.is_unique and set(a.index)==set(counts.index)
 assert a[['patient_key','sample_id','timepoint','treatment','rds_phase7_cluster_res0_8']].equals(src.loc[a.index,['patient_key','sample_id','timepoint','treatment','rds_phase7_cluster_res0_8']]);assert not (a.receiver_label_v1=='NK').any()
 for row in ee[ee.dataset==ds].itertuples():
  u=a[(a.patient_key==row.patient)&(a.sample_id==row.sample_id)&(a.timepoint==row.timepoint)&(a.treatment==row.treatment)];base=src[(src.patient_key==row.patient)&(src.sample_id==row.sample_id)&(src.timepoint==row.timepoint)];nm=int(base.phase7_fine_annotation.isin(['Macrophage_C1QC','Macrophage_SPP1']).sum());assert nm==row.n_mature
  labs=['CD4_Tcon','CD4_Treg'] if row.endpoint=='Treg_CD4' else (['CD8_T'] if row.endpoint.startswith('CD8') else ['NK']);rec=u[u.receiver_label_v1.isin(labs)];assert len(rec)==row.n_receiver
  frac=(~u.receiver_label_v1.isin(['CD8_T','CD4_Tcon','CD4_Treg','NK'])).mean();assert abs(frac-row.unresolved_fraction)<1e-12
  tissue_ok=set(u.rds_source_tissue_status)=={'Tumor'};cov=True
  if row.endpoint in programs:
   ng=int(((counts.loc[rec.index,programs[row.endpoint]]>0).sum()>=5).sum());assert ng==row.covered_genes;cov=ng>=3
  valid=tissue_ok and nm>=50 and len(rec)>=(50 if row.endpoint=='Treg_CD4' else 20) and frac<=0.2 and cov and not (ds=='GSE264586' and row.endpoint=='Treg_CD4');assert bool(row.eligible)==valid
 ext_checks.append({'dataset':ds,'cells':len(a),'exact_ID_mapping':True,'eligibility_recomputed':True})
checks['external_IDs_and_eligibility_recomputed']=True;checks['external_NK_conflict_not_mislabeled']=not ee.loc[ee.endpoint=='NK_cytotoxicity','eligible'].any();checks['external_C4_not_substituted']=not ee.C4_available.any();checks['no_association_tests']=definition['association_tests_run']==False
checks={k:bool(v) for k,v in checks.items()};assert all(checks.values());out={'status':'PASS','checks':checks,'n_checks':len(checks),'discovery_library_endpoint_rows_recomputed':countchecks,'fresh_source_cells_checked':len(probe),'external':ext_checks};(O/'VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
stat=json.loads((F/'STATUS.json').read_text());stat['status']='FROZEN_VALIDATED';stat['validation']=str(O/'VALIDATION.json');(F/'STATUS.json').write_text(json.dumps(stat,indent=2))
print('FREEZE_VALIDATION_PASS')
