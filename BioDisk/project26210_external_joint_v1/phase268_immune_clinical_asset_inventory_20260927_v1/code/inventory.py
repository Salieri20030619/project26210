from pathlib import Path
import os,re,json,csv,hashlib,datetime,collections
import pandas as pd,numpy as np
P=Path(__file__).resolve().parents[1];R=P.parent;H=R.parent
root_patterns=['phase209_*','phase210_*','phase211_*','phase212_*','phase213_*','phase214_*','phase215_*','phase216_*','phase236_*','phase240_*','phase243_*','phase35_*','phase36_*','phase4e_*','phase4f_*','phase4g_*','phase4h_*','phase4i_*','phase144_*','phase168_*']
peer_patterns=['project26210_phase71_*','project26210_phase72_*','project26210_phase81_*','project26210_phase82_*','project26210_phase83_*','project26210_phase84_*','project26210_phase85_*','project26210_external_joint_v1_phase64_*','project26210_external_joint_v1_phase69_*','project26210_immunotherapy_scrna_phase60_*']
roots=[x for pattern in root_patterns for x in R.glob(pattern) if x.is_dir()]+[x for pattern in peer_patterns for x in H.glob(pattern) if x.is_dir()]
tcga=Path('/mnt/BioDisk/ludingyu/project26210/data_public/TCGA_KIRC/kirc_tcga_pan_can_atlas_2018');roots.append(tcga)
roots.extend(Path('/mnt/BioDisk/ludingyu/project26210/results')/n for n in ['part3','part3_3','part3_macro','part3_macro_subtype','part3_liana'])
excluded={'root_record_backups','root_record_backup','record_backup','final_root_record_backups','.git','__pycache__','deps','software','raw_extracted_v1','tf_regulon_v1','tf_regulon_de_novo_v1','raw_cel_v1.part','raw_cel_v2'}
files=[]
for root in roots:
 for base,dirs,names in os.walk(root):
  dirs[:]=[x for x in dirs if x not in excluded]
  for name in sorted(names):
   f=Path(base)/name
   try:
    st=f.stat();row=dict(root=str(root),path=str(f),relative_path=str(f.relative_to(root)),bytes=st.st_size,modified=datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone.utc).isoformat(),sha256='',verification='metadata_only')
    if st.st_size<=1024*1024:
     row['sha256']=hashlib.sha256(f.read_bytes()).hexdigest();row['verification']='current_sha256'
    files.append(row)
   except OSError as e:files.append(dict(root=str(root),path=str(f),relative_path=str(f.relative_to(root)),bytes='',modified='',sha256='',verification=str(e)))
main=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
if main.exists():files.append(dict(root=str(main.parent),path=str(main),relative_path=main.name,bytes=main.stat().st_size,modified=datetime.datetime.fromtimestamp(main.stat().st_mtime,datetime.timezone.utc).isoformat(),sha256='',verification='metadata_only_no_expression_reload'))
with (P/'tables/ASSET_FILE_INVENTORY.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(files[0]),delimiter='\t');w.writeheader();w.writerows(files)
# Recheck current frozen identities and archived analysis outputs against their existing manifests.
hashcheck=[]
for rootname,manifest in [('phase211_immune_receiver_refreeze_20260919_v1/freeze_v2','FREEZE_SHA256.json'),('phase213_immune_restart_execution_20260919_v1','OUTPUT_SHA256.json')]:
 root=R/rootname;m=json.loads((root/manifest).read_text());items=[dict(path=k,sha256=v) for k,v in m.items()] if isinstance(m,dict) else m
 for x in items:
  f=root/x['path'];actual=hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else None
  hashcheck.append(dict(manifest=str(root/manifest),path=str(f),expected=x['sha256'],actual=actual,match=actual==x['sha256']))
(P/'qc/FROZEN_MANIFEST_CHECKS.json').write_text(json.dumps(hashcheck,indent=2))
assert all(x['match'] for x in hashcheck)
base=H/'project26210_external_joint_v1_phase64_large_immune_validation_v1/inputs_v1'
cm=pd.read_csv(base/'CheckMate_Braun2020/processed_v1/checkmate_clinical_rna_subset_v1.tsv',sep='\t')
jv=pd.read_csv(base/'JAVELIN_Renal101/processed_v1/javelin_clinical_v1.tsv',sep='\t')
scores=pd.read_csv(H/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1/results_v1/phase81_sample_scores_v1.tsv',sep='\t')
print('SCORE_DATASETS',scores.dataset.value_counts().to_dict())
print('CM_DUPLICATES',int(cm.SUBJID.duplicated().sum()),int(cm.RNA_ID.duplicated().sum()))
print('JV_UNIQUE',len(jv),jv.ID.nunique())
obs=[]
for arm,d in cm.groupby('Arm'):
 row=dict(cohort='CheckMate',arm=arm,clinical_rows=len(d),unique_patient_ids=d.SUBJID.nunique(),unique_expression_ids=d.RNA_ID.nunique(),coding='event=PFS_CNSR / OS_CNSR')
 for end in ['PFS','OS']:
  t=pd.to_numeric(d[end],errors='coerce');ev=pd.to_numeric(d[end+'_CNSR'],errors='coerce');valid=(t>0)&np.isfinite(t)&ev.isin([0,1]);row[end+'_eligible']=int(valid.sum());row[end+'_events']=int(ev[valid].sum())
 row['MSKCC_missing']=int(d.MSKCC.isna().sum());row['age_missing']=int(d.Age.isna().sum());row['sex_missing']=int(d.Sex.isna().sum());row['response_counts']=d.ORR.fillna('MISSING').value_counts().to_dict();row['prior_treatment_counts']=d.Received_Prior_Therapy.fillna('MISSING').value_counts().to_dict();obs.append(row)
# Determine exact expression overlap using source sample IDs, without choosing outcome-based cutoffs.
for dataset,sc in scores.groupby('dataset'):
 overlap=set(sc.sample_id.astype(str))&set(jv.ID.astype(str))
 if overlap:
  d=jv[jv.ID.astype(str).isin(overlap)]
  for arm,z in d.groupby('TRT01P'):
   t=pd.to_numeric(z.PFS_P,errors='coerce');ev=1-pd.to_numeric(z.PFS_P_CNSR,errors='coerce');valid=(t>0)&np.isfinite(t)&ev.isin([0,1]);obs.append(dict(cohort='JAVELIN',arm=arm,clinical_rows=len(z),unique_patient_ids=z.ID.nunique(),PFS_eligible=int(valid.sum()),PFS_events=int(ev[valid].sum()),coding='event=1-PFS_P_CNSR',age_missing=int(z.AGE.isna().sum()),sex_missing=int(z.SEX.isna().sum()),available_columns=list(jv.columns)))
  assert len(d)==726
cohort=pd.read_csv(R/'phase36_part1_c05_tcga_joint_survival_v1/analysis_v1/c05_tcga_analysis_cohort_v1.tsv',sep='\t')
tc=dict(cohort='TCGA-KIRC C05',clinical_rows=len(cohort),unique_patient_ids=cohort.patient.nunique(),covariate_missing={x:int(cohort[x].isna().sum()) for x in ['age','sex','stage_simple']})
for end in ['OS','PFS']:
 t=pd.to_numeric(cohort[end+'_time'],errors='coerce');ev=pd.to_numeric(cohort[end+'_event'],errors='coerce');ok=(t>0)&np.isfinite(t)&ev.isin([0,1]);tc[end+'_eligible']=int(ok.sum());tc[end+'_events']=int(ev[ok].sum())
obs.append(tc)
# Gene-index presence in the stored full TCGA matrix (no expression modeling).
expr=tcga/'data_mrna_seq_v2_rsem.txt';genes=pd.read_csv(expr,sep='\t',usecols=[0]);targets=['TYMP','SPP1','MST1R','ITGAL','ITGB2','GPNMB','CCL18','CYP27A1','CD8A','CD8B','FOXP3','PDCD1'];presence={g:bool((genes.iloc[:,0]==g).any()) for g in targets}
(P/'qc/CLINICAL_ASSET_COUNTS.json').write_text(json.dumps(dict(cohorts=obs,clinical_total_JAVELIN=len(jv),TCGA_gene_index=presence),indent=2))
# Read actual stored model tables; keep raw-P nominal signals even when FDR does not pass.
model_inventory=[]
paths=[R/'phase213_immune_restart_execution_20260919_v1/I2_main_statistics/PRIMARY_12_HYPOTHESES.tsv',R/'phase243_donor_verified_immune_20260924_v1/tables/IMMUNE_PRIMARY.tsv',H/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1/results_v1/phase81_primary_models_v1.tsv',H/'project26210_phase82_paper_state_bridge_v1/results_v1/phase82_primary_models_v1.tsv',H/'project26210_phase83_bulk_state_projection_v1/results_v1/clinical_v1/phase83_primary_models_v1.tsv',R/'phase36_part1_c05_tcga_joint_survival_v1/analysis_v1/c05_joint_group_logrank_summary_v1.tsv']
for f in paths:
 d=pd.read_csv(f,sep='\t');model_inventory.append(dict(path=str(f),rows=len(d),columns=list(d.columns),records=d.where(pd.notna(d),None).to_dict(orient='records')))
(P/'tables/EXISTING_PRIMARY_MODEL_TABLES.json').write_text(json.dumps(model_inventory,indent=2,default=str))
summary=dict(status='INVENTORY_CHECKS_PASS',roots=len(roots),files=len(files),current_small_file_hashes=sum(bool(x['sha256']) for x in files),frozen_manifest_matches=len(hashcheck),clinical_counts=obs,TCGA_gene_index=presence,no_new_scores=True,no_new_models=True,no_new_figures=True,large_expression_matrices_not_loaded=True,excluded_directories=sorted(excluded))
(P/'qc/INVENTORY_CHECKS.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
