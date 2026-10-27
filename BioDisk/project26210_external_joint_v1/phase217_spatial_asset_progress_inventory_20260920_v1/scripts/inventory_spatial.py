from pathlib import Path
import csv,json,gzip,hashlib,datetime,collections
P=Path(__file__).resolve().parents[1];R=P.parent;BASE=R.parent
read=lambda f:list(csv.DictReader(f.open(),delimiter='\t'))
q=read(R/'phase165_A_input_confirmation_20260915_v1/SECTION_INPUT_QC.tsv')
patient=read(R/'phase168_reference_source_audit_20260916_v1/SPATIAL_26_SECTION_PATIENT_MAP.tsv')
(P/'source_snapshots/PATIENT_MAP_SOURCE.json').write_text(json.dumps(patient,indent=2,ensure_ascii=False))
s54=BASE/'project26210_external_joint_v1_phase54_gse175540_spatial_v1'
s78=BASE/'project26210_phase78_fine_state_spatial_v1'
s55=BASE/'project26210_external_joint_v1_phase55_cell2location_v5_spawn_resume_20260821_v1'
assets=[];sections=[]
def asset(path,kind,dataset):
 f=Path(path);v=dict(path=str(f),kind=kind,dataset=dataset,exists=f.is_file())
 if f.is_file():
  st=f.stat();v.update(bytes=st.st_size,mtime=st.st_mtime)
  if st.st_size<2_000_000 and f.suffix in ['.md','.json','.tsv','.py']:v['sha256']=hashlib.sha256(f.read_bytes()).hexdigest()
 assets.append(v);return v['exists']
for x in q:
 sid=x['section'];d=s54/'run_v1/sections_v1'/sid;e=s78/'run_v1/sections_v1'/sid
 ss=json.loads((d/'section_summary_v1.json').read_text());row={'section':sid,'preparation':x['preparation'],'spots':int(x['spots']),'genes':int(x['genes']),'historical_qc_status':ss['status']}
 for name,path in [('processed',d/'processed_section_v1.h5ad'),('abundance',e/'fine_state_abundance_v1.tsv.gz'),('posterior',e/'fine_state_posterior_v1.h5ad'),('section_pass',d/'SECTION_PASS_v1'),('mapping_pass',e/'MAPPING_PASS_v1')]:row[name]=asset(path,name,'GSE175540')
 with gzip.open(e/'fine_state_abundance_v1.tsv.gz','rt') as f:
  header=f.readline().strip().split('\t');row['abundance_rows']=sum(1 for _ in f)
 row['abundance_rows_match']=row['abundance_rows']==row['spots']
 row['abundance_columns']=header
 for label,pattern in [('raw_h5','*raw_feature_bc_matrix.h5'),('filtered_h5','*filtered_feature_bc_matrix.h5'),('coordinates','*tissue_positions_list.csv.gz'),('scalefactors','*scalefactors_json.json.gz'),('hires_image','*tissue_hires_image.png.gz'),('lowres_image','*tissue_lowres_image.png.gz'),('TLS_annotation','*TLS_annotation.csv.gz')]:
  paths=list((s54/'input_extract_v1').glob(sid+pattern));row[label]=len(paths)
  for f in paths:asset(f,label,'GSE175540')
 sections.append(row)
full=read(R/'phase37_part4_s02_s04_spatial_metabolism_lr_v1/full_expression_manifest_v1/full_expression_input_manifest_v1.tsv')
legacy=[]
for row in full:
 row=dict(row);row['present_now']=asset(row['input_h5ad'],'legacy_full_expression',row['cohort']);legacy.append(row)
roots=[s54,s78,s55,BASE/'project26210_external_joint_v1_phase70_cell2location_spatial_composition_v1',BASE/'project26210_phase75_mta_spatial_transfer_v1',BASE/'project26210_part4_composition_niche_phase18_20260813_v1',BASE/'project26210_tymp_spp1_spatial_phase19_20260814_v2',BASE/'project26210_spp1_vascular_phase20_20260814_v1',R/'phase37_part4_s02_s04_spatial_metabolism_lr_v1',R/'phase12_glmf_spatial_analysis_v1']
for num in [164,165,166,167,168,169,170,172,188,110]:roots+=list(R.glob('phase'+str(num)+'_*'))
for d in roots:
 for f in d.rglob('*'):
  if not f.is_file() or any(x in f.parts for x in ['authority_backups','root_record_backups','__pycache__','source_snapshots']):continue
  if f.suffix in ['.md','.json','.tsv'] and f.stat().st_size<2_000_000:asset(f,'record_or_result',d.name)
unique={a['path']:a for a in assets};assets=list(unique.values())
summary={'created_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'READ_ONLY_ASSET_PROGRESS_AUDIT_NO_NEW_ANALYSIS','sections_gse175540':len(sections),'spots_gse175540':sum(x['spots'] for x in sections),'preparations':dict(collections.Counter(x['preparation'] for x in sections)),'coverage':{k:sum(bool(x[k]) for x in sections) for k in ['processed','abundance','posterior','section_pass','mapping_pass','abundance_rows_match','raw_h5','filtered_h5','coordinates','scalefactors','hires_image','lowres_image','TLS_annotation']},'tls_annotation_absent_sections':[x['section'] for x in sections if not x['TLS_annotation']],'legacy_input_count':len(legacy),'legacy_input_existing':sum(x['present_now'] for x in legacy),'legacy_spots_by_cohort':{c:sum(int(x['n_spots']) for x in legacy if x['cohort']==c) for c in {x['cohort'] for x in legacy}},'inventoried_files':len(assets),'missing_files':[a for a in assets if not a['exists']],'validation_level':'Current file presence/size; current abundance row counts; source JSON and small-table content. Prior H5AD and numerical validations inherited, not rerun. No histology registration QA this audit.'}
for name,obj in [('ASSET_INVENTORY.json',assets),('GSE175540_SECTION_ASSETS.json',sections),('LEGACY14_EXPRESSION_ASSETS.json',legacy),('INVENTORY_SUMMARY.json',summary)]: (P/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False))
assert len(sections)==24 and sum(x['spots'] for x in sections)==72952
assert all(x['processed'] and x['abundance'] and x['posterior'] and x['abundance_rows_match'] for x in sections)
print(json.dumps(summary,ensure_ascii=False))
