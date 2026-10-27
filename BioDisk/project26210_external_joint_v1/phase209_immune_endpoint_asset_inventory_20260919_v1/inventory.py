from pathlib import Path
import csv,gzip,json,hashlib,collections,datetime,h5py
h=Path.home(); root=h/'project26210_external_joint_v1'; out=root/'phase209_immune_endpoint_asset_inventory_20260919_v1'
def write(name,rows):
 if not rows:return
 with (out/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def read(p):
 with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
p=Path('/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/03_malignancy_cnv/cnv_annotated_h5ad_v1/part1_scvi_integrated_cell_annotated_cnv_v1.h5ad')
with h5py.File(p) as f:
 def col(k):
  a=f['obs'][k]
  if isinstance(a,h5py.Group):
   cats=a['categories'].asstr()[:] if a['categories'].dtype.kind in 'OSU' else a['categories'][:]
   return [str(cats[i]) if i>=0 else '' for i in a['codes'][:]]
  return a.asstr()[:].tolist() if a.dtype.kind in 'OSU' else a[:].tolist()
 cols={k:col(k) for k in ['cell_id','curated_patient_id','sample_id','study','source_tissue','cell_type_broad_v1','annotation_confidence_v1','annotation_ambiguous_v1']}
 print('FIELDS_LOADED',len(cols['cell_id']))
 ids=read(root/'part2_rebuild_v1/09_downstream_crosswalk/C4_CANDIDATE_CELL_IDS.tsv.gz');print('C4_FIELDS',list(ids[0]))
 c4=set(x.get('cell_id','') for x in ids);assert len(c4)==2674
 counter=collections.Counter();lineages=collections.defaultdict(lambda: {'n':0,'patients':set(),'libraries':set()})
 for i,cid in enumerate(cols['cell_id']):
  study,patient,lib,tissue,label=[cols[k][i] for k in ['study','curated_patient_id','sample_id','source_tissue','cell_type_broad_v1']]
  for lab in ([label,'C4_frozen'] if cid in c4 else [label]):
   counter[(study,patient,lib,tissue,lab)]+=1
   for scope in ['All',tissue]:
    d=lineages[(scope,lab)];d['n']+=1;d['patients'].add(patient);d['libraries'].add((study,lib))
 write('PATIENT_LIBRARY_LINEAGE_COUNTS.tsv',[dict(study=k[0],patient=k[1],library=k[2],tissue=k[3],lineage=k[4],n_cells=v) for k,v in sorted(counter.items())])
 summary=[dict(scope=k[0],lineage=k[1],cells=v['n'],patients=len(v['patients']),libraries=len(v['libraries'])) for k,v in sorted(lineages.items())]
 write('DISCOVERY_LINEAGE_COVERAGE.tsv',summary)
 match=collections.defaultdict(collections.Counter)
 for (s,pt,lib,tis,lab),n in counter.items():
  if tis=='Tumor':match[(s,pt,lib)][lab]+=n
 pairs=[]
 for labs in [('Macrophage','T cell'),('Macrophage','NK cell'),('C4_frozen','T cell'),('C4_frozen','NK cell')]:
  keys=[k for k,v in match.items() if all(v[x]>0 for x in labs)]
  pairs.append(dict(predictor=labs[0],receiver=labs[1],tumor_libraries_with_presence=len(keys),patients_with_presence=len(set((k[0],k[1]) for k in keys)),qualification='Presence only; no expression or minimum-cell eligibility established'))
 write('DISCOVERY_TUMOR_CO_PRESENCE.tsv',pairs)
 print('SUMMARY',json.dumps(summary));print('PAIRING',json.dumps(pairs))
 print('VAR_FIELDS',list(f['var']))
source_dirs=[root/'part2_rebuild_v1/09_downstream_crosswalk',h/'project26210_phase71_immune_advanced_v1',h/'project26210_phase72_immune_ecosystem_v1',h/'project26210_external_joint_v1_phase64_large_immune_validation_v1/inputs_v1',h/'project26210_external_joint_v1_phase66_immune_program_transfer_v1',h/'project26210_external_joint_v1_phase67_direct_program_replication_v1',h/'project26210_external_joint_v1_phase69_r2_composite_clinical_validation_v1',h/'project26210_phase81_bulk_ici_transcriptome_rebuild_v1',root/'phase109_targeted_gap_analysis_local_v1',root/'phase35_part1_c01_c04_cell_communication_v1']
rows=[]
for d in source_dirs:
 for x in sorted(d.rglob('*')):
  if x.is_file():rows.append(dict(path=str(x),bytes=x.stat().st_size,mtime=datetime.datetime.fromtimestamp(x.stat().st_mtime).isoformat(),sha256=hashlib.sha256(x.read_bytes()).hexdigest() if x.stat().st_size<2000000 else 'not_rehashed_large_file'))
rows.append(dict(path=str(p),bytes=p.stat().st_size,mtime=datetime.datetime.fromtimestamp(p.stat().st_mtime).isoformat(),sha256='historical_e4779c46af13a8643e8dec584eeff4322258181f5edab3e55dcd8c20bc28a9e7_not_rehashed'))
write('ASSET_FILE_INVENTORY.tsv',rows)
for cohort in ['GSE264586','GSE294109']:
 meta=h/f'project26210_external_joint_v1_phase67_direct_program_replication_v1/inputs_v1/metadata_v2/{cohort}_phase67_metadata_v1.tsv.gz'
 data=read(meta);print(cohort,'FIELDS',list(data[0]))
 for key in ['phase7_broad_annotation','phase7_fine_annotation','timepoint','treatment','response']:
  cnt=collections.Counter(x.get(key,'') for x in data);print(cohort,key,dict(cnt))
 write(cohort+'_ANNOTATION_COUNTS.tsv',[dict(field=key,label=k,cells=n) for key in ['phase7_broad_annotation','phase7_fine_annotation','timepoint','treatment','response'] for k,n in collections.Counter(x.get(key,'') for x in data).items()])
(out/'INVENTORY_STATUS.json').write_text(json.dumps({'status':'COMPLETE_ASSET_METADATA_INVENTORY','source_cells':len(cols['cell_id']),'C4_ids':len(c4),'C4_ids_matched':sum(x in c4 for x in cols['cell_id']),'files':len(rows),'no_new_models':True,'scope':'Metadata and existing files only; no new cell identity or endpoint expression validation'},indent=2))
print('INVENTORY_COMPLETE',len(rows))
