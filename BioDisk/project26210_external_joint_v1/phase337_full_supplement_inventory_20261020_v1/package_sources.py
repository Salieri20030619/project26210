from pathlib import Path
import json,re,shutil,hashlib,os,tarfile,csv
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');inv=r/'phase337_full_supplement_inventory_20261020_v1';d=r/'phase338_expanded_supplement_20261020_v2'
for s in ['statistics','qc','code','delivery']:(d/s).mkdir(parents=True,exist_ok=True)
mods=json.loads((inv/'MODULES.json').read_text(encoding='utf-8'));catalog=json.loads((inv/'DEEP_RESULT_CATALOG.json').read_text(encoding='utf-8'))
manifest=[];skipped=[]
for mod in mods:
 if mod['id']=='E30':continue
 folder=d/'statistics'/mod['id'];folder.mkdir(exist_ok=True);seen=set()
 for n in mod['phases']:
  options=[x for x in catalog if x['phase']==n]
  if not options:skipped.append([mod['id'],n,'','未在已收集的详细目录找到统计表']);continue
  # Current root first, then actual server analysis, then migrated snapshot.
  options.sort(key=lambda x:(0 if Path(x['root']).parent==r else 1 if Path(x['root']).parent==r.parent else 2,-len(x['files'])))
  entry=options[0];root=Path(entry['root'])
  if n in [179,180] and (root/'run_v2').exists():walk=root/'run_v2'
  elif n==200 and (root/'delivery/统计表').exists():walk=root/'delivery/统计表'
  elif n==213:walk=root/'delivery_v1'
  elif n in [103,104]:walk=root/'results_v1'
  elif (root/'tables').exists():walk=root/'tables'
  else:walk=root
  for here,dirs,names in os.walk(walk):
   depth=len(Path(here).relative_to(walk).parts)
   dirs[:]=[x for x in dirs if not any(t in x.lower() for t in ['backup','failed','attempt','before','archive','input','raw','cache','script','code','figure','qc','validation','source','resource'])] if depth<3 else []
   if n==74:dirs[:]=[x for x in dirs if 'nominal_target_sensitivity_v1'!=x]
   for name in names:
    p=Path(here)/name;rel=str(p.relative_to(root))
    if not (name.endswith('.tsv') or name.endswith('.tsv.gz')):continue
    if any(x in name.lower() for x in ['manifest','hash','validation','preflight','command','software','run_metadata','status','gene_input','model_membership','interaction_membership','analysis_input','cell_source','cell_metadata','cell_labels']):continue
    if n==292 and name=='ASSOCIATIONS.tsv':continue
    if n==291 and name=='SPECTRAL_SECTION.tsv':continue
    if p.stat().st_size>10000000:skipped.append([mod['id'],n,str(p),'大型逐项辅助表保留服务器']);continue
    data=p.read_bytes();sha=hashlib.sha256(data).hexdigest()
    if sha in seen:continue
    seen.add(sha);clean=re.sub(r'[^A-Za-z0-9_.-]+','_',name)
    target=f'P{n}_{len(manifest)+1:04}_{clean}'
    (folder/target).write_bytes(data)
    manifest.append({'module':mod['id'],'phase':n,'file':str(Path(mod['id'])/target),'source':str(p),'bytes':len(data),'sha256':sha})
(d/'qc/STATISTICAL_SOURCE_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
(d/'qc/TABLE_EXCLUSIONS.json').write_text(json.dumps(skipped,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(inv/'MODULES.json',d/'code/MODULES.json')
with tarfile.open(d/'statistics_bundle.tar.gz','w:gz') as t:
 t.add(d/'statistics',arcname='statistics');t.add(d/'qc/STATISTICAL_SOURCE_MANIFEST.json',arcname='qc/STATISTICAL_SOURCE_MANIFEST.json');t.add(d/'qc/TABLE_EXCLUSIONS.json',arcname='qc/TABLE_EXCLUSIONS.json')
print('Tables',len(manifest),'bytes',sum(x['bytes'] for x in manifest),'bundle',(d/'statistics_bundle.tar.gz').stat().st_size)
print('Modules',[(m['id'],sum(x['module']==m['id'] for x in manifest)) for m in mods]);print('Excluded',skipped)
