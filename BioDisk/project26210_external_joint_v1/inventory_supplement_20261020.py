from pathlib import Path
import re,json,csv
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase337_full_supplement_inventory_20261020_v1';d.mkdir(exist_ok=True)
roots=[r,r/'local_records_migration_20260908_150455/project',r.parent]
entries=[]
for root in roots:
 if not root.exists():continue
 for p in root.iterdir():
  if not p.is_dir():continue
  m=re.search(r'phase[_ ]?(\d+)',p.name,re.I)
  if not m or (root==r.parent and '26210' not in p.name):continue
  n=int(m[1]);files=[]
  for sub in [p]+[p/x for x in ['delivery','reports','report','tables','analysis','qc','results'] if (p/x).is_dir()]:
   for f in sub.iterdir():
    if f.is_file() and f.suffix.lower() in ['.md','.txt','.tsv','.csv','.json','.pdf','.docx']:
     files.append({'path':str(f),'bytes':f.stat().st_size})
  entries.append({'phase':n,'root':str(p),'files':files})
entries.sort(key=lambda x:(x['phase'],x['root']))
(d/'PHASE_FILE_INVENTORY.json').write_text(json.dumps(entries,ensure_ascii=False,indent=2),encoding='utf-8')
with (d/'PHASE_DIRECTORY_INDEX.tsv').open('w',encoding='utf-8',newline='') as f:
 cw=csv.writer(f,delimiter='\t');cw.writerow(['phase','directory','indexed_files']);cw.writerows((e['phase'],e['root'],len(e['files'])) for e in entries)
print('directories',len(entries),'distinct phases',len(set(x['phase'] for x in entries)))
for e in entries:
 if e['phase']<100:print(e['phase'],Path(e['root']).name,','.join(Path(x['path']).name for x in e['files'] if x['path'].endswith('.md'))[:250])
