from pathlib import Path
import json,re,hashlib,tarfile,os
r=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1');d=r/'phase337_full_supplement_inventory_20261020_v1';entries=json.loads((d/'PHASE_FILE_INVENTORY.json').read_text());dest=d/'evidence_reports';dest.mkdir(exist_ok=True)
selected={25,28,35,36,57,58,60,62,63,67,69,71,72,73,74,76,77,79,81,83,84,85,87,88,98,99,103,104,109,114,126,133,150,151,154,157,158,161,163,174,175,178,179,180,182,183,184,186,187,200,213,214,220,222,234,242,245,248,257,269,270,271,275,276,277,278,280,282,284,287,289,290,291,292,293,295,303,305,306,307,308,309,311,312,313,314,317,318,319,320,321}
catalog=[];seen=set()
for e in entries:
 if e['phase'] not in selected:continue
 p=Path(e['root']);files=[]
 for root,dirs,names in os.walk(p):
  depth=len(Path(root).relative_to(p).parts)
  dirs[:]=[x for x in dirs if not any(k in x.lower() for k in ['backup','failed','attempt','raw','cache','input','vendor','env','figure','script','code','image'])] if depth<3 else []
  for n in names:
   f=Path(root)/n
   if f.suffix.lower() not in ['.md','.txt','.tsv','.json'] or f.stat().st_size>5000000:continue
   if f.suffix.lower() in ['.md','.txt'] and f.stat().st_size<150000 and not any(k in n.lower() for k in ['plan','prereg','registry','protocol','readme','failure','before','session','software','manifest']):
    data=f.read_bytes();h=hashlib.sha256(data).hexdigest()
    if h not in seen:
     seen.add(h);target=dest/f'{e["phase"]:03}_{len(catalog):04}_{n}';target.write_bytes(data)
     files.append({'source':str(f),'local':target.name,'sha256':h,'bytes':len(data)})
   elif f.suffix.lower()=='.tsv':files.append({'source':str(f),'bytes':f.stat().st_size})
 catalog.append({'phase':e['phase'],'root':str(p),'files':files})
(d/'DEEP_RESULT_CATALOG.json').write_text(json.dumps(catalog,ensure_ascii=False,indent=2),encoding='utf-8')
with tarfile.open(d/'reports_bundle.tar.gz','w:gz') as t:t.add(dest,arcname='evidence_reports');t.add(d/'DEEP_RESULT_CATALOG.json',arcname='DEEP_RESULT_CATALOG.json')
print('Reports',len(seen),'catalog directories',len(catalog),'bytes',(d/'reports_bundle.tar.gz').stat().st_size)
