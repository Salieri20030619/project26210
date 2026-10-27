from pathlib import Path
import json,hashlib
base=Path('/mnt/BioDisk/wangmingyuan/.codex/sessions/2026/08/05');o=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1/reference_transition_trace');o.mkdir()
needles=['build_patient_balanced_reference_v2','within-patient continuous program','SPP1_PROGRAM','reference_v2_cell_assignments','SPP1/C1QC','巨噬细胞']
def strings(v):
 if isinstance(v,str):yield v
 elif isinstance(v,dict):
  for k,x in v.items():
   if k not in ['encrypted_content']:yield from strings(x)
 elif isinstance(v,list):
  for x in v:yield from strings(x)
hits=[];scanned=[]
for p in sorted(base.rglob('*.jsonl')):
 try:
  with p.open() as h:first=json.loads(h.readline())
 except Exception:continue
 if first.get('payload',{}).get('cwd','') not in ['/mnt/BioDisk/wangmingyuan','/mnt/BioDisk/ludingyu/project26210']:continue
 scanned.append({'path':str(p),'bytes':p.stat().st_size})
 with p.open() as h:
  for i,l in enumerate(h,1):
   if not any(x in l for x in needles):continue
   try:d=json.loads(l)
   except Exception:continue
   texts=[s for s in strings(d.get('payload',{})) if any(x in s for x in needles)]
   if not texts:continue
   k=len(hits)+1;f=o/f'event_{k:03d}.json';record={'source':str(p),'source_line':i,'timestamp':d.get('timestamp'),'type':d.get('type'),'payload_type':d.get('payload',{}).get('type'),'role':d.get('payload',{}).get('role'),'event_sha256':hashlib.sha256(l.encode()).hexdigest(),'matching_texts':texts}
   f.write_text(json.dumps(record,ensure_ascii=False,indent=2));ix={k:v for k,v in record.items() if k!='matching_texts'};ix['event']=f.name;ix['excerpt']=texts[0][:250];hits.append(ix)
(o/'MATCH_INDEX.json').write_text(json.dumps(hits,ensure_ascii=False,indent=2));(o/'SCANNED_PROJECT_SESSIONS.json').write_text(json.dumps(scanned,ensure_ascii=False,indent=2))
print('PROJECT_SESSIONS',len(scanned),'MATCHES',len(hits))
for r in hits:print(r['event'],r['timestamp'],r['payload_type'],r['role'],r['excerpt'][:180].replace('\n',' '))
