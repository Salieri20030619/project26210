from pathlib import Path
import json,hashlib
p=Path('/mnt/BioDisk/ludingyu/.codex/sessions/2026/07/17/rollout-2026-07-17T09-45-41-019f6dc0-35f6-7f61-8004-d55e18fd39d1.jsonl')
o=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase145_macrophage_provenance_trace_remote_v1/session_trace');o.mkdir()
needles=['02_macrophage_final_subtypes','macro_cluster_final','pca_macro_clean','macro_subtype','SPP1+Macro']
def strings(v):
 if isinstance(v,str):yield v
 elif isinstance(v,dict):
  for k,x in v.items():
   if k not in ['encrypted_content']:yield from strings(x)
 elif isinstance(v,list):
  for x in v:yield from strings(x)
hits=[]
with p.open() as h:
 for i,l in enumerate(h,1):
  if not any(x in l for x in needles):continue
  try:d=json.loads(l)
  except Exception:continue
  texts=[s for s in strings(d.get('payload',{})) if any(x in s for x in needles)]
  if not texts:continue
  k=len(hits)+1;f=o/f'event_{k:03d}.json';f.write_text(json.dumps(d,ensure_ascii=False,indent=2))
  excerpts=[]
  for t in texts:
   for n in needles:
    j=t.find(n)
    if j>=0:excerpts.append(t[max(0,j-200):min(len(t),j+750)])
  hits.append({'source_line':i,'event':f.name,'timestamp':d.get('timestamp'),'type':d.get('type'),'payload_type':d.get('payload',{}).get('type'),'role':d.get('payload',{}).get('role'),'excerpts':excerpts[:5]})
(o/'MATCH_INDEX.json').write_text(json.dumps(hits,ensure_ascii=False,indent=2))
print('MATCHES',len(hits))
for r in hits:print(json.dumps({k:v for k,v in r.items() if k!='excerpts'},ensure_ascii=False),r['excerpts'][0][:650])
