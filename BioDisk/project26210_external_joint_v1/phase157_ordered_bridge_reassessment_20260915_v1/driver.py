from pathlib import Path
import sys,json,subprocess
root=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase157_ordered_bridge_reassessment_20260915_v1')
for line in sys.stdin:
 try:
  q=json.loads(line)
  if q.get('action')=='close':break
  p=(root/q['file']).resolve();assert root in p.parents
  p.parent.mkdir(parents=True,exist_ok=True)
  p.write_text(q['code'],encoding='utf-8')
  if q.get('run',True):
   r=subprocess.run([q.get('python',sys.executable),str(p)],capture_output=True,text=True)
   result=dict(id=q['id'],code=r.returncode,stdout=r.stdout,stderr=r.stderr)
  else:result=dict(id=q['id'],code=0,stdout='SAVED',stderr='')
 except Exception as e:result=dict(id=q.get('id','error'),code=1,stdout='',stderr=repr(e))
 payload=json.dumps(result,ensure_ascii=True).encode().hex()
 print('@@START_'+result['id'],flush=True)
 for i in range(0,len(payload),60):print(':'+payload[i:i+60],flush=True)
 print('@@END_'+result['id'],flush=True)
