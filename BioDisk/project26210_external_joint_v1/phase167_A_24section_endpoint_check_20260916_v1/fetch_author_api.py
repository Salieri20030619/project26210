from pathlib import Path
import requests,json,base64,concurrent.futures
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase167_A_24section_endpoint_check_20260916_v1')
D=O/'metadata_api';D.mkdir(exist_ok=True)
t=json.loads((O/'metadata/author_repository_tree.json').read_text())
def fetch(x):
 u=x['url']
 try:
  r=requests.get(u,timeout=15);r.raise_for_status()
  b=base64.b64decode(r.json()['content']);(D/x['path']).write_bytes(b)
  return dict(url=u,path=x['path'],status=r.status_code,bytes=len(b))
 except Exception as e:return dict(url=u,path=x['path'],error=str(e))
items=[x for x in t['tree'] if x['type']=='blob' and x['path'].endswith(('.R','.sh','.md'))]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:logs=list(pool.map(fetch,items))
(D/'FETCH_LOG.json').write_text(json.dumps(logs,indent=2))
for p in D.glob('*.R'):
 for i,s in enumerate(p.read_text().splitlines()):
  if any(x in s.lower() for x in ['patient','clinical','ffpe','frozen','rename']): print(p.name,i+1,s[:220])
print(json.dumps(logs))
