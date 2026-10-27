from pathlib import Path
import requests,json,hashlib,re,datetime
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase167_A_24section_endpoint_check_20260916_v1')
D=O/'metadata';D.mkdir(exist_ok=True)
log=[]
def get(u,path):
 try:
  r=requests.get(u,timeout=25)
  log.append(dict(url=u,status=r.status_code,bytes=len(r.content),retrieved=datetime.datetime.now().astimezone().isoformat()))
  r.raise_for_status();path.write_bytes(r.content)
  return r
 except Exception as e:
  log[-1:]=log[-1:] if log else []
  log.append(dict(url=u,error=str(e)))
  return None
r=get('https://api.github.com/repos/maximemeylan/Meylan_et_al_2022/git/trees/main?recursive=1',D/'author_repository_tree.json')
if r:
 tree=r.json();commit=tree['sha']
 for x in tree['tree']:
  if x['type']=='blob' and x['size']<100000 and x['path'].endswith(('.R','.sh','.md')):
   get('https://raw.githubusercontent.com/maximemeylan/Meylan_et_al_2022/'+commit+'/'+x['path'],D/x['path'])
get('https://api.github.com/repos/maximemeylan/Meylan_et_al_2022/issues?state=all',D/'author_repository_issues.json')
(O/'SOURCE_FETCH_LOG.json').write_text(json.dumps(log,indent=2))
for f in sorted(D.glob('*')):
 if f.suffix in ['.R','.sh','.md']:
  lines=f.read_text().splitlines()
  for i,line in enumerate(lines):
   if re.search(r'patient|clinical|metadata|ffpe|frozen|c_2|[Ss]ample',line):
    print(f.name,i+1,line[:250])
