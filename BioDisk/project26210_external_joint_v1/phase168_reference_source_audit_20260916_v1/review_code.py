from pathlib import Path
import json,requests,base64,concurrent.futures,re
O=Path('/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase168_reference_source_audit_20260916_v1')
tree=json.loads((O/'author_repo_tree.json').read_text());D=O/'author_api';D.mkdir(exist_ok=True)
items=[x for x in tree['tree'] if x['type']=='blob' and x.get('size',0)<1200000 and x['path'].endswith(('.R','.md','.csv','.tsv','.txt','.ipynb'))]
def get(x):
 try:
  r=requests.get(x['url'],timeout=(10,20));r.raise_for_status();b=base64.b64decode(r.json()['content'])
  p=D/x['path'];p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
  return dict(path=x['path'],bytes=len(b),status=r.status_code)
 except Exception as e:return dict(path=x['path'],error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:logs=list(pool.map(get,items))
(O/'AUTHOR_API_FETCH_LOG.json').write_text(json.dumps(logs,indent=2))
hits=[]
for f in D.rglob('*'):
 if not f.is_file():continue
 if f.suffix in ['.R','.md','.txt','.ipynb']:
  for i,l in enumerate(f.read_text(errors='replace').splitlines()):
   if re.search(r'GSM|GSE|ffpe|frozen|patient|c_2|a_17|Load10X|readRDS|read.csv|load\(',l,re.I):
    hits.append(dict(file=str(f.relative_to(D)),line=i+1,text=l[:500]))
 elif f.suffix in ['.csv','.tsv']:
  print('TABLE',f.relative_to(D),'\n','\n'.join(f.read_text(errors='replace').splitlines()[:2])[:800])
(O/'CODE_SOURCE_HITS.json').write_text(json.dumps(hits,ensure_ascii=False,indent=2))
for h in hits:
 if re.search(r'GSM|GSE|ffpe|frozen|patient|Load10X|readRDS|load\(',h['text'],re.I):print(h['file'],h['line'],h['text'])
print('FETCH',len(logs),sum('error' not in x for x in logs))
