from pathlib import Path
import urllib.request,json,hashlib,datetime
P=Path(__file__).resolve().parents[1]
opener=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}))
urls={'article.html':'https://www.nature.com/articles/s41467-022-29366-6','github_tree.json':'https://api.github.com/repos/youqiongye/CRC_scRNAseq/git/trees/V1.0.0?recursive=1'}
rows=[]
for name,url in urls.items():
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
  with opener.open(req,timeout=45) as f: data=f.read(); final=f.url
  (P/'resources'/name).write_bytes(data)
  rows.append({'name':name,'url':url,'final_url':final,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'status':'OK'})
 except Exception as e: rows.append({'name':name,'url':url,'status':'ERROR','error':str(e)})
print(json.dumps(rows,indent=2),flush=True)
(P/'resources/INITIAL_FETCH.json').write_text(json.dumps(rows,indent=2))
