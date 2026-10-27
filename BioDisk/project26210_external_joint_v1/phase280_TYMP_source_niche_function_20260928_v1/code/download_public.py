from pathlib import Path
import urllib.request,json,hashlib,datetime
p=Path(__file__).resolve().parents[1]/'resources/GSE171542'
urls={'GSE171542_RAW.tar':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171542/suppl/GSE171542_RAW.tar','GSE171542_family.soft.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE171nnn/GSE171542/soft/GSE171542_family.soft.gz'}
rows=[]
for name,url in urls.items():
 print('DOWNLOAD',name,flush=True)
 with urllib.request.urlopen(url,timeout=60) as r, (p/name).open('wb') as f:
  while True:
   b=r.read(1024*1024)
   if not b:break
   f.write(b)
 rows.append({'file':name,'url':url,'bytes':(p/name).stat().st_size,'sha256':hashlib.sha256((p/name).read_bytes()).hexdigest()})
(p/'DOWNLOAD_MANIFEST.json').write_text(json.dumps(rows,indent=2))
print('DOWNLOAD_COMPLETE',flush=True)
